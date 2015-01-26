#include "align_eval_aux.h"
#include "file_utils.h"

#include <parallel/algorithm>
#include <ext/algorithm>
#include <zlib.h>
#include <cassert>
#include <limits.h>

//tmp_fhs are open.  samlines chunk_buffer_out allocated.  line_index will be
//appended with a K-sorted chunk range
/*
1. index
2. sort
3. write sorted chunk to buffer
4. write buffer to file
5. append sorted index to main index

Returns a pair of <number bytes, number of lines> in chunk processed
*/

struct init_index_input_t
{
    char **samline;
    sam_index *line_index;
    sam_index *end;
    SAM_INDEX_TYPE itype;
    SAM_QNAME_FORMAT qfmt;
    const contig_dict *cdict;
    index_dict_t *idict;
};

void *init_index_worker(void *input)
{
    init_index_input_t *ii = static_cast<init_index_input_t *>(input);
    while (ii->line_index != ii->end)
    {
        set_sam_index(*ii->samline++, ii->itype, ii->qfmt, ii->cdict, ii->idict, ii->line_index++);
    }
    pthread_exit((void*) 0);
}



void init_index(size_t num_threads,
                sam_index **line_index_loads, // contains num_threads + 1.  (last one is the sentinel
                char **samlines,
                SAM_INDEX_TYPE itype,
                SAM_QNAME_FORMAT qfmt,
                const contig_dict *cdict,
                index_dict_t *work_dict)
{

    pthread_t * threads = new pthread_t[num_threads];
    init_index_input_t * inputs = new init_index_input_t[num_threads];

    // initialize inputs
    for (size_t t = 0; t != num_threads; ++t)
    {
        inputs[t].samline = samlines + (line_index_loads[t] - line_index_loads[0]);
        inputs[t].line_index = line_index_loads[t];
        inputs[t].end = line_index_loads[t+1];
        inputs[t].itype = itype;
        inputs[t].qfmt = qfmt;
        inputs[t].cdict = cdict;
        inputs[t].idict = &work_dict[t];
    }

    for (size_t t = 0; t != num_threads; ++t)
    {
        int rc = pthread_create(&threads[t], NULL, &init_index_worker, static_cast<void *>(&inputs[t]));
        assert(! rc);
    }

    for (size_t t = 0; t != num_threads; ++t) {
        int rc = pthread_join(threads[t], NULL);
        assert(! rc);
    }

    delete threads;
    delete inputs;
}


// initialize <num_threads> identical copies of dict, storing them in work_dict
void init_work_dicts(index_dict_t *dict, index_dict_t *work_dict, size_t num_threads)
{
    index_dict_t::iterator w;
    for (size_t t = 0; t != num_threads; ++t)
    {
        for (w = (*dict).begin(); w != (*dict).end(); ++w)
        {
            char *key = new char[strlen((*w).first) + 1];
            strcpy(key, (*w).first);
            work_dict[t].insert(std::make_pair(key, (*w).second));
        }
    }    
}


// merge each of the <num_threads> work_dict dictionaries into index_dict,
// then generate a remap array for each work chunk.
// appends to dict and initializes work_remap
void merge_and_remap_dict(index_dict_t *work_dict, size_t num_threads, 
                          index_dict_t *dict, unsigned int **work_remap)
{
    index_dict_t::iterator w;
    for (size_t t = 0; t != num_threads; ++t)
    {
        work_remap[t] = new unsigned int[work_dict[t].size()];
        for (w = work_dict[t].begin(); w != work_dict[t].end(); ++w)
        {
            if ((*dict).find((*w).first) == (*dict).end())
            {
                char *key = new char[strlen((*w).first) + 1];
                strcpy(key, (*w).first);
                (*dict).insert(std::make_pair(key, (*dict).size()));
            }
        }
    }
    for (size_t t = 0; t != num_threads; ++t)
    {
        index_dict_t::iterator ww;
        for (w = work_dict[t].begin(); w != work_dict[t].end(); ++w)
        {
            ww = (*dict).find((*w).first);
            assert(ww != (*dict).end());
            work_remap[t][(*w).second] = (*ww).second;
        }
    }
}


struct apply_remap_input_t
{
    sam_index *line_index;
    sam_index *end;
    unsigned int *remap;
    SAM_INDEX_TYPE itype;
};

void *apply_remap_worker(void *input)
{
    apply_remap_input_t *ii = static_cast<apply_remap_input_t *>(input);
    while (ii->line_index != ii->end)
    {
        sam_update_flowcell_id(ii->remap, ii->itype, ii->line_index++);
    }
    pthread_exit((void*) 0);
}

//transform a set of <num_threads> loads of index, remapping the
//flowcell ids according to work_remap
void apply_remap(size_t num_threads, sam_index **line_index_loads, 
                 SAM_INDEX_TYPE itype, 
                 SAM_QNAME_FORMAT qfmt,
                 unsigned int **work_remap)
{
    if (qfmt != SAM_ILLUMINA &&
        qfmt != SAM_CASAVA18)
        return;

    pthread_t * threads = new pthread_t[num_threads];
    apply_remap_input_t * inputs = new apply_remap_input_t[num_threads];

    // initialize inputs
    for (size_t t = 0; t != num_threads; ++t)
    {
        inputs[t].line_index = line_index_loads[t];
        inputs[t].end = line_index_loads[t+1];
        inputs[t].itype = itype;
        inputs[t].remap = work_remap[t];
    }

    for (size_t t = 0; t != num_threads; ++t)
    {
        int rc = pthread_create(&threads[t], NULL, &apply_remap_worker, static_cast<void *>(&inputs[t]));
        assert(! rc);
    }

    for (size_t t = 0; t != num_threads; ++t) {
        int rc = pthread_join(threads[t], NULL);
        assert(! rc);
    }

    delete threads;
    delete inputs;
}


// de-allocate work_dict keys
void clean_work_dicts(index_dict_t *work_dict, size_t num_threads)
{
    index_dict_t::iterator w;
    for (size_t t = 0; t != num_threads; ++t)
    {
        for (w = work_dict[t].begin(); w != work_dict[t].end(); ++w)
        {
            delete (*w).first;
            work_dict[t].erase(w);
        }
    }    
}

sam_index **create_load_ranges(sam_index *line_index, size_t num_threads, size_t num_lines)
{
    sam_index **loads = new sam_index*[num_threads + 1];
    for (size_t t = 0; t != num_threads + 1; ++t)
    {
        loads[t] = line_index + (num_lines * t) / num_threads;
    }
    return loads;
}

// initialize line_index_chunk, and augment flowcell_dict
void samlines_to_index(size_t num_threads, 
                       char **samlines,
                       size_t num_lines,
                       SAM_INDEX_TYPE itype,
                       SAM_QNAME_FORMAT qfmt,
                       const contig_dict *cdict,
                       sam_index *line_index_chunk,
                       index_dict_t *flowcell_dict)
{
    sam_index **line_index_loads = create_load_ranges(line_index_chunk, num_threads, num_lines);

    index_dict_t *work_dict = new index_dict_t[num_threads];
    init_work_dicts(flowcell_dict, work_dict, num_threads);

    // Let each thread traverse the input, updating its local copy
    // with new entries as it sees new flowcells.
    init_index(num_threads, line_index_loads, samlines, itype, qfmt, cdict, work_dict);

    /* Integrate each of the T flowcell id dictionaries into the
       working copy. */
    unsigned int **work_remap = new unsigned int*[num_threads];
    
    // Compare each of the T flowcell id dictionaries to the new
    // working copy, creating a remap array
    merge_and_remap_dict(work_dict, num_threads, flowcell_dict, work_remap);
    
    clean_work_dicts(work_dict, num_threads);
    
    // Apply the remap array to the same workload, updating the mappings.
    apply_remap(num_threads, line_index_loads, itype, qfmt, work_remap);
    
    // clean up
    delete [] work_dict;
    for (size_t t = 0; t != num_threads; ++t)
        delete work_remap[t];
    
    delete work_remap;
}



std::pair<size_t, size_t> 
process_chunk(std::vector<char *> & samlines,
              char * chunk_buffer_in,
              char * chunk_buffer_out,
              size_t num_threads,
              SAM_INDEX_TYPE itype,
              SAM_QNAME_FORMAT qfmt,
              const contig_dict *cdict,
              index_dict_t * flowcell_dict, // this will accumulate each time process_chunk is called
              std::vector<sam_index> * line_index)
{
    size_t S = samlines.size();

    if (S == 0)
    {
        return std::pair<size_t, size_t>(0, 0);
    }

    sam_index * line_index_chunk = new sam_index[S];

    samlines_to_index(num_threads, samlines.data(), S, itype, qfmt, 
                      cdict, line_index_chunk, flowcell_dict);

    set_start_offsets(line_index_chunk, line_index_chunk + S, 0);

    __gnu_parallel::sort(line_index_chunk, line_index_chunk + S, samidx_less_key);

    

    char * write_pointer = chunk_buffer_out;
    sam_index * si;
    for (si = line_index_chunk; si != line_index_chunk + S; ++si)
    {
        memcpy(write_pointer,
               chunk_buffer_in + (*si).start_offset, (*si).line_length);

        write_pointer[(*si).line_length - 1] = '\n';
        write_pointer += (*si).line_length;
    }

    size_t nbytes_written = std::distance(chunk_buffer_out, write_pointer);

    (*line_index).insert((*line_index).end(), line_index_chunk, line_index_chunk + S);
    delete [] line_index_chunk;

    return std::make_pair(nbytes_written, S);
}


// find 'num_chunks' quantiles based on a key ordering. assume
// line_index is ordered by (K,.). construct sam_index sentinels to be
// used to query each sub-range during catenate_subchunks.  since the
// constructed sentinels will be real entries in line_index, exactly
// one subchunk will contain the sentinel, and the 'lower_bound' call
// in catenate_subchunks will not 
void
get_key_quantiles(std::vector<sam_index> const& line_index,
                  size_t num_chunks,
                  size_t * key_quantile_sizes,
                  size_t * key_quantile_nlines,
                  std::vector<sam_index> * key_quantile_sentinels)
{
    sam_index const** line_index_ptr = new sam_index const*[line_index.size()];
    sam_index const** p, ** end = line_index_ptr + line_index.size();
    std::vector<sam_index>::const_iterator lit;
    for (p = line_index_ptr, lit = line_index.begin(); p != end; ++p, ++lit)
    {
        *p = &(*lit);
    }

    sam_index const*** key_quantiles =
        get_quantiles(line_index_ptr, line_index.size(), samidx_less_key, num_chunks);
    
    std::fill(key_quantile_sizes, key_quantile_sizes + num_chunks, 0);

    sam_index const **kit, **kit_end;

    for (size_t k = 0; k != num_chunks; ++k)
    {
        kit = key_quantiles[k];
        kit_end = key_quantiles[k+1];
        key_quantile_nlines[k] = kit_end - kit;

        sam_index sentinel = (kit_end == end) ? samidx_make_max() : **kit_end;
        (*key_quantile_sentinels).push_back(sentinel);

        size_t * kq_size = key_quantile_sizes + k;
        for (kit = key_quantiles[k]; *kit != *kit_end; ++kit)
        {
            (*kq_size) += (**kit).line_length;
        }
    }
    delete line_index_ptr;
    delete key_quantiles;
}

/*
// the original
void
get_key_quantiles(std::vector<sam_index> const& line_index,
                  size_t num_chunks,
                  size_t * key_quantile_sizes,
                  size_t * key_quantile_nlines,
                  std::vector<sam_index> * key_quantile_sentinels)
{
    std::vector<sam_index> line_index_copy(line_index);

    std::vector<INDEX_ITER> key_quantiles =
        get_quantiles(&line_index_copy, samidx_less_key, num_chunks);
    
    std::fill(key_quantile_sizes, key_quantile_sizes + num_chunks, 0);

    for (size_t k = 0; k != num_chunks; ++k)
    {
        INDEX_ITER kit, kit_end = key_quantiles[k+1];
        kit = key_quantiles[k];
        key_quantile_nlines[k] = std::distance(kit, kit_end);

        sam_index sentinel = (kit_end == line_index_copy.end())
            ? samidx_make_max()
            : *kit_end;

        (*key_quantile_sentinels).push_back(sentinel);

        size_t * kq_size = key_quantile_sizes + k;
        for (kit = key_quantiles[k]; kit != kit_end; ++kit)
        {
            (*kq_size) += (*kit).line_length;
        }
    }
}
*/


// structure to pass to for_each loop in catenate_subchunks.  defines
// the information needed to gzread a subchunk and write it to the
// proper section of the buffer, in parallel.
struct gzread_target
{
    gzread_target();
    gzread_target(gzFile, char *, size_t);
    gzFile source_file;
    char * dest_buffer;
    size_t bytes_to_read;
};

gzread_target::gzread_target() : source_file(NULL), dest_buffer(NULL), bytes_to_read(0) { }
gzread_target::gzread_target(gzFile s, char *d, size_t b) : source_file(s), dest_buffer(d), bytes_to_read(b) { }

struct gzread_and_copy
{
    gzread_and_copy();
    void operator()(gzread_target & item);
};

gzread_and_copy::gzread_and_copy() { }

// since gzread cannot read more than UINT_MAX uncompressed bytes,
// we need to loop until it does so.
void gzread_and_copy::operator()(gzread_target & item)
{

    size_t bytes_remaining = item.bytes_to_read;
    unsigned int bytes_to_read;
    unsigned int bytes_read;
    char * write_pointer = item.dest_buffer;

    while (bytes_remaining > 0)
    {
        bytes_to_read = bytes_remaining < UINT_MAX ? bytes_remaining : UINT_MAX;
        bytes_read = gzread(item.source_file, write_pointer, bytes_to_read);
        if ((size_t)bytes_read != bytes_to_read) {
            fprintf(stderr, "Error, tried to read %i bytes, could only read %i bytes, from tmp file\n",
                    bytes_read, bytes_to_read);
            exit(23);
        }
        bytes_remaining -= bytes_read;
        write_pointer += bytes_read;
    }
}




// gathers next sub-chunks from each temp file, starting at
// 'prev_chunk_starts', and using 'query_key_quantile' to find
// the end of each subchunk.
// postconditions: 
// catenated_index populated and ordered as (O, K)
// chunk_buffer populated to match catenated_index
// prev_chunk_starts updated to next set of chunk starts for next call
// returns number of bytes written
size_t catenate_subchunks(sam_index const& query_key_quantile,
                          // std::vector<FILE *> const& tmp_fhs,
                          std::vector<gzFile> const& tmp_fhs,
                          std::vector<INDEX_ITER> * prev_chunk_starts,
                          std::vector<INDEX_ITER> const& chunk_ends, /* ends of each chunk ordered by (O, K) */
                          char ** chunk_buffer,
                          std::vector<sam_index> * catenated_index,
                          size_t ** subrange_sizes)
{
    size_t num_chunks = tmp_fhs.size();

    std::pair<INDEX_ITER, INDEX_ITER> * subrange_iters =
        new std::pair<INDEX_ITER, INDEX_ITER>[num_chunks];

    char * write_pointer = (*chunk_buffer);
    size_t total_subrange_size = 0;

    INDEX_ITER beg, end, sub_k;
    size_t S;
    std::vector<gzread_target> subrange_bytes(num_chunks);

    // initialize subrange_iters, subrange_sizes and total_subrange_size
    for (size_t o = 0; o != num_chunks; ++o)
    {
        beg = (*prev_chunk_starts)[o];
        end = chunk_ends[o+1];

        assert(std::is_sorted(beg, end, samidx_less_key));

        sub_k = std::lower_bound(beg, end, query_key_quantile, samidx_less_key);

        // fprintf(stdout, "key %Zu, offset %Zu, num_lines: %Zu\n",
        //         k, o, part);

        S = 0;
        for (INDEX_ITER kit = beg; kit != sub_k; ++kit)
        {
            S += (*kit).line_length;
        }

        // the information we need to perform parallel gzreads below in the for_each loop
        subrange_bytes[o] = gzread_target(tmp_fhs[o], write_pointer, S);

        // gzread(tmp_fhs[o], write_pointer, S);

        write_pointer += S;

        // last_merged_size = merged.size();
        subrange_iters[o] = std::make_pair(beg, sub_k);
        (*subrange_sizes)[o] = std::distance(beg, sub_k);
        total_subrange_size += (*subrange_sizes)[o];

        // fprintf(stdout, "key %Zu, offset %Zu, added: %Zu, num_lines: %Zu\n",
        //         k, o, part, merged.size());

        (*prev_chunk_starts)[o] = sub_k;
    }

    //initialize data pointed to by chunk_buffer by performing
    //multiple gzread operations at the same time.  Unfortunately,
    //for_each doesn't allow passing the iterator itself.
    //so we use a dummy 'index'
    __gnu_parallel::for_each(subrange_bytes.begin(), subrange_bytes.end(),
                             gzread_and_copy());

    (*catenated_index).reserve(total_subrange_size);
    (*catenated_index).resize(0);

    for (size_t o = 0; o != num_chunks; ++o)
    {
        catenated_index->insert(catenated_index->end(), 
                                subrange_iters[o].first, 
                                subrange_iters[o].second);
    }

    delete [] subrange_iters;
    return total_subrange_size;
}


// merge an (O, K) ordered index.
// precondition: ok_index is an (O, K) ordered index in 'num_ranges' pieces,
// of size 'subrange_sizes'
// postcondition: k_index is the (K) ordered index
// warning: destroys the original ok_index

void merge_ok_index(size_t const* subrange_sizes,
                    size_t num_ranges,
                    std::vector<sam_index> ** pre_merge_index,
                    std::vector<sam_index> ** post_merge_index)
{

    size_t * sz_tmp = new size_t[num_ranges + 1];
    std::copy(subrange_sizes, subrange_sizes + num_ranges, sz_tmp);
    sz_tmp[num_ranges] = 0;

    size_t total_size = std::accumulate(subrange_sizes, subrange_sizes + num_ranges, 0);

    (*post_merge_index)->reserve(total_size);
    (*post_merge_index)->resize(total_size);

    //here's the tricky part.
    //iteratively merge consecutive ranges until there is only one.
    //invariant: subrange_sizes valid, num_ranges valid, k_index valid
    size_t num_merge_ranges = num_ranges + (num_ranges % 2); // square it off

    std::vector<sam_index>::iterator pre, post; //before and after a pairwise range merge.
    
    //precondition: pre_merge_index is (O, K) sorted, in sz_tmp partitions
    //postcondition: post_merge_index is (O, K) sorted in half as many partitions as pre_merge_index
    
    //this makes it appear we had done at least one iteration.
    std::swap((*post_merge_index), (*pre_merge_index));
    while (num_merge_ranges > 1)
    {
        std::swap((*post_merge_index), (*pre_merge_index));
        pre = (*pre_merge_index)->begin();
        post = (*post_merge_index)->begin();

        for (size_t o = 0; o != num_merge_ranges; o += 2)
        {
            __gnu_parallel::merge(pre, pre + sz_tmp[o],
                                  pre + sz_tmp[o],
                                  pre + sz_tmp[o] + sz_tmp[o + 1],
                                  post, samidx_less_key);

            pre += sz_tmp[o] + sz_tmp[o + 1];
            post += sz_tmp[o] + sz_tmp[o + 1];

            sz_tmp[o / 2] = sz_tmp[o] + sz_tmp[o + 1];
        }
        num_merge_ranges /= 2;

        if (num_merge_ranges > 1 && (num_merge_ranges % 2) != 0)
        {
            num_merge_ranges += (num_merge_ranges % 2);
            sz_tmp[num_merge_ranges - 1] = 0;
        }
    }
    delete [] sz_tmp;

}


void check_unique_index(std::vector<sam_index> const& line_index)
{
    sam_index max_index = samidx_make_max();
    sam_index const* prior_index = &max_index;

    bool qualified;
    std::vector<sam_index>::const_iterator li;
    for (li = line_index.begin(); li != line_index.end(); ++li)
    {
        qualified = 
            samidx_less_key(*prior_index, (*li)) 
            || prior_index == &max_index;

        if (! qualified)
        {
            fprintf(stderr, "Error: check_index: non-ascending or duplicate line indices.\n");
            exit(1);
        }
        prior_index = &(*li);
    }
    return;
}



//write blocks from 'unordered' to 'ordered' according to 'ordering'
//index.  return number of bytes written
size_t write_new_ordering(char const* unordered, 
                          std::vector<sam_index> * ordering,
                          char * ordered)
{
    char * write_pointer = ordered;
    
    std::vector<sam_index>::iterator mit;
    size_t S = 0;

    for (mit = (*ordering).begin(); mit != (*ordering).end(); ++mit)
    {
        memcpy(write_pointer, unordered + (*mit).start_offset,
               (*mit).line_length);
        
        write_pointer += (*mit).line_length;
        (*mit).start_offset = S;
        S += (*mit).line_length;
    }
    return std::distance(ordered, write_pointer);
}


size_t set_start_offsets(sam_index *beg, sam_index *end, size_t off)
{
    while (beg != end)
    {
        beg->start_offset = off;
        off += beg->line_length;
        ++beg;
    }
    return off;
}
                         


// Preconditions: 
// ok_index: sorted by (O, K).  start_offsets globally set
// tmp_fhs: opened and at beginning. hold K-sorted chunks consistent with ok_index
// offset_quantiles define positions in ok_index corresponding to tmp_fhs chunks
// Postconditions:
// out_dat_fh has final K-sorted records
// out_ind_fh has final K-sorted index entries

void write_final_merge(std::vector<sam_index> const& ok_index,
                       std::vector<INDEX_ITER> const& offset_quantiles,
                       std::vector<gzFile> const& tmp_fhs,
                       bool do_check_unique_index,
                       FILE * out_dat_fh,
                       FILE * out_ind_fh)
{
    size_t num_chunks = tmp_fhs.size();

    std::vector<sam_index> key_quantile_sentinels;
    size_t * key_quantile_sizes = new size_t[num_chunks];
    size_t * key_quantile_nlines = new size_t[num_chunks];

    // 2. Find N quantile key values in index.  does not change
    // order of line_index.
    get_key_quantiles(ok_index, num_chunks, key_quantile_sizes, 
                      key_quantile_nlines,
                      & key_quantile_sentinels);

    size_t kq_size_max = 
        *std::max_element(key_quantile_sizes, key_quantile_sizes + num_chunks);

    // 4. Allocate chunk of memory to hold biggest key_quantile
    char * chunk_buffer_in = new char[kq_size_max];
    char * chunk_buffer_out = new char[kq_size_max];

    size_t * subrange_sizes = new size_t[num_chunks];

    std::vector<INDEX_ITER> prev_sub_k(offset_quantiles.size() - 1);
    std::copy(offset_quantiles.begin(), offset_quantiles.end() - 1, prev_sub_k.begin());

    fprintf(stderr, "Merging chunks [0-%zu]:", num_chunks - 1);
    fflush(stderr);

    size_t chunk_offset = 0;

    for (size_t k = 0; k != num_chunks; ++k)
    {

        std::vector<sam_index> buf1, buf2;
        std::vector<sam_index> * ok_index_ptr = & buf1;
        std::vector<sam_index> * k_index_ptr = & buf2;

        size_t lines_written = 
            catenate_subchunks(key_quantile_sentinels[k],
                               tmp_fhs,
                               & prev_sub_k,
                               offset_quantiles,
                               & chunk_buffer_in,
                               ok_index_ptr,
                               & subrange_sizes);

        // now initialize start offsets
        set_start_offsets(&(*ok_index_ptr)[0], &(*ok_index_ptr)[0] + (*ok_index_ptr).size(), 0);
        
        if (lines_written != key_quantile_nlines[k])
        {
            fprintf(stderr, "lines_written: %zu, key_quantile_nlines[%zu]: %zu\n",
                    lines_written, k, key_quantile_nlines[k]);
            assert(false);
        }

        size_t bytes_written;
        //assert(bytes_written == key_quantile_sizes[k]);

        merge_ok_index(subrange_sizes, num_chunks, & ok_index_ptr, & k_index_ptr);

        if (do_check_unique_index)
        {
            check_unique_index(* k_index_ptr);
        }

        bytes_written = 
            write_new_ordering(chunk_buffer_in, k_index_ptr, chunk_buffer_out);

        //assert(bytes_written == key_quantile_sizes[k]);
        fwrite(chunk_buffer_out, 1, key_quantile_sizes[k], out_dat_fh);
        fflush(out_dat_fh);

        if (out_ind_fh != NULL)
        {
            for (std::vector<sam_index>::const_iterator i = (*k_index_ptr).begin();
                 i != (*k_index_ptr).end(); ++i)
            {
                fprintf(out_ind_fh, "%zu:%zu\t%zu\t%zu\n",
                        (*i).key.raw[0], 
                        (*i).key.raw[1], 
                        (*i).start_offset + chunk_offset, (*i).line_length);
            }
        }
        fflush(out_ind_fh);

        chunk_offset += bytes_written;

        fprintf(stderr, " %zu", k);
        fflush(stderr);
    }

    delete [] chunk_buffer_in;
    delete [] chunk_buffer_out;
    delete [] key_quantile_sizes;
    delete [] key_quantile_nlines;
    delete [] subrange_sizes;

}


/*
There is a gotcha in using nth_element.  The postcondition states:

"There exists no iterator i in the range [first, nth) such that *nth <
*i, and there exists no iterator j in the range [nth + 1, last) such
that *j < *nth."

But, if [first, end) happens to have two identical elements, which
would occupy positions nth and (n+1)th, then it is unspecified which
of the two will be assigned to nth.

 */
/*
  std::vector<INDEX_ITER> 
  get_quantiles(std::vector<sam_index> * line_index,
  bool (less_fcn)(sam_index const&, sam_index const&),
  size_t num_chunks)
  {
  
  size_t lines_per_chunk = (*line_index).size() / num_chunks;
  std::vector<INDEX_ITER> quantiles(num_chunks);
  
  INDEX_ITER iter = (*line_index).begin();
  INDEX_ITER end = (*line_index).end();
  
  quantiles[0] = iter;
  for (size_t n = 0; n != num_chunks - 1; ++n)
  {
  std::advance(iter, lines_per_chunk);
  std::nth_element(quantiles[n], iter, end, less_fcn);
  //__gnu_parallel::nth_element(quantiles[n], iter, end, less_fcn);
  quantiles[n + 1] = iter;
  }
  quantiles.push_back(end);
  return quantiles;
  }
*/

struct less_wrap
{
    bool (*less)(sam_index const&, sam_index const&);
    bool operator()(sam_index const*a, sam_index const*b)
    {
        return this->less(*a, *b);
    }
};


// std::vector<INDEX_ITER> 
// returns an array of pointers into line_index, after partially sorting line_index
// j
sam_index const***
get_quantiles(sam_index const** line_index,
              size_t index_size,
              bool (less_fcn)(sam_index const&, sam_index const&),
              size_t num_chunks)
{
    
    less_wrap less_ptr;
    less_ptr.less = less_fcn;

    size_t lines_per_chunk = index_size / num_chunks;
    sam_index const*** quantiles = new sam_index const**[num_chunks + 1];
    
    sam_index const**iter = line_index;
    sam_index const**end = line_index + index_size;

    quantiles[0] = iter;
    for (size_t n = 0; n != num_chunks - 1; ++n)
    {
        iter += lines_per_chunk;
        // std::nth_element(quantiles[n], iter, end, less_ptr);
        __gnu_parallel::nth_element(quantiles[n], iter, end, less_ptr);
        quantiles[n + 1] = iter;
    }
    quantiles[num_chunks] = end;
    return quantiles;
}
