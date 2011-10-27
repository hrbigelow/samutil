#include <vector>
#include <cstdio>

#include <omp.h>
#include <parallel/algorithm>

#include "align_eval_aux.h"
#include "file_utils.h"
#include "sam_score_aux.h"
#include "sam_aux.h"

int sam_rejoin_usage(size_t mdef, size_t qdef,
                     size_t rdef, size_t Sdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil rejoin [OPTIONS] seq.fqi seq.fqd input.rsam output.sam [output.unmapped.fq]\n\n"
            "Options:\n\n"
            "-m     INT     number bytes of memory to use [%zu]\n"
            "-t     INT     number of threads to use [1]\n"
            "-C     STRING  work in the directory named here [.]\n"
            "-f     STRING  symbolic record filtering string [<no filter>]\n"
            "-q     INT     only output records with this or greater mapping quality [%Zu]\n"
            "-A     STRING  only output records from these alignment spaces [<all spaces>]\n"
            "-r     INT     only output records with this or smaller stratum rank [%zu]\n"
            "-S     INT     only output records with this or smaller stratum size [%zu]\n"
            "\n\n"
            "NOTES:\n"
            "       seq.{fqi,fqd} are produced from 'samutil seqindex' on the fastq\n"
            "       data used to generate input.rsam\n"
            "\n"
            "       The symbolic filtering string has format [Mm][Pp][Qq][Dd].\n"
            "       M (m) means output only mapped (unmapped) fragments (all_fragments_mapped flag)\n"
            "       P (p) means output only primary (non-primary) alignments (inverse of 'alignment_not_primary' flag)\n"
            "       Q (q) means output only passing (failing) quality check (inverse of 'failed_quality_check' flag)\n"
            "       D (d) means output only non-duplicate (duplicate) entries (inverse of 'pcr_or_optical_duplicate' flag)\n"
            "       Missing letters denote that no filtering occurs on that field and both categories are output\n"
            "\n\n"
            "       If present on command line, <output.unmapped.fq> is an intercalated fastq formatted\n"
            "       file of all unmapped records. Entries do not have trailing /1, /2 etc, and entries\n"
            "       are printed in intercalated fashion.  Use 'deal_fastq' to unmix them\n"
            ""
            , mdef, qdef, rdef, Sdef);
    return 1;
}


struct less_fragment_id
{
    bool operator()(SamLine const* a, SamLine const* b)
    {
        return a->fragment_id < b->fragment_id;
    }
};


struct less_index
{
    bool operator()(LineIndex const& a, LineIndex const& b)
    {
        return a.index < b.index;
    }
};

std::vector<SamLine *>::iterator
find_rsam_upper_bound(std::vector<SamLine *>::iterator beg,
                      std::vector<SamLine *>::iterator end,
                      size_t index);

int main_sam_rejoin(int argc, char ** argv)
{
    char c;

    size_t mdef = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    char const* fdef = "";
    size_t qdef = 0;
    size_t rdef = 1000;
    size_t Sdef = 1000;

    size_t max_mem = mdef;
    char const* tag_filter = fdef;

    size_t min_mapping_quality = qdef;

    char const* alignment_space_filter = NULL;

    size_t max_stratum_rank = rdef;

    size_t max_stratum_size = Sdef;

    size_t num_threads = 1;
    char const* working_dir = ".";


    while ((c = getopt(argc, argv, "m:t:C:f:q:A:r:S:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        case 'f': tag_filter = optarg; break;
        case 'q': min_mapping_quality = static_cast<size_t>(atof(optarg)); break;
        case 'A': alignment_space_filter = optarg; break;
        case 'r': max_stratum_rank = static_cast<size_t>(atof(optarg)); break;
        case 'S': max_stratum_size = static_cast<size_t>(atof(optarg)); break;
        default: return sam_rejoin_usage(mdef, qdef, rdef, Sdef); break;
        }
    }

    if ((argc != optind + 4) && (argc != optind + 5))
    {
        return sam_rejoin_usage(mdef, qdef, rdef, Sdef);
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    char * seq_index_file = argv[optind];
    char * seq_data_file = argv[optind + 1];
    char * input_rsam_file = argv[optind + 2];
    char * output_sam_file = argv[optind + 3];
    char * output_fastq_file = argc == optind + 5 ? argv[optind + 4] : NULL;


    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * seq_index_fh = open_or_die(seq_index_file, "r", "Input seq index file");
    FILE * seq_data_fh = open_or_die(seq_data_file, "r", "Input seq data file");
    FILE * input_rsam_fh = open_or_die(input_rsam_file, "r", "Input rsam-formatted file");
    FILE * output_sam_fh = open_or_die(output_sam_file, "w", "Output sam-formatted file");
    FILE * output_fastq_fh = open_if_present(output_fastq_file, "w");

    // assume we have rSAM format.  By definition, it requires
    // 'SAM_NUMERIC' and empty expected layout.
    SamLine::SetGlobalFlags(SAM_NUMERIC, "", "", 0, false);

    size_t index, start_offset;
    int line_length;

    std::vector<LineIndex> seq_index;
    seq_index.reserve(1000000); // one extra for a sentinel value
    while (! feof(seq_index_fh))
    {
        int nfields = 
            fscanf(seq_index_fh, "%zu\t%zu\t%i\n", &index, &start_offset, &line_length);

        if (nfields != 3)
        {
            fprintf(stderr, "Error: line %zu does not have three integers "
                    "(index, start_offset, line_length)\n", seq_index.size());
        }
        seq_index.push_back(LineIndex(index, start_offset, line_length));
    }
    size_t total_seq_size = (*seq_index.rbegin()).start_offset
        + (*seq_index.rbegin()).line_length;

    seq_index.push_back(LineIndex(SIZE_MAX, total_seq_size, 0));

    fclose(seq_index_fh);

    // aim to parse record-matching amounts of rsam and seq in
    // general, each rsam record is about 1/4 the size of the data
    // line, but there may be about 4 lines, so perhaps they should be
    // the same size
    size_t rsam_chunk_size = max_mem / 2;
    size_t seq_chunk_size = max_mem / 2;
    
    char * rsam_buffer = new char[rsam_chunk_size + 1];
    char * seq_buffer = new char[seq_chunk_size + 1];
    char * read_pointer = rsam_buffer;

    std::vector<SamLine *> rsam_records;

    // beg and end define the range within rsam_records that is to be
    // populated 
    std::vector<SamLine *>::iterator beg = rsam_records.end();
    std::vector<SamLine *>::iterator end = rsam_records.end();

    std::vector<char *> sam_strings;
    std::vector<char const*> fqd_strings;

    std::vector<char *> sam_lines;
    size_t nbytes_read, nbytes_unused = 0;
    size_t n_databytes_read = 0;
    size_t data_buffer_offset = 0;
    char * last_fragment;
    
    // iterator bounding last loaded record
    INDEX_ITER last_index_with_data = seq_index.end(); // symbolizes 'uninitialized'
    INDEX_ITER index_with_data_bound = seq_index.end();

    char * header_buf = ReadAllocSAMHeader(input_rsam_fh);
    size_t header_length = strlen(header_buf);
    fwrite(header_buf, 1, header_length, output_sam_fh);
    delete [] header_buf;

    SamFilter sam_filter(tag_filter, 
                         min_mapping_quality,
                         max_stratum_rank,
                         max_stratum_size,
                         alignment_space_filter);

    while (! rsam_records.empty() || ! feof(input_rsam_fh))
    {
        // 1) All records in rec have been processed and purged. Load
        // new records

        // precondition: beg == end == rec.end()
        // postcondition: beg = rec.begin(),end = upper_bound(fqd)
        if (beg == rsam_records.end() && end == rsam_records.end())
        {
            nbytes_read = fread(read_pointer, 1, rsam_chunk_size - nbytes_unused, input_rsam_fh);
            read_pointer[nbytes_read] = '\0';
            sam_lines = FileUtils::find_complete_lines_nullify(rsam_buffer, & last_fragment);

            rsam_records.resize(sam_lines.size());
            sam_strings.resize(sam_lines.size());

            __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                      rsam_records.begin(), parse_sam_unary());

            beg = rsam_records.begin();
            end = (last_index_with_data == seq_index.end())
                ? rsam_records.begin()
                : find_rsam_upper_bound(rsam_records.begin(), rsam_records.end(), (*last_index_with_data).index);

            nbytes_unused = strlen(last_fragment);
            memmove(rsam_buffer, last_fragment, nbytes_unused);
            read_pointer = rsam_buffer + nbytes_unused;
            
        }
        
        // 2) [beg, end) is an empty range, but not at rec.end().
        // means previous attempt to convert ran out of fqd data, so
        // we must load more.

        // precondition: beg == end && end != rec.end()
        // postcondition: end = upper_bound(fqd)
        if (beg == end && end != rsam_records.end())
        {
            data_buffer_offset += n_databytes_read;

            // last_index_with_data should be the highest iterator such that
            // start_offset <= seq_chunk_size.  the query index
            // element, with offset 'seq_chunk_size +
            // data_buffer_offset' is the next data item that would be
            // loaded.  so, the range [seq_index.begin(),
            // last_index_with_data) is guaranteed to have last_index_with_data
            // should be the last real index element with loaded seq
            // data.
            assert(! seq_index.empty());
            index_with_data_bound = 
                std::lower_bound(seq_index.begin(), seq_index.end(),
                                 LineIndex(0, seq_chunk_size + data_buffer_offset, 0),
                                 &less_offset);
            
            last_index_with_data = index_with_data_bound - 1;
            
            size_t data_bytes_to_read = 
                (*last_index_with_data).start_offset
                + (*last_index_with_data).line_length
                - data_buffer_offset;

            //assert(data_bytes_to_read <= seq_chunk_size);

            // update data_buffer_offset *before* loading new data.
            n_databytes_read = fread(seq_buffer, 1, data_bytes_to_read, seq_data_fh);

            end = find_rsam_upper_bound(beg, rsam_records.end(), (*last_index_with_data).index);
        }
                
        // 3) The range [beg, end) needs to be converted to SAM strings
        // precondition: beg != end
        // convert to SAM strings
        // postcondition: beg == end
        if (beg != end)
        {
            size_t beg_elem = std::distance(rsam_records.begin(), beg);

            // serially construct a corresponding vector of index iterators
            std::vector<INDEX_ITER> index_iters(std::distance(beg, end));
            std::vector<INDEX_ITER>::iterator ii_iter;
            std::vector<SamLine *>::iterator rec_iter;

            // find the line index corresponding to the beg iterator
            INDEX_ITER lit = 
                std::lower_bound(seq_index.begin(), index_with_data_bound, 
                                 LineIndex((*beg)->fragment_id, 0, 0), 
                                 less_index());

            assert((*lit).index == (*beg)->fragment_id);

            for (ii_iter = index_iters.begin(), rec_iter = beg; 
                 ii_iter != index_iters.end(); 
                 ++ii_iter, ++rec_iter)
            {
                while ((*lit).index < (*rec_iter)->fragment_id
                       && lit != index_with_data_bound)
                {
                    ++lit;
                }
                if ((*lit).index != (*rec_iter)->fragment_id)
                {
                    fprintf(stderr, "Error: couldn't find read index %Zu in fastq index file\n",
                            (*rec_iter)->fragment_id);
                    exit(1);
                }
                (*ii_iter) = lit;
            }

            rsam_to_sam_binary rsb(seq_buffer, data_buffer_offset, & sam_filter);

            __gnu_parallel::transform(beg, end, 
                                      index_iters.begin(),
                                      sam_strings.begin() + beg_elem,
                                      rsb);

            unmapped_rsam_to_fastq usb(seq_buffer, data_buffer_offset);

            if (output_fastq_fh != NULL)
            {
                fqd_strings.resize(std::distance(beg, end));
                __gnu_parallel::transform(beg, end,
                                          index_iters.begin(),
                                          fqd_strings.begin(),
                                          usb);

                std::vector<char const*>::iterator sit;
                for (sit = fqd_strings.begin(); sit != fqd_strings.end(); ++sit)
                {
                    char const* fqd_string = *sit;
                    if (fqd_string != NULL)
                    {
                        print_fqd_as_fastq(fqd_string, output_fastq_fh);
                        // do not delete fastq string -- it's owned by the index
                    }
                }
                fflush(output_fastq_fh);
            }

            __gnu_parallel::for_each(beg, end, delete_samline());
            
            beg = end;
        }

        // 4) all rsam_records have been populated, time to purge them
        // precondition: beg == end == rec.end(), 
        // purge converted SAM strings
        // postcondition: rec.empty(), beg == rec.end(), end == rec.end()
        if (beg == end && end == rsam_records.end())
        {
            std::vector<char *>::iterator sit;
            for (sit = sam_strings.begin(); sit != sam_strings.end(); ++sit)
            {
                char * sam_string = *sit;
                if (sam_string != NULL)
                {
                    fputs(sam_string, output_sam_fh);
                    delete sam_string;
                }
            }
            fflush(output_sam_fh);

            rsam_records.clear();
            beg = rsam_records.end();
            end = rsam_records.end();
        }
        
    }

    fclose(seq_data_fh);
    fclose(input_rsam_fh);
    fclose(output_sam_fh);
    close_if_present(output_fastq_fh);

    delete [] rsam_buffer;
    delete [] seq_buffer;

    return 0;
}


// find the iterator rit in [beg, end) such that for all it in [beg, rit)
// (*it).index <= index
std::vector<SamLine *>::iterator
find_rsam_upper_bound(std::vector<SamLine *>::iterator beg,
                      std::vector<SamLine *>::iterator end,
                      size_t index)
{
    if (beg == end)
    {
        return end;
    }
    else
    {
        SamLine rsam_dummy(**beg);
        rsam_dummy.fragment_id = index;
        std::vector<SamLine *>::iterator rit =
            std::upper_bound(beg, end, &rsam_dummy, less_fragment_id());
        
        return rit;
    }
}    
