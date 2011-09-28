#include <vector>
#include <cstdio>

#include <omp.h>
#include <parallel/algorithm>

#include "align_eval_aux.h"
#include "file_utils.h"
#include "sam_score_aux.h"
#include "sam_aux.h"

int sam_rejoin_usage(size_t mdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil rejoin [OPTIONS] seq.fqi seq.fqd input.rsam output.sam\n\n"
            "Options:\n\n"
            "-m     INT     number bytes of memory to use [%zu]\n"
            "-t     INT     number of threads to use [1]\n"
            "-C     STRING  work in the directory named here [.]\n"
            "\n\n"
            "NOTES: seq.{fqi,fqd} are produced from 'samutil seqindex' on the fastq\n"
            "       data used to generate input.rsam\n\n",
            mdef);
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
find_rsam_by_id(std::vector<SamLine *>::iterator beg,
                      std::vector<SamLine *>::iterator end,
                      size_t index);

int main_sam_rejoin(int argc, char ** argv)
{
    char c;

    size_t mdef = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = mdef;
    size_t num_threads = 1;
    char const* working_dir = ".";

    while ((c = getopt(argc, argv, "m:t:C:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        default: return sam_rejoin_usage(mdef); break;
        }
    }

    if (argc != optind + 4)
    {
        return sam_rejoin_usage(mdef);
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    char * seq_index_file = argv[optind];
    char * seq_data_file = argv[optind + 1];
    char * input_rsam_file = argv[optind + 2];
    char * output_sam_file = argv[optind + 3];

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

    // assume we have rSAM format.  By definition, it requires
    // 'SAM_NUMERIC' and empty expected layout.
    SamLine::SetGlobalFlags(SAM_NUMERIC, "", "", 0);

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

    std::vector<SamLine *> rsam_records;

    // beg and end define the range within rsam_records that is to be
    // populated 
    std::vector<SamLine *>::iterator beg = rsam_records.end();
    std::vector<SamLine *>::iterator end = rsam_records.end();

    std::vector<char *> sam_strings;

    std::vector<char *> sam_lines;
    size_t nbytes_read, nbytes_unused;
    size_t n_databytes_read = 0;
    size_t data_buffer_offset = 0;
    
    // iterator bounding last loaded record
    INDEX_ITER loaded_seq_end = seq_index.begin();

    SetToFirstDataLine(& input_rsam_fh);

    size_t header_length = ftell(input_rsam_fh);
    rewind(input_rsam_fh);

    char * header_buf = new char[header_length];
    fread(header_buf, 1, header_length, input_rsam_fh);
    fwrite(header_buf, 1, header_length, output_sam_fh);
    delete [] header_buf;


    while (! feof(input_rsam_fh))
    {
        // 1) All records in rec have been processed and purged. Load
        // new records

        // precondition: beg == end == rec.end()
        // postcondition: beg = rec.begin(),end = upper_bound(fqd)
        if (beg == rsam_records.end() && end == rsam_records.end())
        {
            nbytes_read = fread(rsam_buffer, 1, rsam_chunk_size, input_rsam_fh);
            rsam_buffer[nbytes_read] = '\0';
            sam_lines = FileUtils::find_complete_lines_nullify(rsam_buffer, & nbytes_unused);
            if (nbytes_unused > 0)
            {
                //only do this if we need to. this resets a flag that affects feof()
                fseek(input_rsam_fh, -1 * static_cast<off_t>(nbytes_unused), SEEK_CUR);
            }
            
            rsam_records.resize(sam_lines.size());
            sam_strings.resize(sam_lines.size());

            __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                      rsam_records.begin(), parse_sam_unary());

            beg = rsam_records.begin();
            end = find_rsam_by_id(rsam_records.begin(), rsam_records.end(), (*loaded_seq_end).index);
        }
        
        // 2) [beg, end) is an empty range, but not at rec.end().
        // means previous attempt to convert ran out of fqd data, so
        // we must load more.

        // precondition: beg == end && end != rec.end()
        // postcondition: end = upper_bound(fqd)
        if (beg == end && end != rsam_records.end())
        {
            data_buffer_offset += n_databytes_read;

            // loaded_seq_end should be the highest iterator such that
            // start_offset <= seq_chunk_size. The range
            // [seq_index.begin(), loaded_seq_end) is guaranteed to
            // have seq data loaded for it.
            assert(! seq_index.empty());
            loaded_seq_end = 
                std::lower_bound(seq_index.begin(), seq_index.end(),
                                 LineIndex(0, seq_chunk_size + data_buffer_offset, 0),
                                 &less_offset) - 1;
            
            size_t data_bytes_to_read = 
                (*loaded_seq_end).start_offset - data_buffer_offset;

            assert(data_bytes_to_read <= seq_chunk_size);

            // update data_buffer_offset *before* loading new data.
            n_databytes_read = fread(seq_buffer, 1, data_bytes_to_read, seq_data_fh);

            end = find_rsam_by_id(beg, rsam_records.end(), (*loaded_seq_end).index);
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
                std::lower_bound(seq_index.begin(), loaded_seq_end, 
                                 LineIndex((*beg)->fragment_id, 0, 0), 
                                 less_index());

            assert((*lit).index == (*beg)->fragment_id);

            for (ii_iter = index_iters.begin(), rec_iter = beg; 
                 ii_iter != index_iters.end(); 
                 ++ii_iter, ++rec_iter)
            {
                while ((*lit).index < (*rec_iter)->fragment_id
                       && lit != loaded_seq_end)
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

            __gnu_parallel::transform(beg, end, 
                                      index_iters.begin(),
                                      sam_strings.begin() + beg_elem,
                                      rsam_to_sam_binary(seq_buffer, data_buffer_offset));

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
                fputs(sam_string, output_sam_fh);
                delete sam_string;
            }
            rsam_records.clear();
            beg = rsam_records.end();
            end = rsam_records.end();
        }
        
    }

    fclose(seq_data_fh);
    fclose(input_rsam_fh);
    fclose(output_sam_fh);

    delete [] rsam_buffer;
    delete [] seq_buffer;

    return 0;
}


// finds the rsam iterator corresponding to index or 'end' if not
// exists.  This is because index comes from a superset of all indices
std::vector<SamLine *>::iterator
find_rsam_by_id(std::vector<SamLine *>::iterator beg,
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
            std::lower_bound(beg, end, &rsam_dummy, less_fragment_id());

        assert(rit == end || (*rit)->fragment_id == index);

        return rit;
    }
}    
