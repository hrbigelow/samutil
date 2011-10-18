#include <vector>
#include <utility>
#include <string>
#include <cstdio>

#include <omp.h>

#include "dep/tools.h"
#include "align_eval_aux.h"
#include "file_utils.h"

#include <zlib.h>

int sam_seqindex_usage(size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil seqindex [OPTIONS] read1.fastq[.gz] [read2.fastq[.gz] ...] -- reads_out.fqi reads_out.fqd\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-C  STRING    work in the directory named here [.]\n"
            "-T  STRING    path/and/prefix of temp files [name of output fqd data file]\n\n",
            mdef);

    return 1;
}


int main_sam_seqindex(int argc, char ** argv)
{

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;
    char const* working_dir = ".";

    char const* tmp_file_prefix = NULL;

    std::vector<std::string> fastq_files;

    char c;
    while ((c = getopt(argc, argv, "-m:t:C:T:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        case 'T': tmp_file_prefix = optarg; break;
        case 1: fastq_files.push_back(std::string(optarg)); break;
        default: return sam_seqindex_usage(max_mem_def); break;
        }
    }

    if (optind + 2 != argc || fastq_files.empty())
    {
        return sam_seqindex_usage(max_mem_def);
    }
    
    char const* out_index_file = argv[optind];
    char const* out_data_file = argv[optind + 1];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    if (tmp_file_prefix == NULL)
    {
        tmp_file_prefix = out_data_file;
    }

    FILE * out_index_fh = open_if_present(out_index_file, "w");
    FILE * out_data_fh = open_if_present(out_data_file, "w");

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);


    size_t nfq = fastq_files.size();
    gzFile * fastq_fhs = new gzFile[nfq];
    size_t chunk_size = max_mem / nfq / 3;
    char ** fq_chunks = new char *[nfq];
    char * fq_dat_in = new char[chunk_size * nfq];
    char * fq_dat_out = new char[chunk_size * nfq];
    
    for (size_t f = 0; f != nfq; ++f)
    {
        fastq_fhs[f] = gzopen(fastq_files[f].c_str(), "r");
        if (fastq_fhs[f] == NULL)
        {
            fprintf(stderr, "Error: couldn't open fastq file %s for reading\n",
                    fastq_files[f].c_str());
            exit(1);
        }
        fq_chunks[f] = new char[chunk_size + 1];
    }

    std::vector<std::vector<char *> > line_starts(nfq);
    std::vector<char *>::const_iterator * sit = new std::vector<char *>::const_iterator[nfq];


    std::vector<LineIndex> line_index;

    std::pair<size_t, size_t> chunk_stats;

    char * tmp_file_template = new char[strlen(tmp_file_prefix) + 8];
    strcpy(tmp_file_template, tmp_file_prefix);
    strcat(tmp_file_template, ".XXXXXX");

    std::vector<char *> tmp_files;
    std::vector<FILE *> tmp_fhs;

    // main initializion
    size_t num_complete_records = 0;
    for (size_t f = 0; f != nfq; ++f)
    {
        fq_chunks[f][0] = '\0';
        line_starts[f].push_back(fq_chunks[f]);
        sit[f] = line_starts[f].begin();
    }

    //this is hard-coded because we only ever want to sort by read id.
    SamOrder sam_order(SAM_RID, "FRAGMENT");
    sam_order.InitFromFastqFile(fastq_files[0].c_str());

    std::vector<size_t> chunk_num_lines;

    size_t nbytes_read;
    size_t nbytes_unused;
    char * last_fragment;

    do
    {
        num_complete_records = UINT64_MAX;

        for (size_t f = 0; f != nfq; ++f)
        {
            if (! gzeof(fastq_fhs[f]))
            {
                nbytes_read = gzread(fastq_fhs[f], fq_chunks[f], chunk_size);
            }
            else
            {
                nbytes_read = 0;
            }
            fq_chunks[f][nbytes_read] = '\0';
            line_starts[f] = FileUtils::find_complete_lines(fq_chunks[f], & last_fragment);
            nbytes_unused = strlen(last_fragment);

            if (nbytes_unused > 0)
            {
                gzseek(fastq_fhs[f], - nbytes_unused, SEEK_CUR);

                //setting this new NULL byte is needed later to
                //recover orphaned bytes
                fq_chunks[f][nbytes_read - nbytes_unused] = '\0';
            }
            assert(line_starts[f].empty() || line_starts[f][0][0] == '@');

            //will restore unused bytes plus orphaned bytes later

            size_t this_num_complete_records = line_starts[f].size() / 4;

            num_complete_records = std::min(this_num_complete_records, num_complete_records);
        }

        // now, restore all partial records (< 4 lines) or
        // non-matching records (not present in all fq files)
        size_t num_complete_lines = num_complete_records * 4;
        for (size_t f = 0; f != nfq; ++f)
        {
            std::vector<char *>::iterator orphan = 
                line_starts[f].begin() + num_complete_lines;

            if (orphan != line_starts[f].end())
            {
                // this call to strlen should hit the null byte set 
                off_t orphaned_bytes = strlen((*orphan));
                assert(orphaned_bytes > 0);

                gzseek(fastq_fhs[f], - orphaned_bytes, SEEK_CUR);
                line_starts[f].resize(num_complete_lines);
            }
        }        

        std::vector<char *>::iterator ls_iter;
        for (size_t f = 0; f != nfq; ++f)
        {
            sit[f] = line_starts[f].begin();

            //replace all newlines with nulls, up until the unused portion
            char * line_end = fq_chunks[f];
            for (size_t r = 0; r != num_complete_lines; ++r)
            {
                line_end = strchr(line_end, '\n');
                *line_end++ = '\0';
            }
        }


        char * write_pointer = fq_dat_in;
        write_pointer[0] = '\0';

        for (size_t r = 0; r != num_complete_records; ++r)
        {
            //write the fastq ID
            // there is a ++ because we want to ignore the '@' at the beginning of the fastq ID
            char * fwdslash = strchr(*sit[0], '/');
            fwdslash == NULL ?  : *fwdslash = '\0';


            write_pointer += sprintf(write_pointer, "%s", (*sit[0]) + 1); 

            for (size_t f = 0; f != nfq; ++f)
            {
                write_pointer += 
                    sprintf(write_pointer, "\t%s\t%s", *(sit[f] + 1), *(sit[f] + 3));

                sit[f] += 4;
            }
            write_pointer += sprintf(write_pointer, "\n");
        }
        // now sit[] is set to iterators pointing to first unused portion

        if (num_complete_records > 0)
        {

            char * tmp_file = new char[strlen(tmp_file_template) + 1];
            strcpy(tmp_file, tmp_file_template);
            int fdes = mkstemp(tmp_file);
            
            tmp_files.push_back(tmp_file);
            FILE * ftmp = fdopen(fdes, "w+");
            if (ftmp == NULL)
            {
                fprintf(stderr, "Error: couldn't open temporary chunk file %s for reading/writing.\n",
                        tmp_files[c]);
                exit(1);
            }
            
            tmp_fhs.push_back(ftmp);

            char * last_fragment;

            std::vector<char *> lines = 
                FileUtils::find_complete_lines_nullify(fq_dat_in, & last_fragment);

            assert(strlen(last_fragment) == 0);

            std::pair<size_t, size_t> chunk_info = 
                process_chunk(lines, fq_dat_in, fq_dat_out, sam_order, ftmp, & line_index);
            
            //offset_quantile_sizes.push_back(chunk_info.first);
            chunk_num_lines.push_back(chunk_info.second);
        }
    }
    while (num_complete_records > 0);

    delete sit;

    for (size_t f = 0; f != nfq; ++f)
    {
        if (! gzeof(fastq_fhs[f]))
        {
            fprintf(stderr, "Error: Unequal number of records in provided fastq files\n");
            exit(1);
        }
        gzclose(fastq_fhs[f]);
    }

    delete fq_dat_in;
    delete fq_dat_out;
    for (size_t f = 0; f != nfq; ++f)
    {
        delete fq_chunks[f];
    }
    

    // now we have sorted, reformatted chunk tmp files, and an (O, K)
    // sorted line_index.  offset_quantile_sizes and chunk_num_lines
    // describe the partitioning.
    size_t num_sort_chunks = chunk_num_lines.size();

    set_start_offsets(line_index.begin(), line_index.end(), 0);

    std::vector<INDEX_ITER> offset_quantiles;
    INDEX_ITER iit = line_index.begin();
    offset_quantiles.push_back(iit);
    for (size_t c = 0; c != num_sort_chunks; ++c)
    {
        std::advance(iit, chunk_num_lines[c]);
        offset_quantiles.push_back(iit);
    }
    assert(iit == line_index.end());

    for (size_t c = 0; c != num_sort_chunks; ++c)
    {
        fseek(tmp_fhs[c], 0, std::ios::beg);
    }

    // Now do the merge
    
    write_final_merge(line_index, offset_quantiles, tmp_fhs, 
                      out_data_fh, out_index_fh);
    
    fclose(out_data_fh);
    fclose(out_index_fh);

    fprintf(stderr, ".\n");
    fprintf(stderr, "Cleaning up...");
    fflush(stderr);

    for (size_t t = 0; t != tmp_files.size(); ++t)
    {
        fclose(tmp_fhs[t]);
        remove(tmp_files[t]);
        delete tmp_files[t];
    }

    fprintf(stderr, "done.\n");
    fflush(stderr);

    delete fastq_fhs;
    delete fq_chunks;

    return 0;
}
