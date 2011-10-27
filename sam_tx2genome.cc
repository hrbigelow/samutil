#define __STDC_LIMIT_MACROS
#include <pthread.h>
#include <stdint.h>
#include <algorithm>
#include <map>
#include <set>

#include "readsim_aux.h"
#include "cisortho/nested.h"
#include "dep/tools.h"
#include "cisortho/region.h"
#include "sam_buffer.h"
#include "sam_helper.h"
#include "sam_score_aux.h"
#include "sam_aux.h"
#include "seq_projection.h"

#include <parallel/algorithm>
#include <omp.h>

int tx2genome_usage(size_t mdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil tx2genome [OPTIONS] transcripts.gtf genome.samheader reads.vs.tx.psort.rsam reads.vs.genome.sam\n\n"
            "Options:\n\n"
            "-n     FLAG     If present, use CIGAR 'N' to represent introns.  Otherwise use 'D'\n"
            "-x     FLAG     If present, add Cufflinks XS:A: tag (+/-) denotes source RNA strand\n"
            "-u     FLAG     If present, retain inferred alignment of unsequenced portion [false]\n"
            "-t     INT      Number of threads to use\n"
            "-m     INT      number bytes of memory to use [%Zu]\n"
            "-C     STRING   work in the directory named here [.]\n"
            "\n"
            "Only one of each group of identical fragment alignments is output.\n"
            "These arise from congruent subsets of isoforms in transcripts.gtf.\n"
            "SAM Records successfully projected will have the 'XP:A:T' tag added.\n"
            "reads.vs.tx.psort.rsam must be sorted by PROJALIGN (see samutil sort)\n"
            "\n"
            "expected read layout is a sequenced of 'f' and 'r' denoting\n"
            "whether each read is on the forward or reverse strand of the\n"
            "template molecule.\n"
            "\n",
            mdef
            );
    return 1;
}


// PROJ_MAP 
// LoadProjectionMap(std::set<SequenceProjection, less_seq_projection> const& tx_to_genome,
//                   SamOrder const& genome_sam_order);



int main_tx2genome(int argc, char ** argv)
{
    char c;

    bool inserts_are_introns = false;
    bool add_cufflinks_xs_tag = false;
    bool retain_unsequenced_projection = false;

    char const* working_dir = ".";

    size_t num_threads = 1;
    size_t mdef = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = mdef;


    while ((c = getopt(argc, argv, "nxut:m:C:")) >= 0)
    {
        switch(c)
        {
        case 'n': inserts_are_introns = true; break;
        case 'x': add_cufflinks_xs_tag = true; break;
        case 'u': retain_unsequenced_projection = true; break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        default: return tx2genome_usage(mdef); break;
        }
    }

    if (argc != 4 + optind)
    {
        return tx2genome_usage(mdef);
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    char * gtf_file = argv[optind];
    char * input_genome_sam_header_file = argv[optind + 1];
    char * input_tx_sam_file = argv[optind + 2];
    char * output_genome_sam_file = argv[optind + 3];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * input_tx_sam_fh = open_or_die(input_tx_sam_file, "r", "Input transcript-based SAM file");
    FILE * input_genome_sam_header_fh = open_or_die(input_genome_sam_header_file, "r", "Input genome SAM header file");
    FILE * output_genome_sam_fh = open_or_die(output_genome_sam_file, "w", "Output genome-projected SAM file");

    SamOrder genome_sam_order(SAM_POSITION_RID, "PROJALIGN");

    // just to slurp up the tx header.
    char * tx_header_buf = ReadAllocSAMHeader(input_tx_sam_fh);
    delete [] tx_header_buf;

    char * genome_header_buf = ReadAllocSAMHeader(input_genome_sam_header_fh);
    fclose(input_genome_sam_header_fh);

    size_t genome_header_length = strlen(genome_header_buf);
    genome_sam_order.AddHeaderContigStats(genome_header_buf);
    fwrite(genome_header_buf, 1, genome_header_length, output_genome_sam_fh);
    fflush(output_genome_sam_fh);
    delete [] genome_header_buf;
    genome_sam_order.InitProjection(gtf_file);

    fprintf(stderr, "Starting Projection...\n");
    fflush(stderr);

    std::vector<char *> sam_lines;
    std::vector<char *>::iterator sit;

    size_t nbytes_unused = 0;
    size_t nbytes_read;
    char * last_fragment;

    // at any one point, we will have:  a chunk_buffer, a buffer of converted chunks,
    size_t chunk_size = std::min(1024UL * 1024UL * 1024UL, max_mem / 2);
    char * chunk_buffer_in = new char[chunk_size + 1];
    char * read_pointer = chunk_buffer_in;

    std::vector<SamLine *> sam_records;

    // number of bytes loaded that haven't been purged
    size_t cumul_nbytes_read = 0;
    size_t num_records_processed = 0;
    size_t num_records_discarded = 0;
    size_t num_records_retained = 0;

    while (! feof(input_tx_sam_fh))
    {
        /*
          itinerary:
          1. read a raw chunk
          2. find complete newlines, recycle unused.
          3. convert to std::vector<SamLine *> (transform)
          4. make ~1000 chunks at non-overlapping transcripts
          5. project, load to buffer, string_alloc, print, delete (transform)
          6. delete sam_buffers and alloc strings
         */

        // 1. read a raw chunk
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, input_tx_sam_fh);
        read_pointer[nbytes_read] = '\0';
        cumul_nbytes_read += nbytes_read;

        // 2. find complete newlines, recycle unused.
        sam_lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        num_records_processed += sam_lines.size();

        if (! SamLine::initialized && ! sam_lines.empty())
        {
            // initialize tx_sam_order, genome_sam_order, SamLine, tx_to_genome, tx_projections
            SAM_QNAME_FORMAT tx_qname_fmt = QNAMEFormat(sam_lines[0]);
            SamLine::SetGlobalFlags(tx_qname_fmt, "", "", 0, false);
            genome_sam_order.InitFromChoice(tx_qname_fmt);
        }

        // 3. convert to std::vector<SamLine *> (transform)
        size_t prev_num_records = sam_records.size();
        sam_records.resize(prev_num_records + sam_lines.size());

        __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                  sam_records.begin() + prev_num_records,
                                  parse_sam_unary());

        //compute average printed line length of first 100 records
        size_t n_sampled_records = 1000;
        size_t average_rsam_line_length = 
            average_line_length(sam_records.begin(), sam_records.end(), n_sampled_records);


        // 4. make chunks at non-overlapping transcripts.
        size_t min_range = 10000; // arbitrarily set the single work unit size to 1000
        size_t est_num_work_units = 1 + (sam_lines.size() / min_range);
        std::vector<std::pair<SAMIT, SAMIT> > ranges;
        ranges.reserve(est_num_work_units);

        SAMIT cur = sam_records.begin();
        SAMIT pre;

        size_t n_records_left = sam_records.size();

        // traverse 'sam_records'. find all ranges [start, end) that
        // are non-overlapping and thus safe to process individually
        // in SAM buffers. A range bound can either be a record having
        // a greater flat_pos than the previous one, or it can have
        // the same flat_pos but be unmapped. must tweak PROJALIGN
        // indexing scheme so it can partially order unmapped records
        // by fragment_id.
        
        while (cur != sam_records.end())
        {
            // jump minimum distance
            pre = cur;
            size_t jump = std::min(min_range, n_records_left - 1);
            std::advance(cur, jump);

            assert(n_records_left == static_cast<size_t>(std::distance(pre, sam_records.end())));

            SamLine * cur_rec = new SamLine(**cur);

            char cur_rec_string[1024];

            size_t last_rec_index = SIZE_MAX;
            size_t cur_rec_index = SIZE_MAX;
            
            while (last_rec_index >= cur_rec_index)
            {
                last_rec_index = cur_rec_index;

                if (++cur == sam_records.end())
                {
                    // by definition, the 'end' is a valid bound
                    break;
                }

                delete cur_rec;
                cur_rec = new SamLine(**cur);
                cur_rec->sprint(cur_rec_string);
                cur_rec_index = (genome_sam_order.*(genome_sam_order.sam_index))(cur_rec_string);

            }
            fprintf(stderr, "last: %zu, cur: %zu, jump: %zu\n",
                    last_rec_index, cur_rec_index, cur_rec_index - last_rec_index);
            
            delete cur_rec;

            // update with the valid range
            if (cur != sam_records.end() || feof(input_tx_sam_fh))
            {
                // either cur is a valid boundary or it's the end
                ranges.push_back(std::make_pair(pre, cur));
                n_records_left -= std::distance(pre, cur);
            }

        }        

        // postcondition: pre is the last valid bound found.

        // unless we're at the end of the file, we don't want to use
        // up all remaining records. (because we can't find a valid
        // bound then)
        assert(feof(input_tx_sam_fh) || n_records_left > 0);

        // at this point, recover the partially read line back to the
        // start of the buffer.  sam_records is not affected by this.
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

        // Where did we account for unused SAM records??? This is n_records_left
        if(ranges.empty())
        {
            assert(pre == sam_records.begin());
            // we need to load more records.
            if (cumul_nbytes_read + chunk_size > max_mem)
            {
                char first[1024];
                char last[1024];
                (*sam_records.begin())->sprint(first);
                (*sam_records.rbegin())->sprint(last);
                fprintf(stderr, 
                        "Error: Run cannot be completed with allotted memory %zu bytes\n"
                        "due to too a clump of too many entries on overlapping genes.\n"
                        "First entry in clump: %s"
                        "Last entry in clump: %s"
                        "Number entries in clump: %zu\n"
                        "Please run with more memory.\n",
                        max_mem,
                        first, last, sam_records.size());
                exit(1);
            }
            else
            {
                continue;
            }
        }

        __gnu_parallel::_Settings psettings;
        psettings.transform_minimal_n = 2;
        __gnu_parallel::_Settings::set(psettings);

        std::vector<pdp_result> projected_sam_records(ranges.size());

        project_dedup_print pdp_aux(&genome_sam_order,
                                    inserts_are_introns,
                                    retain_unsequenced_projection,
                                    average_rsam_line_length);

        // 5. project, load to buffer, string_alloc, string_print, delete (transform)
        __gnu_parallel::transform(ranges.begin(), ranges.end(),
                                  projected_sam_records.begin(), pdp_aux);
                                  
        // 6. delete allocated vectors
        for (size_t r = 0; r != ranges.size(); ++r)
        {
            fwrite(&(*projected_sam_records[r].lines)[0], 1, 
                   projected_sam_records[r].lines->size(), 
                   output_genome_sam_fh);

            delete projected_sam_records[r].lines;
            num_records_retained += projected_sam_records[r].num_records_retained;
            num_records_discarded += projected_sam_records[r].num_records_discarded;
        }
        // now we want to copy the range [pre, end) into sam_records
        sam_records.swap(std::vector<SamLine *>(pre, sam_records.end()));

        cumul_nbytes_read = 0;
    }

    delete [] chunk_buffer_in;

    fclose(input_tx_sam_fh);
    fclose(output_genome_sam_fh);

    fprintf(stderr, 
            "Records processed: %zu\n"
            "Records discarded: %zu\n"
            "Records retained : %zu\n",
            num_records_processed,
            num_records_discarded,
            num_records_retained);
    
    return 0;
}

/*
PROJ_MAP LoadProjectionMap(std::set<SequenceProjection, less_seq_projection> const& tx_to_genome,
                           SamOrder const& genome_sam_order)
{
    // 1. tx_to_genome  (projections of transcripts to the genome)
    // 2. tx_projections.  (tx names pointing to tx_to_genome iterators)
    PROJ_MAP tx_projections;

    CONTIG_OFFSETS::const_iterator genome_offset_iter = 
        genome_sam_order.contig_offsets.begin();

    size_t prev_flat_start_pos = 0;
    for (SEQ_PROJ_ITER t = tx_to_genome.begin(); t != tx_to_genome.end(); ++t)
    {
        char * tx_name = new char[(*t).source_dna.size() + 1];
        strcpy(tx_name, (*t).source_dna.c_str());
        std::pair<NP_ITER, bool> ins = tx_projections.insert(std::make_pair(tx_name, t));

        SequenceProjection const& p2 = (*t);
                
        if (! ins.second)
        {
            SequenceProjection const& p1 = *(*ins.first).second;
            SequenceProjection const& p2 = (*t);
                    
            fprintf(stderr, "Error: duplicate transcript named %s encountered.\n"
                    "First instance: %s\t%s\n"
                    "Second instance: %s\t%ss\n", 
                    tx_name,
                    p1.source_dna.c_str(), (p1.same_strand ? "+" : "-"),
                    p2.source_dna.c_str(), (p2.same_strand ? "+" : "-"));
            exit(1);
        }


        size_t flat_start_pos = 
            flattened_position_aux(p2.target_dna.c_str(), p2.target_start_pos(),
                                   genome_sam_order.contig_offsets,
                                   & genome_offset_iter);
                
        if (flat_start_pos < prev_flat_start_pos)
        {
            fprintf(stderr, "Error: previous projection %s has target start position\n"
                    "before current one (%s).\n", 
                    (*t).source_dna.c_str(),
                    (*(--SEQ_PROJ_ITER(t))).source_dna.c_str());
            exit(1);
        }
        prev_flat_start_pos = flat_start_pos;
    }

    return tx_projections;
}
*/
