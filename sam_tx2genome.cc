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

#include <parallel/algorithm>
#include <omp.h>

/*

projects transcriptome-aligned reads to the genome.

transcript -> genome, transcript -> read :    genome -> read

NEEDED: 

1. tx_to_genome  (projections of transcripts to the genome)
2. tx_projections.  (tx names pointing to tx_to_genome iterators)

0.  source = transcript, target = genome
1.  parse a samline
2.  retrieve the source_to_target cigar from the tx_to_genome 
3.  merge the source_to_target and source_to_read CIGARs to create a target_to_read CIGAR

If our input is samlines aligned to transcripts, and the transcripts are ordered
by start on the genome, then we just keep track of new transcripts in the SamLine input.

*/
int tx2genome_usage(size_t mdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil tx2genome [OPTIONS] transcripts.gtf genome.samheader reads.vs.tx.asort.rsam reads.vs.genome.sam\n\n"
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
            "reads.vs.tx.rsam must be sorted by [rname, read_pair_flag, pos].\n"
            "\n"
            "reads.vs.tx.asort.sam must be sorted by alignment position, with contigs (transcript)\n"
            "ordered by projection order according to genomic contig order as specified.\n"
            "expected read layout is a sequenced of 'f' and 'r' denoting whether each read is\n"
            "on the forward or reverse strand of the template molecule.\n",
            mdef
            );
    return 1;
}


PROJ_MAP 
LoadProjectionMap(std::set<SequenceProjection, less_seq_projection> const& tx_to_genome,
                  SamOrder const& genome_sam_order);


//call at the end of loading the input_buffer with all alignments to a
//given transcript.




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

    SamOrder genome_sam_order(SAM_POSITION_RID, "ALIGN");
    SamOrder tx_sam_order(SAM_POSITION_RID, "ALIGN");

    char * tx_header_buf = ReadAllocSAMHeader(input_tx_sam_fh);
    tx_sam_order.AddHeaderContigStats(tx_header_buf);
    delete [] tx_header_buf;

    char * genome_header_buf = ReadAllocSAMHeader(input_genome_sam_header_fh);
    fclose(input_genome_sam_header_fh);

    size_t genome_header_length = strlen(genome_header_buf);
    genome_sam_order.AddHeaderContigStats(genome_header_buf);
    fwrite(genome_header_buf, 1, genome_header_length, output_genome_sam_fh);
    fflush(output_genome_sam_fh);
    delete [] genome_header_buf;


    char const* species = "dummy";

    std::set<SequenceProjection> tx_to_genome_tmp = 
        gtf_to_sequence_projection(gtf_file, species);

    less_seq_projection lsp(& genome_sam_order);
    std::set<SequenceProjection, less_seq_projection> tx_to_genome(lsp);
    std::set<SequenceProjection, less_seq_projection>::iterator tx_iter;
    PROJ_MAP tx_projections;

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

        if (! SamLine::initialized && ! sam_lines.empty())
        {
            // initialize tx_sam_order, genome_sam_order, SamLine, tx_to_genome, tx_projections
            SAM_QNAME_FORMAT tx_qname_fmt = QNAMEFormat(sam_lines[0]);
            SamLine::SetGlobalFlags(tx_qname_fmt, "", "", 0);
            tx_sam_order.InitFromChoice(tx_qname_fmt);
            genome_sam_order.InitFromChoice(tx_qname_fmt);

            tx_to_genome.insert(tx_to_genome_tmp.begin(), tx_to_genome_tmp.end());
            tx_projections = LoadProjectionMap(tx_to_genome, genome_sam_order);
        }

        // 3. convert to std::vector<SamLine *> (transform)
        size_t prev_num_records = sam_records.size();
        sam_records.resize(prev_num_records + sam_lines.size());

        __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                  sam_records.begin() + prev_num_records,
                                  parse_sam_unary());

        __gnu_parallel::for_each(sam_records.begin() + prev_num_records, 
                                 sam_records.end(), 
                                 set_flattened_pos_unary(& tx_sam_order));

        //compute average printed line length of first 100 records
        size_t average_rsam_line_length = 200; // reasonable value in absence of any other information
        size_t n_sampled_records = std::min(100UL, sam_records.size());
        char write_buffer[1024];
        size_t sum_of_lengths = 0;
        for (size_t r = 0; r != n_sampled_records; ++r)
        {
            sam_records[r]->sprint(write_buffer);
            sum_of_lengths += strlen(write_buffer);
        }
        if (n_sampled_records > 0)
        {
            average_rsam_line_length = sum_of_lengths / n_sampled_records;
        }

        // 4. make chunks at non-overlapping transcripts.
        size_t min_range = 1000; // arbitrarily set the single work unit size to 1000
        size_t est_num_work_units = 1 + (sam_lines.size() / min_range);
        std::vector<std::pair<SAMIT, SAMIT> > ranges;
        ranges.reserve(est_num_work_units);

        SAMIT cur = sam_records.begin();
        SAMIT pre;

        SEQ_PROJ_ITER cur_proj_iter;
        size_t n_records_left = sam_records.size();

        CONTIG_OFFSETS::const_iterator dummy = 
            tx_sam_order.contig_offsets.begin();

        while (cur != sam_records.end())
        {
            pre = cur;
            size_t jump = std::min(min_range, n_records_left);
            std::advance(cur, jump);

            // find the projection map iterator corresponding to this
            // record, if it exists
            PROJ_MAP::const_iterator p = 
                (cur == sam_records.end())
                ? tx_projections.end()
                : tx_projections.find((*cur)->rname);

            // find the projection iterator, if it exists
            cur_proj_iter = 
                (p == tx_projections.end()) 
                ? tx_to_genome.end()
                : (*p).second;
            
            // find the next valid bound projection iterator if it exists
            cur_proj_iter = non_overlapping_range(cur_proj_iter, tx_to_genome.end());
            
            // compute the lowest flattened position for this valid bound

            if (cur_proj_iter == tx_to_genome.end())
            {
                // do nothing
            }
            else
            {
                size_t flat_lowbound = 
                    flattened_position_aux((*cur_proj_iter).source_dna.c_str(),
                                           0,
                                           tx_sam_order.contig_offsets,
                                           & dummy);
                                           
                for ( ; cur != sam_records.end()
                          && (*cur)->flattened_pos < flat_lowbound; ++cur) { }

                assert(cur == sam_records.end() || (*cur)->flattened_pos >= flat_lowbound);
            }
            // postcondition: cur represents the valid bound
            // corresponding to the projection bound.


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

        // Where did we account for unused SAM records??? This is n_records_left

        if(ranges.empty())
        {
            if (cumul_nbytes_read + chunk_size > max_mem)
            {
                fprintf(stderr, "Error: Run cannot be completed with allotted memory: %zu\n"
                        "due to too many entries on overlapping genes. "
                        "Please run with more memory.\n",
                        max_mem);
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

        std::vector<std::vector<char> *> projected_sam_records(ranges.size());

        project_dedup_print pdp_aux(&tx_projections,
                                    &tx_sam_order,
                                    &genome_sam_order,
                                    inserts_are_introns,
                                    retain_unsequenced_projection,
                                    average_rsam_line_length);

        // 5. project, load to buffer, string_alloc, string_print, delete (transform)
        __gnu_parallel::transform(ranges.begin(), ranges.end(),
                                  projected_sam_records.begin(), pdp_aux);
                                  
        // 6. delete allocated vectors
        for (size_t r = 0; r != ranges.size(); ++r)
        {
            fwrite(&(*projected_sam_records[r])[0], 1, 
                   projected_sam_records[r]->size(), output_genome_sam_fh);

            delete projected_sam_records[r];
        }
        // now we want to copy the range [pre, end) into sam_records
        sam_records.swap(std::vector<SamLine *>(pre, sam_records.end()));

        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

        cumul_nbytes_read = 0;
    }

    delete [] chunk_buffer_in;

    fclose(input_tx_sam_fh);
    fclose(output_genome_sam_fh);

    for (NP_ITER tn = tx_projections.begin(); tn != tx_projections.end(); ++tn)
    {
        delete (*tn).first;
    }

    
    return 0;
}


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


//check that the genome contig order implied by the combination of
//transcript projections and transcript contig order is the same as
//that given by genome_contig_order
bool CheckProjectionOrder(CONTIG_OFFSETS const& genome_contig_order,
                          CONTIG_OFFSETS const& tx_contig_order,
                          std::set<SequenceProjection> const& tx_to_genome)
{
    //1. Compute a flattened genome coordinate for each transcript
    std::multimap<size_t, char const*> tx_implied_genomic_offsets;
    std::multimap<size_t, char const*>::const_iterator tit;

    std::set<SequenceProjection>::const_iterator pit;
    OFFSETS_ITER oit;
    for (pit = tx_to_genome.begin(); pit != tx_to_genome.end(); ++pit)
    {
        char const* tx_contig = (*pit).source_dna.c_str();
        char const* tx_proj_contig = (*pit).target_dna.c_str();
        oit = genome_contig_order.find(tx_proj_contig);
        if (oit == genome_contig_order.end())
        {
            fprintf(stderr, "Error: CheckProjectionOrder: didn't find genome contig %s "
                    "that was found in transcript-to-genome projection from gtf file\n",
                    tx_proj_contig);
            exit(1);
        }
        size_t genome_contig_offset = (*oit).second;
        size_t tx_local_offset = (*pit).transformation[0].jump_length;

        tx_implied_genomic_offsets.insert(std::make_pair(genome_contig_offset + tx_local_offset, tx_contig));
    }

    size_t last_implied_offset = SIZE_MAX;
    size_t last_given_offset = 0;
    char const* last_tx_contig = "";

    //tx implied offsets are partially ordered.  given offsets are
    //fully ordered.  here, we merely check that if a consecutive pair
    //of contigs are implied strictly increasing, they are increasing
    //by given offset as well. 
    for (tit = tx_implied_genomic_offsets.begin(); tit != tx_implied_genomic_offsets.end(); ++tit)
    {
        size_t implied_offset = (*tit).first;
        char const* tx_contig = (*tit).second;
        oit = tx_contig_order.find(tx_contig);
        if (oit == tx_contig_order.end())
        {
            fprintf(stderr, "Error: CheckProjectionOrder: didn't find transcript contig %s "
                    "in given transcript ordering that was found in transcript-to-genome "
                    "projection from gtf file\n",
                    tx_contig);
            exit(1);
        }
        size_t given_offset = (*oit).second;
        if (last_implied_offset < implied_offset
            && ! (last_given_offset < given_offset))
        {
            fprintf(stderr, "Error: CheckProjectionOrder: transcript contig ordering implied from "
                    "gtf file and given genome contig ordering differs from "
                    "given transcript ordering. Please re-order transcript contigs in SAM header\n"
                    "Last contig: %s (%Zu implied offset) (%Zu given offset)\n"
                    "Current ctg: %s (%Zu implied offset) (%Zu given offset)\n",
                    last_tx_contig, last_implied_offset, last_given_offset, 
                    tx_contig, implied_offset, given_offset);

            return false;
        }
        else
        {
            last_implied_offset = implied_offset;
            last_given_offset = given_offset;
            last_tx_contig = tx_contig;
        }
    }
    return true;
}



