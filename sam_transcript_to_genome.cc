#define __STDC_LIMIT_MACROS
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

/*

projects transcriptome-aligned reads to the genome.

transcript -> genome, transcript -> read :    genome -> read

NEEDED: 

1. tx_to_genome  (projections of transcripts to the genome)
2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)

0.  source = transcript, target = genome
1.  parse a samline
2.  retrieve the source_to_target cigar from the tx_to_genome 
3.  merge the source_to_target and source_to_read CIGARs to create a target_to_read CIGAR

If our input is samlines aligned to transcripts, and the transcripts are ordered
by start on the genome, then we just keep track of new transcripts in the SamLine input.

*/
int transcript_to_genome_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil tx2genome [OPTIONS] species annotations.gtf transcript_alignment.sam genome_header.sam projected_genome_alignment.sam\n\n"
            "Options:\n\n"
            "-h      STRING    genome-space SAM header file to be included in output\n"
            "\n"
            );
    return 1;
}

int main_transcript_to_genome(int argc, char ** argv)
{
    char c;

    char const* input_genome_sam_header_file = "";

    while ((c = getopt(argc, argv, "h:")) >= 0)
    {
        switch(c)
        {
        case 'h': input_genome_sam_header_file = optarg; break;
        default: return transcript_to_genome_usage(); break;
        }
    }

    //printf("optind: %i\n", optind);
    if (argc != 4 + optind)
    {
        return transcript_to_genome_usage();
    }

    char * species = argv[optind];
    char * gtf_file = argv[optind + 1];
    char * input_tx_sam_file = argv[optind + 2];
    char * output_genome_sam_file = argv[optind + 3];

    FILE * input_tx_sam_fh = open_or_die(input_tx_sam_file, "r", "Input transcript-based SAM file");
    FILE * input_genome_sam_header_fh = open_if_present(input_genome_sam_header_file, "r");
    FILE * output_genome_sam_fh = open_or_die(output_genome_sam_file, "w", "Output genome-projected SAM file");


    typedef std::set<SequenceProjection>::const_iterator SP_ITER;

    bool ones_based_pos = true;

    // input SAM file is allowed to have or not have seq / qual.
    // basically, they are just payload, and unaffected by the transformation
    bool allow_absent_seq_qual = true;

    //bool flip_query_strand_flag = false;
    SAM_ORDER sam_order = SAM_POSITION_RID;
    bool paired_reads_are_same_stranded = false;

    std::set<SequenceProjection> genome_to_tx = 
        gtf_to_sequence_projection(gtf_file, species);

    //needed if we are expanding reads that are originally mapped to a
    //transcript to genomic coordinates
    std::set<SequenceProjection> tx_to_genome;

    std::set<cis::dna_t> contigs;

    typedef std::set<cis::dna_t>::const_iterator DNAS_ITER;
    DNAS_ITER contig_iter, prev_contig_iter;

    for (SP_ITER gt_iter = genome_to_tx.begin();
         gt_iter != genome_to_tx.end(); ++gt_iter)
    {
        contigs.insert(cis::dna_t((*gt_iter).species, (*gt_iter).source_dna,
                                  0, INT64_MAX));
    }


    cis::TREE_MAP transcript_extent_trees;
    cis::REGIONS_MULTI transcript_extents;

    //this maps tx_extent regions to sequence projections.  when in
    //expansion mode, the sequence projection will be
    //transcript->genome 't2g', and be iterators into tx_to_genome.
    //otherwise, these will be iterators into genome_to_tx
    typedef std::map<char const*, SP_ITER, less_char_ptr> PROJ_MAP;
    PROJ_MAP tx_name_to_projection;
    typedef PROJ_MAP::const_iterator NP_ITER;

    std::map<cis::region const*, SP_ITER> tx_extent_to_projection;
    typedef std::map<cis::region const*, SP_ITER>::const_iterator MP_ITER;

    typedef std::pair<cis::REGIONS_MULTI::iterator, bool> REG_INS;

    // 1. tx_to_genome  (projections of transcripts to the genome)
    // 2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)
    for (SP_ITER gt_iter = genome_to_tx.begin();
         gt_iter != genome_to_tx.end(); ++gt_iter)
    {
        SP_ITER tg_iter = tx_to_genome.insert(tx_to_genome.begin(), 
                                              InvertProjection((*gt_iter)));
        
        //transcript-to-genome projection for this transcript
        SequenceProjection const& tx_t2g_projection = *tg_iter;
        
        char * tx_name = new char[tx_t2g_projection.source_dna.size() + 1];
        strcpy(tx_name, tx_t2g_projection.source_dna.c_str());
        tx_name_to_projection[tx_name] = tg_iter;
    }

    SP_ITER first_proj = tx_to_genome.end();
    SP_ITER second_proj = tx_to_genome.end();

    bool ignore_duplicate_mapped_pairs = false;

    char prev_rname[1024] = "";

    SamLine * samline;
    SamBuffer input_buffer(sam_order,
                           paired_reads_are_same_stranded,
                           ones_based_pos,
                           ignore_duplicate_mapped_pairs);
    

    SamBuffer output_buffer(sam_order,
                            paired_reads_are_same_stranded,
                            ones_based_pos,
                            ignore_duplicate_mapped_pairs);
    
    
    contig_iter = contigs.end();

    std::map<std::string, size_t> tx_lengths = ContigLengths(input_tx_sam_fh);
    CONTIG_OFFSETS tx_offsets = ContigOffsets(tx_lengths);
    
    if (input_genome_sam_header_fh != NULL)
    {
        PrintSAMHeader(&input_genome_sam_header_fh, output_genome_sam_fh);
    }
    close_if_present(input_genome_sam_header_fh);

    char fake_samline[1024];

    size_t prev_pos_index = 0;
    size_t cur_pos_index = 0;

    bool inserted_to_input;

    while (! feof(input_tx_sam_fh))
    {
        // parse samline
        samline = new SamLine(input_tx_sam_fh, ones_based_pos, allow_absent_seq_qual);

        switch (samline->parse_flag)
        {
        case END_OF_FILE: 
        case HEADER:
            delete samline;
            break;

        case PARSE_ERROR:
            fprintf(stderr, "Parse error in input sam file %s", input_tx_sam_file);
            exit(1);
            break;

        case DATA_LINE:

            sprintf(fake_samline, "fake\t%i\t%s\t%Zu", samline->flag, samline->rname, samline->pos);

            cur_pos_index = samline_position_align(fake_samline, tx_offsets);
            if (cur_pos_index < prev_pos_index)
            {
                fprintf(stderr, "Error: input is not sorted by [rname, flag, pos] fields\n");
                exit(1);
            }
            prev_pos_index = cur_pos_index;

            if (samline->query_unmapped())
            {
                samline->print(output_genome_sam_fh, true);
                delete samline;
                inserted_to_input = false;
            }
            else
            {
                inserted_to_input = input_buffer.insert(samline);
            }

            if (inserted_to_input && strcmp(prev_rname, samline->rname) != 0)
            {
                //on a new contig.  
                // if (input_buffer.low_bound != NULL)
                // {
                //     delete input_buffer.low_bound;
                // }
                // input_buffer.low_bound = 
                //     new SamLine(DATA_LINE, "0:", 0, samline->rname, 0, 0, "", "*", 0, 0, "", "", "");

                if (! input_buffer.unique_entry_pairs.empty())
                {
                    //find projection by name
                    //assume all unique entry pairs in input buffer are on the same
                    //transcript, namely 'prev_rname'.
                    NP_ITER nit = tx_name_to_projection.find(prev_rname);
                    if (nit == tx_name_to_projection.end())
                    {
                        fprintf(stderr, "Error: couldn't find transcript %s in projections "
                                "derived from GTF file\n", prev_rname);
                        exit(1);
                    }
                    else
                    {
                        first_proj = (*nit).second;
                    }
                    
                    // PAIRED_READ_SET::const_iterator paired_bound =
                    //     input_buffer.unique_entry_pairs.lower_bound
                    //     (std::make_pair(input_buffer.low_bound, input_buffer.low_bound));
                    
                    for (PAIRED_READ_SET::iterator pit = input_buffer.unique_entry_pairs.begin();
                         pit != input_buffer.unique_entry_pairs.end(); ++pit)
                    {

                        SamLine * first = const_cast<SamLine *>((*pit).first);
                        SamLine * second = const_cast<SamLine *>((*pit).second);

                        assert(strcmp(prev_rname, first->rname) == 0);

                        //but, their mates may not be
                        NP_ITER nit_mate = strcmp(first->rname, first->mate_ref_name()) == 0 ? nit 
                            : tx_name_to_projection.find(first->mate_ref_name());
                    
                        if (nit_mate == tx_name_to_projection.end())
                        {
                            fprintf(stderr, "Error: couldn't find transcript %s in projections "
                                    "derived from GTF file\n", samline->mate_ref_name());
                            exit(1);
                        }
                        else
                        {
                            second_proj = (*nit_mate).second;
                        }
                    
                        //extract all completed input pairs
                        bool projection_applied = 
                            ApplyProjectionToSAM(*first_proj, *second_proj, first, second);
                        
                        assert(projection_applied);

                        assert(AreMappedMatePairs(*first, *second));

                        output_buffer.insert(first);
                        output_buffer.insert(second);
                        char const* target_contig = (*first_proj).target_dna.c_str();
                        size_t target_min_pos = (*first_proj).cigar[0].length;
                        // if (! first->query_unmapped())
                        // {
                        //     assert(target_min_pos <= first->pos);
                        // }
                    
                    }
                    input_buffer.unique_entry_pairs.clear();
                
                    //at this point, the low bound
                    char const* target_contig = (*first_proj).target_dna.c_str();
                    size_t target_min_pos = (*first_proj).cigar[0].length;

                    if (output_buffer.low_bound != NULL)
                    {
                        delete output_buffer.low_bound;
                    }
                    output_buffer.low_bound =
                        new SamLine(DATA_LINE, "0:", 0, target_contig, target_min_pos, 0, "", "*", 0, 0, "", "", "");

                } // finished processing unique entry pairs up to the low bound
                strcpy(prev_rname, samline->rname);
            } // if on new contig
            
            output_buffer.purge(output_genome_sam_fh, NULL, NULL, false);
        }
    }
    output_buffer.purge(output_genome_sam_fh, NULL, NULL, true);

    delete output_buffer.low_bound;

    fclose(input_tx_sam_fh);
    fclose(output_genome_sam_fh);

    for (NP_ITER tn = tx_name_to_projection.begin();
         tn != tx_name_to_projection.end(); ++tn)
    {
        delete (*tn).first;
    }
    
    for (MP_ITER tex = tx_extent_to_projection.begin();
         tex != tx_extent_to_projection.end(); ++tex)
    {
        delete (*tex).first;
    }

    return 0;
}
