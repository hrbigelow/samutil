//Project a set of SAM-formatted alignments from one space to another,
//defined by a coord transform file

//in gcc 4.4.4 this doesn't work to define INT64_MAX.  Using INT_MAX instead
//#define __STDC_LIMIT_MACROS
//#include <stdint.h>

#include <limits.h>
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

Takes a set of reads aligned to one coordinate system and
projects them to the other.

By default, projects reads aligned to a transcriptome to the genome.

If the -c flag is set, 
MERGE:


transcript -> genome, transcript -> read :    genome -> read
genome -> transcript, genome -> read     :    transcript -> read


In expansion mode:
input samlines are alignments to individual transcripts
output samlines are alignments to the genome
NEEDED: 

1. tx_to_genome  (projections of transcripts to the genome)
2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)

0.  source = transcript, target = genome
1.  parse a samline
2.  retrieve the source_to_target cigar from the tx_to_genome 
3.  merge the source_to_target and source_to_read CIGARs to create a target_to_read CIGAR

For output buffering, we want to keep track of a low_bound, below which
we will not see any new reads from a sorted input.
If our input is samlines aligned to transcripts, and the transcripts are ordered
by start on the genome, then we just keep track of new transcripts in the SamLine input.

In contraction mode:
input samlines are alignments to the genome
output samlines are alignments to the transcriptome
NEEDED:
1. genome_to_tx  (projections of the genome to transcripts)
2. transcript_extent_tree  (nclist of regions which are trascript extents on the genome)
3. tx_extent_to_projection   (regions pointing to genome_to_tx iterators)


0.  source = genome, target = transcript
0a. build the transcript_extent_tree on the genome, with dnas being the genome.
1.  parse a samline
2.  look up the set of containing transcripts in the tree
3.  for each containing transcript, retrieve the source_to_target CIGAR from tx_extent_to_projection
4.  merge the source_to_target and source_to_read CIGARs to create a target_to_read CIGAR
    (source here is 'genome', target is 'transcript')

Here, from the set of containing transcripts, find the lowest, and set the current_transcript
to this.  (When the current_transcript changes, the low_bound can be updated.



Separately, need some way to detect advancement to a new region.

  
 */

int usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "project_sam [OPTIONS] genome_to_transcript.txt alignment.sam projected_alignment.sam\n\n"
            "Options:\n\n"
            "-c  FLAG    if set, perform an expansion from transcript to genome.\n"
            "            if not set, perform a contraction from genome to transcript.\n"
            "\n\n"
            );
    return 1;
}

int main(int argc, char ** argv)
{
    bool expansion_mode = false;
    char c;
    while ((c = getopt(argc, argv, "c")) >= 0)
    {
        switch(c)
        {
        case 'c': expansion_mode = true; break;
        }
    }

    //printf("optind: %i\n", optind);
    if (argc != 3 + optind)
    {
        return usage();
    }

    char * coord_transform_file = argv[optind];
    char * input_sam_file = argv[optind + 1];
    char * output_sam_file = argv[optind + 2];

    FILE * input_sam_fh = open_if_present(input_sam_file, "r");
    FILE * output_sam_fh = open_if_present(output_sam_file, "w");

    typedef std::set<SequenceProjection>::const_iterator SP_ITER;

    bool ones_based_pos = true;
    //bool flip_query_strand_flag = false;
    bool use_read_ids_for_uniqueness = true;
    bool paired_reads_are_same_stranded = false;


    std::set<SequenceProjection> genome_to_tx =
        load_sequence_projection(coord_transform_file);

    //needed if we are expanding reads that are originally mapped to a
    //transcript to genomic coordinates
    std::set<SequenceProjection> tx_to_genome;

    std::set<cis::dna_t> contigs;

    typedef std::set<cis::dna_t>::const_iterator DNAS_ITER;
    DNAS_ITER contig_iter, prev_contig_iter;

    for (SP_ITER gt_iter = genome_to_tx.begin();
         gt_iter != genome_to_tx.end(); ++gt_iter)
    {
        contigs.insert(cis::dna_t((*gt_iter).species(), (*gt_iter).source_dna(),
                                  0, INT_MAX));
    }


    cis::TREE_MAP transcript_extent_trees;
    cis::REGIONS_MULTI transcript_extents;

    SP_ITER min_transcript_proj, prev_min_transcript_proj;

    //this maps tx_extent regions to sequence projections.  when in
    //expansion mode, the sequence projection will be
    //transcript->genome 't2g', and be iterators into tx_to_genome.
    //otherwise, these will be iterators into genome_to_tx
    std::map<char const*, SP_ITER> tx_name_to_projection;
    typedef std::map<char const*, SP_ITER>::const_iterator NP_ITER;

    std::map<cis::region const*, SP_ITER> tx_extent_to_projection;
    typedef std::map<cis::region const*, SP_ITER>::const_iterator MP_ITER;

    //assume we are using only one species, since a SAMfile doesn't support
    //multiple species. (unless you encode them in the 'SN' field.
    std::string species = (*genome_to_tx.begin()).species();

    typedef std::pair<cis::REGIONS_MULTI::iterator, bool> REG_INS;

    if (expansion_mode)
    {
        // 1. tx_to_genome  (projections of transcripts to the genome)
        // 2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)
        for (SP_ITER gt_iter = genome_to_tx.begin();
             gt_iter != genome_to_tx.end(); ++gt_iter)
        {
            SP_ITER tg_iter = tx_to_genome.insert(tx_to_genome.begin(), 
                                                  InvertProjection((*gt_iter)));

            SequenceProjection const& tx_g2t_projection = *tg_iter;

            char * tx_name = new char[tx_g2t_projection.source_dna().size() + 1];
            strcpy(tx_name, tx_g2t_projection.source_dna().c_str());
            tx_name_to_projection[tx_name] = gt_iter;
        }
    }
    else
    {
        /* in contraction mode. initialize the following:
           1. genome_to_tx  (projections of the genome to transcripts)
           2. transcript_extent_tree  (nclist of regions which are trascript extents on the genome)
           3. tx_extent_to_projection   (regions pointing to genome_to_tx iterators)
        */
        Cigar::CIGAR_ITER trimmed_left, trimmed_right;
        for (SP_ITER gt_iter = genome_to_tx.begin();
             gt_iter != genome_to_tx.end(); ++gt_iter)
        {
            
            SequenceProjection const& tx_g2t_projection = *gt_iter;
            contig_iter = contigs.find(cis::dna_t(species, tx_g2t_projection.source_dna()));
            assert(contig_iter != contigs.end());
            
            int64_t tx_start_on_genome = Cigar::LeftOffset(tx_g2t_projection.cigar, true);
            Cigar::Trim(tx_g2t_projection.cigar, false, &trimmed_left, &trimmed_right);
            int64_t tx_span_on_genome = Cigar::Length(trimmed_left, trimmed_right, true);
            
            cis::dna_strand tx_strand = 
                tx_g2t_projection.same_strand ? cis::POS : cis::NEG;
            
            cis::region const* transcript_extent = 
                new cis::region(*contig_iter,
                                std::string(tx_g2t_projection.target_dna()), 
                                tx_start_on_genome, 
                                tx_start_on_genome + tx_span_on_genome, 
                                tx_strand, 0, 0, 0, NULL);
            
            transcript_extents.insert(transcript_extent);
            tx_extent_to_projection[transcript_extent] = gt_iter;
        }
        transcript_extent_trees = cis::BuildTrees(transcript_extents);
    }


    SamLine * samline;
    SamBuffer sam_buffer(use_read_ids_for_uniqueness,
                         paired_reads_are_same_stranded,
                         ones_based_pos);

    
    contig_iter = contigs.end();

    while (! feof(input_sam_fh))
    {
        // parse samline
        samline = new SamLine(input_sam_fh, ones_based_pos);
        
        switch (samline->parse_flag)
        {
        case END_OF_FILE: 
            delete samline;
            break;

        case PARSE_ERROR:
            fprintf(stderr, "Parse error in input sam file %s", input_sam_file);
            exit(1);
            break;
        case HEADER:
            samline->print(output_sam_fh, ones_based_pos, true);
            delete samline;
            break;
        case DATA_LINE:
            prev_min_transcript_proj = min_transcript_proj;
            Cigar::CIGAR_VEC source_to_read = 
                Cigar::FromString(samline->cigar, samline->pos);

            if (expansion_mode)
            {
                //find projection by name
                min_transcript_proj = tx_name_to_projection[samline->rname];
                bool projection_applied = 
                    ApplyProjectionToSAM(*min_transcript_proj, samline);

                if (projection_applied)
                {
                    sam_buffer.insert(samline);
                }
            }
            else
            {
                //find all projections by alignment bounds
                Cigar::CIGAR_ITER trimmed_left, trimmed_right;
                Cigar::Trim(source_to_read, false, &trimmed_left, &trimmed_right);
                size_t read_span_on_source = Cigar::Length(trimmed_left, trimmed_right, true);
                cis::dna_strand strand = samline->query_on_pos_strand() ? cis::POS : cis::NEG;

                contig_iter = contigs.find(cis::dna_t(species, samline->rname));
                assert(contig_iter != contigs.end());

                cis::region sam_source_region(*contig_iter, std::string(""), 
                                              samline->pos, samline->pos + read_span_on_source, 
                                              strand, -1, 0, 0, NULL);
                
                
                cis::dna_t const* contig = &(*contig_iter);
                cis::region_tree const* tree_on_contig = transcript_extent_trees[contig];
                cis::REGIONS_MULTI overlapping_sources = 
                    cis::IntervalOverlap(*tree_on_contig, sam_source_region);
                
                if (! overlapping_sources.empty())
                {
                    cis::region const* first_overlapping_tx_extent =
                        *(overlapping_sources.begin());
                    min_transcript_proj = tx_extent_to_projection[first_overlapping_tx_extent];
                }

                for (cis::REGIONS_MULTI::const_iterator rit = overlapping_sources.begin();
                     rit != overlapping_sources.end(); ++rit)
                {
                    SP_ITER tg_iter = tx_extent_to_projection[*rit];
                    SequenceProjection const& cur_transcript_g2t = *tg_iter;
                    bool projection_applied = ApplyProjectionToSAM(cur_transcript_g2t, samline);

                    if (projection_applied)
                    {
                        sam_buffer.insert(samline);
                    }
                }
            }
            //purge if necessary
            if (prev_min_transcript_proj != min_transcript_proj)
            {
                //we're on a new transcript
                sam_buffer.update_lowbound(samline);
                sam_buffer.purge(output_sam_fh, NULL, NULL, false);
            }


        }
    }
    sam_buffer.purge(output_sam_fh, NULL, NULL, true);

    fclose(input_sam_fh);
    fclose(output_sam_fh);

    if (! expansion_mode)
    {
        for (cis::TREE_MAP::iterator tit = transcript_extent_trees.begin();
             tit != transcript_extent_trees.end(); ++tit)
        {
            delete (*tit).second;
        }
    }
    
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
}
