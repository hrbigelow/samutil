/*
  Preparation of Auxiliary data structures:

  1. parse GTF file, creating:
  std::map<contig, std::set<size_t, annot_feature_jump> > jumps

  2. process the jumps, creating, simultaneously:
  CIGAR_VEC meta_exons
  std::map<gene, std::pair<start_chunk, end_chunk> > gene_cigar_bounds
  std::map<gene, std::pair<start_chunk, end_chunk> > gene_coord_bounds

  3. from the gene_bounds and meta_exons, create:
  TREE_MAP gene_extents

  Main loop:

  1. parse GTF file:
  for each exon record:
  a. find containing gene(s) by querying gene_extents
  b. for each gene:
  query/store transcript, get transcript_number.
  i. calculate adjusted coords using CIGAR_VEC[gene_bounds[gene].first] ...
  ii. output adjusted, annotated exon record

  2. parse loci file:
  for each locus:
  a. find containing gene(s) by querying gene_extents
  b. for each gene:
  i. calculate adjusted coordinates
  ii. output adjusted locus, with found gene annotation

  
*/

//this must come first for some reason
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <map>
#include <set>



#include "nclist.h"
#include "readsim_aux.h"
#include "cigar_ops.h"


int pretty_plot_annot_usage(size_t ldef)
{
    fprintf(stderr,
            "pretty_plot_annot [OPTIONS] exons.gtf loci.txt output_exons.txt output_loci.txt\n\n"
            "Options:\n\n"
            "-l  INT   length for pseudo-introns [%Zu]\n",
            ldef
            );
    return 1;

}


struct less_char_ptr
{
    bool operator()(char const* a, char const* b) const
    {
        return strcmp(a, b) < 0;
    }
};


typedef std::set<std::string> CHAR_SET;
typedef std::map<std::string, CHAR_SET> GENE_TRANSCRIPT_MAP;


size_t TranscriptNumber(GTFEntry const& gtf_entry, 
                        GENE_TRANSCRIPT_MAP const& transcripts)
{
    CHAR_SET::iterator transcript_iter;
    GENE_TRANSCRIPT_MAP::const_iterator transcripts_iter;

    transcripts_iter = transcripts.find(gtf_entry.gene_id);
    assert(transcripts_iter != transcripts.end());
   
    transcript_iter = (*transcripts_iter).second.find(gtf_entry.transcript_id);
    assert(transcript_iter != (*transcripts_iter).second.end());
    
    return std::distance(transcript_iter, (*transcripts_iter).second.end());
}
    







size_t pseudo_intron_length_def = 20;


int main(int argc, char **argv)
{
    size_t pseudo_intron_length = pseudo_intron_length_def;

    char c;
    while ((c = getopt(argc, argv, "l:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'l': pseudo_intron_length = static_cast<size_t>(atof(optarg)); break;
        default: return pretty_plot_annot_usage(pseudo_intron_length_def); break;
        }
    }

    int num_required_args = 4;

    if (argc != num_required_args + optind)
    {
        return pretty_plot_annot_usage(pseudo_intron_length_def);
    }


    //INPUT: contig gene transcript strand position cigar_string
    char * exon_input_gtf_file = argv[optind];
    char * loci_input_file = argv[optind + 1];
    char * exon_output_file = argv[optind + 2];
    char * loci_output_file = argv[optind + 3];

    FILE * exon_input_gtf_fh = fopen(exon_input_gtf_file, "r");
    FILE * exon_output_fh = fopen(exon_output_file, "w");

    FILE * loci_input_fh = fopen(loci_input_file, "r");
    FILE * loci_output_fh = fopen(loci_output_file, "w");



    //Process the jumps, creating CIGAR_VEC and gene_bounds structures
    int64_t cumul_height = 0;
    size_t start = 0;
    size_t end = 0;

    MetaGene meta_gene;
    meta_gene.InitializeBounds(exon_input_gtf_fh, pseudo_intron_length);
    meta_gene.InitializeTree();
    //print_tree(meta_gene.gene_extent_tree, 0);

    rewind(exon_input_gtf_fh);

    //Main loop
    /*
      1. parse GTF file:
      for each exon record:
      a. find containing gene(s) by querying the gene_extent_tree
      b. for each gene:
         i. retrieve transcript number from 'transcripts'
         ii. calculate adjusted coords using CIGAR_VEC[gene_bounds[gene].first] ...
         iii. output adjusted coordinate exon record, annotated with transcript number
              and containing gene (in addition to gene and transcript of origin)
    */
    std::vector<IntervalTree const*> overlaps;
    std::vector<IntervalTree const*>::iterator extent_iter;
    GTFEntry gtf_entry;

    if (do_plot)
    {
        
    }
    else
    {

        //parse and transform GTF
        while (gtf_entry.get_next_record(exon_input_gtf_fh))
        {
            assert(meta_gene.exon_offsets.find(gtf_entry.seqname) !=
                   meta_gene.exon_offsets.end());

            Cigar::CIGAR_VEC gtf_cigar(Cigar::Unit(Cigar::M, gtf_entry.end_boundary() -
                                                   gtf_entry.start_boundary()));

            Cigar::CIGAR_VEC projected =
                Cigar::TransitiveMerge(meta_gene.meta_exon[gtf_entry.seqname],
                                       meta_gene.meta_exon_offsets[gtf_entry.seqname],
                                       gtf_cigar, gtf_entry.start_boundary(), false);
            
            gtf_entry.start = Cigar::LeftOffset(projected, true) + 1;
            gtf_entry.end = gtf_entry.start_boundary() + Cigar::Length(projected, false);

            gtf_entry.print(exon_output_fh);
        }

        fclose(exon_input_gtf_fh);
        fclose(exon_output_fh);

        //for SAM file
        bool ones_based_pos = true;
        std::string transformed_cigar;

        //parse and transform SAM file
        while (! feof(input_sam_fh))
        {

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
                assert(meta_gene.exon_offsets.find(samline->rname) !=
                       meta_gene.exon_offsets.end());

                Cigar::CIGAR_VEC sam_cigar = Cigar::FromString(samline->cigar);
                Cigar::CIGAR_VEC sam_cigar_projected = 
                    Cigar::TransitiveMerge(meta_gene.meta_exon[samline->rname],
                                           meta_gene.meta_exon_offsets[samline->rname],
                                           sam_cigar, samline->pos, false);

                transformed_cigar = Cigar::ToString(sam_cigar_projected);
                if (strlen(samline->cigar) < transformed_cigar.size() + 1)
                {
                    delete samline->cigar;
                    samline->cigar = new char[transformed_cigar.size() + 1];
                }
                strcpy(samline->cigar, transformed_cigar.c_str());

                samline->print(output_sam_fh, ones_based_pos, false);
            }
        }

        fclose(input_sam_fh);
        fclose(output_sam_fh);

    }

    while (gtf_entry.get_next_record(exon_input_gtf_fh))
    {
        
        if (strcmp(gtf_entry.feature, "exon") != 0)
        {
            continue;
        }



        //1. find all overlapping gene intervals on current contig.
        overlaps = meta_gene.gene_extent_tree->overlap(gtf_entry.seqname,
                                                       gtf_entry.start_bound(),
                                                       gtf_entry.end_bound());

        //by definition, every annotated exon must reside in at least one gene extent
        assert(! overlaps.empty());

        size_t start_bound, end_bound;
        bool start_is_missing, end_is_missing;

        //for each gene overlapping this exon, project it into the gene's
        //collapsed coordinate meta-exon structure.
        /* Steps involved:
           1. find overlapping gene extents
           2. for each containing gene, find *contained* genes
           3. for each contained gene, find its meta-exon structure
         */
        for (extent_iter = overlaps.begin(); extent_iter != overlaps.end(); ++extent_iter)
        {

            IntervalTree const* overlap = *extent_iter;
            assert(overlap->has_interval);

            GENE_ITERS * genes = static_cast<GENE_ITERS *>(overlap->payload);

            for (size_t i = 0; i != (*genes).size(); ++i)
            {
                GENE_INDEX::value_type & gene_data = *(*genes)[i];

                unique_gene_description const& gene_desc = gene_data.first;
                BoundsInfo const& gene_bounds = gene_data.second;

                assert(meta_exons.find(gene_desc.contig) != meta_exons.end());

                std::map<std::string, Cigar::CIGAR_VEC>::const_iterator mit;

                mit = meta_exons.find(gene_desc.contig);
                Cigar::CIGAR_VEC const& contig_cigar = (*mit).second;

                size_t start_index = gene_bounds.start_cigar_index;
                assert(start_index < contig_cigar.size());

                Cigar::CIGAR_ITER cigar_start_iter = contig_cigar.begin() + start_index;
                assert((*cigar_start_iter).op == Cigar::M);
            
                start_bound = 
                    Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                                        gtf_entry.start_bound() - gene_bounds.start_boundary, 
                                        true, &start_is_missing);
            
                end_bound = 
                    Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                                        gtf_entry.end_bound() - gene_bounds.start_boundary, 
                                        true, &end_is_missing);
            
                assert(! start_is_missing);
                assert(! end_is_missing);
             
                char const* containing_gene = gene_desc.gene.c_str();

                size_t transcript_num = TranscriptNumber(gtf_entry, transcripts);

                //containing_gene, gene, transcript, transcript_num, start, end
                fprintf(exon_output_fh, "%s\t%s\t%s\t%Zu\t%s\t%Zu\t%Zu\n",
                        containing_gene, gtf_entry.gene_id, gtf_entry.transcript_id,
                        transcript_num, gene_desc.contig.c_str(), 
                        start_bound + 1, end_bound);
            
            }
        }
    }

    fclose(exon_output_fh);
    fclose(exon_input_gtf_fh);

    /*
      2. parse loci file:
      for each locus:
        a. find containing gene(s) by querying gene_extents
        b. for each gene:
           i. calculate adjusted coordinates
           ii. output adjusted locus, with found gene annotation
    */

    char locus_contig[1000];
    size_t locus_start;
    size_t const extra_info_chars = 1000000;
    char * extra_info = new char[extra_info_chars];
    
    while (! feof(loci_input_fh))
    {
        int parsed_fields =
            fscanf(loci_input_fh, "%s\t%zu\t%[^\n]\n", locus_contig, &locus_start, extra_info);

        if (strlen(extra_info) > extra_info_chars)
        {
            fprintf(stderr, "Error: loci extra info field exceeds "
                    "max space of %Zu characters\n", 
                    extra_info_chars);
            exit(1);
        }

        if (parsed_fields != 3)
        {
            fprintf(stderr, "Error: loci input file bad format\n");
            exit(1);
        }

        size_t locus_start_boundary = locus_start - 1;
        size_t locus_end_boundary = locus_start_boundary + 1;

        //1. find all overlapping gene intervals on current contig.
        overlaps = meta_gene.gene_extent_tree->overlap(locus_contig, 
                                                       locus_start_boundary, 
                                                       locus_end_boundary);

        if (overlaps.empty())
        {
            //this locus falls outside of any annotated exon
            continue;
        }


        //projects the locus into the collapsed meta-exon structure of each
        //containing gene
        else
        {
            for (extent_iter = overlaps.begin(); extent_iter != overlaps.end(); ++extent_iter)
            {
                IntervalTree const* overlap = *extent_iter;

                GENE_ITERS * genes = static_cast<GENE_ITERS *>(overlap->payload);

                size_t start_bound, end_bound;
                bool start_is_missing, end_is_missing;

                for (size_t i = 0; i != (*genes).size(); ++i)
                {
                    GENE_INDEX::value_type & gene_data = *(*genes)[i];
                    unique_gene_description const& gene_desc = gene_data.first;
                    BoundsInfo const& gene_bounds = gene_data.second;

                    assert(meta_exons.find(gene_desc.contig) != meta_exons.end());
                    Cigar::CIGAR_VEC const& contig_cigar = meta_exons[gene_desc.contig];

                    size_t start_index = gene_bounds.start_cigar_index;
                    Cigar::CIGAR_ITER cigar_start_iter = contig_cigar.begin() + start_index;
                    assert((*cigar_start_iter).op == Cigar::M);

                    start_bound = 
                        Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                                            locus_start_boundary - gene_bounds.start_boundary, 
                                            true, &start_is_missing);

                    end_bound = 
                        Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                                            locus_end_boundary - gene_bounds.start_boundary, 
                                            true, &end_is_missing);

                    if (start_is_missing || end_is_missing)
                    {
                        //this locus falls outside a meta-exon.  ignore.
                        continue;
                    }
             
                    char const* containing_gene = gene_desc.gene.c_str();

                    fprintf(loci_output_fh, "%s\t%s\t%Zu\t%s\n",
                            containing_gene, locus_contig, start_bound, extra_info);
                }
            }
        }
    }
    fclose(loci_input_fh);
    fclose(loci_output_fh);

    delete extra_info;

    return 0;
}
