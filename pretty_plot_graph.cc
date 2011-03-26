//this must come first for some reason
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <map>
#include <set>

#include "pretty_plot.h"
#include "nclist.h"
#include "readsim_aux.h"
#include "cigar_ops.h"
#include "meta_gene.h"

int pretty_plot_graph_usage(size_t ldef)
{
    fprintf(stderr,
            "\npretty_plot graph [OPTIONS] exons.gtf loci.txt output_exons.txt output_loci.txt\n\n"
            "Options:\n\n"
            "-l  INT   length for pseudo-introns [%Zu]\n\n",
            ldef
            );
    return 1;

}

int main_pretty_plot_graph(int argc, char **argv)
{
    size_t pseudo_intron_length = pseudo_intron_length_def;

    char c;
    while ((c = getopt(argc, argv, "l:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'l': pseudo_intron_length = static_cast<size_t>(atof(optarg)); break;
        default: return pretty_plot_graph_usage(pseudo_intron_length_def); break;
        }
    }

    int num_required_args = 4;

    if (argc != num_required_args + optind)
    {
        return pretty_plot_graph_usage(pseudo_intron_length_def);
    }


    //INPUT: contig gene transcript strand position cigar_string
    char * exon_input_gtf_file = argv[optind];
    char * loci_input_file = argv[optind + 1];
    char * exon_output_file = argv[optind + 2];
    char * loci_output_file = argv[optind + 3];

    FILE * exon_input_gtf_fh = open_if_present(exon_input_gtf_file, "r");
    FILE * exon_output_fh = open_if_present(exon_output_file, "w");

    FILE * loci_input_fh = open_if_present(loci_input_file, "r");
    FILE * loci_output_fh = open_if_present(loci_output_file, "w");

    MetaGene meta_gene;
    meta_gene.Initialize(exon_input_gtf_fh, pseudo_intron_length);
    //print_tree(meta_gene.gene_extent_tree, 0);

    rewind(exon_input_gtf_fh);

    std::vector<IntervalTree const*> overlaps;
    std::vector<IntervalTree const*>::iterator ovit;
    GTFEntry gtf_entry;

    size_t gtf_start_bound, gtf_end_bound;
    size_t gene_projected_bound;

    bool add_padding = false;
    bool inserts_are_introns = true;

    bool projection_applied;
                

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

        Cigar::CIGAR_VEC gtf_cigar(1, Cigar::Unit(Cigar::D, gtf_entry.start_bound()));
        gtf_cigar.push_back(Cigar::Unit(Cigar::M, gtf_entry.length()));

        Cigar::CIGAR_VEC gtf_projected =
            Cigar::TransitiveMerge(meta_gene.meta_exon[gtf_entry.seqname],
                                   meta_gene.meta_exon_offsets[gtf_entry.seqname],
                                   gtf_cigar, add_padding, inserts_are_introns);
        
        gtf_start_bound = Cigar::LeftOffset(gtf_projected, true);
        gtf_end_bound = gtf_start_bound + Cigar::Length(gtf_projected, false);

        for (ovit = overlaps.begin(); ovit != overlaps.end(); ++ovit)
        {
            IntervalTree const* overlap = *ovit;
            assert(overlap->has_interval);

            GENE_ITERS * genes = static_cast<GENE_ITERS *>(overlap->payload);

            for (size_t i = 0; i != (*genes).size(); ++i)
            {
                GENE_INDEX::value_type & gene_data = *(*genes)[i];

                unique_gene_description const& gene_desc = gene_data.first;
                BoundsInfo const& gene_bounds = gene_data.second;

                gene_projected_bound = 
                    Cigar::ProjectCoord(meta_gene.meta_exon[gtf_entry.seqname],
                                        meta_gene.meta_exon_offsets[gtf_entry.seqname],
                                        gene_bounds.start_boundary,
                                        & projection_applied);

                assert(projection_applied);

                char const* containing_gene = gene_desc.gene.c_str();

                size_t transcript_num = meta_gene.transcript_number(gtf_entry);

                //containing_gene, gene, transcript, transcript_num, start, end
                fprintf(exon_output_fh, "%s\t%s\t%s\t%Zu\t%s\t%Zu\t%Zu\n",
                        containing_gene, gtf_entry.gene_id, gtf_entry.transcript_id,
                        transcript_num, gene_desc.contig.c_str(), 
                        gtf_start_bound - gene_projected_bound + 1, 
                        gtf_end_bound - gene_projected_bound);
            
            }
        }
    }

    fclose(exon_input_gtf_fh);
    fclose(exon_output_fh);

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

    //for representing a single-base locus ;)
    Cigar::CIGAR_VEC unit_match(1, Cigar::Unit(Cigar::D, 0));
    unit_match.push_back(Cigar::Unit(Cigar::M, 1));
    
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
            //not an error.  this locus falls outside of any annotated exon
            continue;
        }


        //projects the locus into the collapsed meta-exon structure of
        //each containing gene
        else
        {
            size_t start_bound, end_bound;

            unit_match[0].length = locus_start_boundary;
            Cigar::CIGAR_VEC locus_cigar =
                Cigar::TransitiveMerge(meta_gene.meta_exon[locus_contig],
                                       meta_gene.meta_exon_offsets[locus_contig],
                                       unit_match, add_padding, inserts_are_introns);

            start_bound = Cigar::LeftOffset(locus_cigar, true);
            end_bound = start_bound + Cigar::Length(locus_cigar, false);
            
            if (start_bound == end_bound)
            {
                //this locus falls in a deleted region
                continue;
            }
            
            for (ovit = overlaps.begin(); ovit != overlaps.end(); ++ovit)
            {
                IntervalTree const* overlap = *ovit;

                GENE_ITERS * genes = static_cast<GENE_ITERS *>(overlap->payload);


                for (size_t i = 0; i != (*genes).size(); ++i)
                {
                    GENE_INDEX::value_type & gene_data = *(*genes)[i];
                    unique_gene_description const& gene_desc = gene_data.first;
                    BoundsInfo const& gene_bounds = gene_data.second;

                    gene_projected_bound = 
                        Cigar::ProjectCoord(meta_gene.meta_exon[gtf_entry.seqname],
                                            meta_gene.meta_exon_offsets[gtf_entry.seqname],
                                            gene_bounds.start_boundary,
                                            & projection_applied);

                    assert(projection_applied);
                    
                    char const* containing_gene = gene_desc.gene.c_str();

                    fprintf(loci_output_fh, "%s\t%s\t%Zu\t%s\n",
                            containing_gene, locus_contig, 
                            start_bound - gene_projected_bound + 1, 
                            extra_info);
                }
            }
        }
    }
    fclose(loci_input_fh);
    fclose(loci_output_fh);

    delete extra_info;

    return 0;
}
