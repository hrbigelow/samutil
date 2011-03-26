#include "pretty_plot.h"

#include <cstdio>
#include <cstdlib>

int pretty_loci_usage(size_t ldef)
{
    fprintf(stderr,
            "pretty_plot loci [OPTIONS] genome_to_metaexon.txt locus_input.txt locus_output.txt\n\n"
            "Options:\n\n"
            "-l  INT   length for pseudo-introns [%Zu]\n",
            ldef
            );
    return 1;
}


int main_pretty_loci(int argc, char **argv)
{

    size_t pseudo_intron_length = pseudo_intron_length_def;

    char c;
    while ((c = getopt(argc, argv, "l:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'l': pseudo_intron_length = static_cast<size_t>(atof(optarg)); break;
        default: return pretty_loci_usage(pseudo_intron_length_def); break;
        }
    }

    int num_required_args = 3;

    if (argc != num_required_args + optind)
    {
        return pretty_loci_usage(pseudo_intron_length_def);
    }

    //INPUT: contig contig_start strand guide_depth correct_depth error_depth gene
    char * genome_to_metaexon_transform = argv[optind];
    char * locus_input_file = argv[optind + 1];
    char * locus_output_file = argv[optind + 2];

    FILE * locus_input_fh = fopen(locus_input_file, "r");
    FILE * locus_output_fh = fopen(locus_output_file, "w");

    size_t contig_start;
    size_t guide_depth, correct_depth, error_depth;

    char contig[1000];
    char gene[1000];
    char strand;

    std::map<unique_gene_description, Cigar::CIGAR_VEC> 
        meta_gene_projections = parse_transformation(genome_to_metaexon_transform,
                                                     pseudo_intron_length);

    bool is_missing_projection;

    while (! feof(locus_input_fh))
    {
        fscanf(locus_input_fh, "%s\t%zu\t%c\t%zu\t%zu\t%zu\t%s\n",
               contig, &contig_start, &strand, &guide_depth, &correct_depth, &error_depth, gene);
        
        unique_gene_description gene_desc(contig, gene, strand);
        if (meta_gene_projections.find(gene_desc) == meta_gene_projections.end())
        {
            // fprintf(skipped_loci_fh, "Missing_projection\t%s\t%Zu\t%c\t%Zu\t%Zu\t%Zu\t%s\n",
            //         contig, contig_start, strand, guide_depth, 
            //         correct_depth, error_depth, gene);
            continue;
        }

        Cigar::CIGAR_VEC const& genome_to_metaexon = meta_gene_projections[gene_desc];

        int64_t metaexon_start = 
            Cigar::ProjectCoord(genome_to_metaexon.begin(), 
                                genome_to_metaexon.end(),
                                contig_start, true, & is_missing_projection);

        if (is_missing_projection)
        {
            // fprintf(skipped_loci_fh, "Intronic\t%s\t%Zu\t%c\t%Zu\t%Zu\t%Zu\t%s\n",
            //         contig, contig_start, strand, guide_depth, 
            //         correct_depth, error_depth, gene);
        }
        else
        {
            fprintf(locus_output_fh, "%s\t%s\t%c\t%Zu\t%Zu\t%Zu\t%Zu\n",
                    contig, gene, strand, 
                    metaexon_start, guide_depth, correct_depth, error_depth);
        }
    }
    fclose(locus_input_fh);
    fclose(locus_output_fh);

    return 0;
}
