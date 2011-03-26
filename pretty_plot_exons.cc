#include "pretty_plot.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include "cisortho/nested.h"
#include "cisortho/region.h"
#include "cisortho/dna.h"

int pretty_exons_usage(size_t ldef)
{
    fprintf(stderr,
            "pretty_plot exons [OPTIONS] genome_to_metaexon.txt exons_input.txt exons_output.txt\n\n"
            "Options:\n\n"
            "-l  INT   length for pseudo-introns [%Zu]\n",
            ldef
            );
    return 1;

}


int main_pretty_exons(int argc, char **argv)
{
    size_t pseudo_intron_length = pseudo_intron_length_def;

    char c;
    while ((c = getopt(argc, argv, "l:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'l': pseudo_intron_length = static_cast<size_t>(atof(optarg)); break;
        default: return pretty_exons_usage(pseudo_intron_length_def); break;
        }
    }

    int num_required_args = 3;

    if (argc != num_required_args + optind)
    {
        return pretty_exons_usage(pseudo_intron_length_def);
    }


    char contig_name[1000];
    char gene[1000];
    char transcript[1000];
    char strand_char;

    size_t contig_start, contig_end, transcript_number;

    std::string fake_species(""); //not used

    //INPUT: contig gene transcript strand position cigar_string
    char * genome_to_metaexon_transform = argv[optind];
    char * exon_input_file = argv[optind + 1];
    char * exon_output_file = argv[optind + 2];

    FILE * exon_input_fh = fopen(exon_input_file, "r");
    FILE * exon_output_fh = fopen(exon_output_file, "w");

    METAGENE_CIGARS meta_gene_projections = 
    parse_transformation(genome_to_metaexon_transform, pseudo_intron_length);
    
    typedef std::set<cis::dna_t>::const_iterator DNAS_ITER;
    DNAS_ITER contig_iter;
    std::set<cis::dna_t> contigs = build_contigs(meta_gene_projections, fake_species);
    cis::TREE_MAP meta_gene_trees = 
    build_transcript_trees(meta_gene_projections, contigs);

    while (! feof(exon_input_fh))
    {
        fscanf(exon_input_fh, "%s\t%zu\t%zu\t%c\t%s\t%s\t%zu\n",
               contig_name, &contig_start, &contig_end, &strand_char, gene, transcript,
               &transcript_number);

        cis::dna_strand strand = strand_char == '+' ? cis::POS : cis::NEG;

        contig_iter = contigs.find(cis::dna_t(fake_species, contig_name));
        assert(contig_iter != contigs.end());

        cis::dna_t const* contig = &(*contig_iter);

        cis::region query_meta_gene(*contig, std::string(""), contig_start, 
                                    contig_end, strand, -1, 0, 0, NULL);
        

        cis::region_tree const* tree_on_contig = meta_gene_trees[contig];

        cis::REGIONS_MULTI overlapping_meta_genes = 
        cis::IntervalOverlap(*tree_on_contig, query_meta_gene);

        if (overlapping_meta_genes.empty())
        {
            fprintf(stderr, "Error: Did not find metaexon transform for exon gene entity "
                    "%s\n", contig_name);
            exit(1);
        }

        bool start_missing_projection;
        bool end_missing_projection;

        for (cis::REGIONS_MULTI::const_iterator rit = overlapping_meta_genes.begin();
             rit != overlapping_meta_genes.end(); ++rit)
        {
            METAGENE_CIGARS::value_type const& meta_gene_data =
            * static_cast<METAGENE_CIGARS::value_type const*>((*rit)->payload);

            unique_gene_description const& associated_gene_desc = meta_gene_data.first;
            Cigar::CIGAR_VEC const& genome_to_metaexon = meta_gene_data.second;

            int64_t metaexon_start = 
                Cigar::ProjectCoord(genome_to_metaexon.begin(), genome_to_metaexon.end(),
                                    contig_start, true, & start_missing_projection);

            int64_t metaexon_end = 
                Cigar::ProjectCoord(genome_to_metaexon.begin(), 
                                    genome_to_metaexon.end(),
                                    contig_end, true, & end_missing_projection);

            bool same_gene = strcmp(associated_gene_desc.gene.c_str(), gene) == 0;

            //since some of these may be
            if (same_gene && (start_missing_projection || end_missing_projection))
            {
                //the projected coordinates on the actual gene should NOT
                //be missing.
                fprintf(stderr, "Error: %s %Zu %Zu %s %s %Zu has coordinates "
                        "out of bounds of meta-exons\n",
                        contig_name, contig_start, contig_end, gene, transcript, transcript_number);
                exit(1);
            }

            // display any non-existing coordinates in the middle of meta-introns
            size_t start_adjust = start_missing_projection ? (pseudo_intron_length - 2) : 0;
            size_t end_adjust = end_missing_projection ? 2 : 0;
            
            fprintf(exon_output_fh, "%s\t%s\t%s\t%Zu\t%Zi\t%Zi\ts: %Zi\te: %Zi\n",
                    associated_gene_desc.gene.c_str(), gene, transcript, 
                    transcript_number, metaexon_start - start_adjust, 
                    metaexon_end - end_adjust, start_adjust, end_adjust);
        }
        
    }
    fclose(exon_input_fh);
    fclose(exon_output_fh);
    return 0;
}
