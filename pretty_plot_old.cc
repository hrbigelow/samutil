#include "pretty_plot.h"

#include "dep/tools.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <limits.h>
/*
  Input a set of per-locus values, and a set of exon annotations,
  Apply a cigar_string transformation of meta_exons of each gene, 
  producing a set of data easy to plot.
 */


size_t pseudo_intron_length_def = 10;


Cigar::CIGAR_VEC add_pseudo_spacers(Cigar::CIGAR_VEC const& meta_gene_cigar,
                                    size_t pseudo_intron_length)
{
    Cigar::CIGAR_VEC cigar_spaced;

    for (std::vector<Cigar::Unit>::iterator it = meta_gene_cigar.begin();
         it != meta_gene_cigar.end(); ++it)
    {
        if ((*it).op == Cigar::D && it != meta_gene_cigar.begin())
        {
            cigar_spaced.push_back(Cigar::Unit(Cigar::I, pseudo_intron_length));
            if ((*it).length > pseudo_intron_length)
            {
                //(*it).length -= pseudo_intron_length;
            }
        }
        cigar_spaced.push_back(*it);
    }
    return cigar_spaced;
}

METAGENE_CIGARS parse_transformation(char const* metaexon_transformation_file,
                                     size_t pseudo_intron_length)
{
    METAGENE_CIGARS meta_gene_projections;
    //parse coord_transform file
    FILE * genome_to_metaexon_transform_fh = 
        open_if_present(metaexon_transformation_file, "r");

    char contig[1000];
    char gene[1000];
    char cigar_string[10000];
    size_t position;
    char strand;

    while (! feof(genome_to_metaexon_transform_fh))
    {
        fscanf(genome_to_metaexon_transform_fh, 
               "%s\t%s\t%zu\t%c\t%s\n", contig, gene, &position, &strand, cigar_string);

        Cigar::CIGAR_VEC cigar = Cigar::FromString(cigar_string, position);
        Cigar::CIGAR_VEC cigar_spaced;

        for (std::vector<Cigar::Unit>::iterator it = cigar.begin();
             it != cigar.end(); ++it)
        {
            if ((*it).op == Cigar::D && it != cigar.begin())
            {
                cigar_spaced.push_back(Cigar::Unit(Cigar::I, pseudo_intron_length));
                if ((*it).length > pseudo_intron_length)
                {
                    //(*it).length -= pseudo_intron_length;
                }
            }
            cigar_spaced.push_back(*it);
        }
        unique_gene_description gene_desc(contig, gene, strand);
        meta_gene_projections[gene_desc] = cigar_spaced;
    }
    fclose(genome_to_metaexon_transform_fh);

    return meta_gene_projections;
}


cis::TREE_MAP build_transcript_trees(METAGENE_CIGARS const& transcript_extents,
                                     std::set<cis::dna_t> const& contigs)
{
    
    cis::REGIONS_MULTI transcript_extent_regions;
    cis::TREE_MAP transcript_extent_trees;

    METAGENE_CIGARS::const_iterator metagene_iter;
    std::set<cis::dna_t>::const_iterator contig_iter;

    assert(! contigs.empty());
    std::string const& species = (*contigs.begin()).species();

    Cigar::CIGAR_ITER trimmed_left, trimmed_right;
    for (metagene_iter = transcript_extents.begin(); 
         metagene_iter != transcript_extents.end(); ++metagene_iter)
    {
        std::string const& contig_name = metagene_iter->first.contig;
        std::string const& gene_name = metagene_iter->first.gene;
        cis::dna_strand strand = 
        metagene_iter->first.strand == '+' ? cis::POS : cis::NEG;

        contig_iter = contigs.find(cis::dna_t(species, contig_name));
        assert(contig_iter != contigs.end());
        Cigar::CIGAR_VEC const& metagene_cigar = metagene_iter->second;

        int64_t tx_start_on_genome = Cigar::LeftOffset(metagene_cigar, true);
        Cigar::Trim(metagene_cigar, false, &trimmed_left, &trimmed_right);
        int64_t tx_span_on_genome = Cigar::Length(trimmed_left, trimmed_right, true);
        
        //tx_g2t_projection.same_strand ? cis::POS : cis::NEG;
        
        cis::region const* transcript_extent = 
        new cis::region(*contig_iter, gene_name, tx_start_on_genome, 
                        tx_start_on_genome + tx_span_on_genome, 
                        strand, 0, 0, 0, &(*metagene_iter));
        
        transcript_extent_regions.insert(transcript_extent);
    }
    return cis::BuildTrees(transcript_extent_regions);
}


std::set<cis::dna_t> 
build_contigs(METAGENE_CIGARS const& meta_gene_projections, std::string species)
{
    METAGENE_CIGARS::const_iterator mit;
    std::set<cis::dna_t> contigs;
    for (mit = meta_gene_projections.begin(); mit != meta_gene_projections.end(); ++mit)
    {
        contigs.insert(cis::dna_t(species, (*mit).first.contig, 0, INT_MAX));
    }
    return contigs;
}

int pretty_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "pretty_plot\n"
            "Author: Henry Bigelow (hbigelow@amgen.com)\n\n"
            "pretty_plot loci           Plot individual loci in condensed transcript coordinates\n"
            "pretty_plot exons          Plot exons in condensed transcript coordinates\n"
            );
    return 1;
}


int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        return pretty_usage();
    }
    else if (strcmp(argv[1], "loci") == 0)
    {
        return main_pretty_loci(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "exons") == 0)
    {
        return main_pretty_exons(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
    
}
