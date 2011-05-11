//this must come first for some reason

#include <map>
#include <set>


#include "pretty_plot.h"
#include "nclist.h"
#include "readsim_aux.h"
#include "cigar_ops.h"
#include "meta_gene.h"

#include "dep/tools.h"

int pretty_plot_gtf_sam_usage(size_t ldef)
{
    fprintf(stderr,
            "\npretty_plot gtf_sam [OPTIONS] reference.gtf input.gtf output.gtf input.sam output.sam\n\n"
            "Options:\n\n"
            "-l  INT   length for pseudo-introns [%Zu]\n\n",
            ldef
            );

    fprintf(stderr,
            "Use reference.gtf to define the coordinate transformation that collapses empty space (introns).\n"
            "Using this coordinate transform, collapse input.gtf and input.sam, outputting\n"
            "output.gtf and output.sam\n"
            "Commonly, input.gtf will be the same as reference.gtf, but need not be.\n\n");
    return 1;

}

int main_pretty_plot_gtf_sam(int argc, char **argv)
{
    size_t pseudo_intron_length = pseudo_intron_length_def;

    char c;
    while ((c = getopt(argc, argv, "l:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'l': pseudo_intron_length = static_cast<size_t>(atof(optarg)); break;
        default: return pretty_plot_gtf_sam_usage(pseudo_intron_length_def); break;
        }
    }

    int num_required_args = 5;

    if (argc != num_required_args + optind)
    {
        return pretty_plot_gtf_sam_usage(pseudo_intron_length_def);
    }


    //INPUT: contig gene transcript strand position cigar_string
    char * gtf_reference_file = argv[optind];
    char * gtf_input_file = argv[optind + 1];
    char * gtf_output_file = argv[optind + 2];

    char * sam_input_file = argv[optind + 3];
    char * sam_output_file = argv[optind + 4];

    FILE * gtf_reference_fh = open_or_die(gtf_reference_file, "r", "Reference gtf file");

    FILE * gtf_input_fh = open_if_present(gtf_input_file, "r");
    FILE * gtf_output_fh = open_if_present(gtf_output_file, "w");

    FILE * sam_input_fh = open_if_present(sam_input_file, "r");
    FILE * sam_output_fh = open_if_present(sam_output_file, "w");

    //Process the jumps, creating CIGAR_VEC and gene_bounds structures
    MetaGene meta_gene;
    meta_gene.Initialize(gtf_reference_fh, pseudo_intron_length);

    fclose(gtf_reference_fh);

    GTFEntry gtf_entry;

    bool add_padding = false;
    bool inserts_are_introns = true;

    //parse and transform GTF
    while (gtf_input_fh != NULL && gtf_entry.get_next_record(gtf_input_fh))
    {
        if (! gtf_entry.is_data_line)
        {
            continue;
        }
        assert(meta_gene.meta_exon.find(gtf_entry.seqname) !=
               meta_gene.meta_exon.end());

        Cigar::CIGAR_VEC gtf_cigar(1, Cigar::Unit(Cigar::D, gtf_entry.start_bound()));
        gtf_cigar.push_back(Cigar::Unit(Cigar::M, gtf_entry.length()));

        Cigar::CIGAR_VEC projected =
            Cigar::TransitiveMerge(meta_gene.meta_exon[gtf_entry.seqname],
                                   meta_gene.meta_exon_offsets[gtf_entry.seqname],
                                   gtf_cigar, add_padding, inserts_are_introns);
            
        gtf_entry.start = Cigar::LeftOffset(projected, true) + 1;
        gtf_entry.end = gtf_entry.start_bound() + Cigar::Length(projected, false);

        gtf_entry.print(gtf_output_fh);
    }

    close_if_present(gtf_input_fh);
    close_if_present(gtf_output_fh);

    //for SAM file
    bool allow_absent_seq_qual = false; //this is false because

    std::string transformed_cigar;
    SamLine * samline;

    //parse and transform SAM file
    while (sam_input_fh != NULL && ! feof(sam_input_fh))
    {

        samline = new SamLine(sam_input_fh, allow_absent_seq_qual);
        switch (samline->parse_flag)
        {
        case END_OF_FILE: 
            delete samline;
            break;

        case PARSE_ERROR:
            fprintf(stderr, "Parse error in input sam file %s", sam_input_file);
            exit(1);
            break;
        case HEADER:
            samline->print(sam_output_fh, true);
            delete samline;
            break;
        case DATA_LINE:
            std::map<std::string, Cigar::CIGAR_VEC>::iterator mit =
                meta_gene.meta_exon.find(samline->rname);

            if (mit == meta_gene.meta_exon.end())
            {
                fprintf(stderr, "Error: Found data line in SAM on contig %s"
                        " which is not found in GTF file\n", samline->rname);
                exit(1);
            }

            Cigar::CIGAR_VEC sam_cigar = Cigar::FromString(samline->cigar, samline->zero_based_pos());
            Cigar::CIGAR_VEC sam_cigar_projected = 
                Cigar::TransitiveMerge((*mit).second,
                                       meta_gene.meta_exon_offsets[samline->rname],
                                       sam_cigar, add_padding, inserts_are_introns);

            transformed_cigar = Cigar::ToString(Cigar::Trim(sam_cigar_projected, false));
            char * tmp_cigar = new char[transformed_cigar.size() + 1];
            strcpy(tmp_cigar, transformed_cigar.c_str());

            samline->pos = Cigar::LeftOffset(sam_cigar_projected, true);
            samline->cigar = tmp_cigar;
            samline->print(sam_output_fh, false);
            delete samline;
            delete tmp_cigar;

        }
    }

    close_if_present(sam_input_fh);
    close_if_present(sam_output_fh);

    return 0;
}
