//Generate a projection-order transcript header from GTF and genome header.

#include "sam_order.h"
#include "cigar_ops.h"
#include "readsim_aux.h"
#include "gtf.h"
#include "dep/tools.h"

#include <set>

int generate_projection_header_usage()
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil gen_header genome.header.sam transcripts.gtf transcriptome.header.sam\n\n"
            "Generate a valid SAM header of transcript contigs in order consistent with\n"
            "their projections to the genome as specified in genome.header.sam\n"
            );
    return 1;
}

int main_generate_projection_header(int argc, char ** argv)
{
    if (argc != 4)
    {
        return generate_projection_header_usage();
    }

    char * genome_header_file = argv[1];
    char * gtf_file = argv[2];
    char * output_header_file = argv[3];

    FILE * genome_header_fh = open_or_die(genome_header_file, "r", "Input genome SAM header file");
    FILE * gtf_fh = open_or_die(gtf_file, "r", "Input GTF file");
    FILE * output_header_fh = open_or_die(output_header_file, "w", "Output genome-projected SAM header file");

    SamOrder genome_order(SAM_POSITION_RID, "ALIGN");
    genome_order.InitFromFile(genome_header_fh);
    genome_order.AddHeaderContigStats(genome_header_fh);

    char const* species = "empty";


    std::set<SequenceProjection> genome_to_tx = 
        gtf_to_sequence_projection(gtf_file, species);

    //needed if we are expanding reads that are originally mapped to a
    //transcript to genomic coordinates
    std::set<SequenceProjection> tx_to_genome;

    // 1. tx_to_genome  (projections of transcripts to the genome)
    std::set<SequenceProjection>::const_iterator gt_iter;
    for (gt_iter = genome_to_tx.begin(); gt_iter != genome_to_tx.end(); ++gt_iter)
    {
        tx_to_genome.insert(tx_to_genome.begin(), 
                            InvertProjection((*gt_iter)));
    }

    GenerateProjectionHeader(genome_order.contig_offsets,
                             tx_to_genome,
                             output_header_fh);

    fclose(genome_header_fh);
    fclose(gtf_fh);
    fclose(output_header_fh);
    
    return 0;
}
