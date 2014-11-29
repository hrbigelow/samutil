#include "readsim_aux.h"

#include <cstring>
#include <fstream>
#include <getopt.h>
#include <sys/timeb.h>

#include "cisortho/memory.h"
#include "dep/tools.h"
#include "cisortho/dnacol.h"
#include "cisortho/dna.h"
#include "sam_line.h"
#include "sam_buffer.h"

int get_tx_sequence_usage(char const* ddef, size_t ldef)
{
    fprintf(stderr, 
            "\nUsage:\n\n"
            "samutil get_tx_sequence [OPTIONS] species transcripts.gtf genome.dnas transcripts.fa \n"
            "\n"
            "Options:\n\n"
            "-d  STRING   dna directory for finding pieces in genome.dnas file ['%s']\n"
            "-l  INT      maximum line length for fasta output (0 for unlimited) ['%Zu']\n"
            "-r  FLAG     reverse complement negative-stranded transcripts\n\n"
            "required file genome.dnas is produced from make_dnas_file and fasta2cisfasta\n"
            ,
            ddef, ldef);

    return 1;
}


int main_get_tx_sequence(int argc, char ** argv)
{
    char c;
    char dna_directory_def[] = ".";
    char dna_directory[1000];
    strcpy(dna_directory, dna_directory_def);

    size_t line_length_def = 0;
    size_t line_length = line_length_def;

    bool revcomp_neg_strand_tx = false;
    
    while ((c = getopt(argc, argv, "d:l:r")) >= 0)
    {
        switch(c)
        {
        case 'd': strcpy(dna_directory, optarg); break;
        case 'l': line_length = static_cast<size_t>(atoi(optarg)); break;
        case 'r': revcomp_neg_strand_tx = true; break;
        default: return get_tx_sequence_usage(dna_directory_def, line_length_def); break;
        }
    }

    int min_req_args = 4;
    int min_arg_count = optind + min_req_args;

    if ((argc != min_arg_count) && argc != (min_arg_count + 1))
    {
        return get_tx_sequence_usage(dna_directory_def, line_length_def);
    }

    char * species = argv[optind];
    char * gtf_file = argv[optind + 1];
    char * dna_index_file = argv[optind + 2];
    char * output_fasta_file = argv[optind + 3];

    FILE * output_fasta_fh = open_if_present(output_fasta_file, "w");

    std::set<SequenceProjection> genome_to_tx = 
        gtf_to_sequence_projection(gtf_file, species);

    //parse dna file
    cis::dna_collection dnacollection;
    dnacollection.make_ready(dna_index_file, std::string(dna_directory));

    //CIGAR to describe the relationship between genome dna and matepair
    std::set<SequenceProjection>::iterator sit;
    LOCUS_SET::const_iterator dummy_iter;
    Cigar::CIGAR_ITER low_trimmed, high_trimmed;

    const bool zero_terminate = true;

    for (sit = genome_to_tx.begin(); sit != genome_to_tx.end(); ++sit)
    {
        SequenceProjection const& sp = *sit;
        cis::dna_t const* dna = 
            GetDNAByName(dnacollection, sp.species, sp.source_dna);

        if (dna == NULL)
        {
            fprintf(stderr, "Skipping contig %s missing from sequence input\n",
                    sp.source_dna.c_str());
            continue;
        }

        size_t transcript_length = Cigar::Length(sp.cigar, false);

        char * transcript_sequence = new char[transcript_length + 1];
        Cigar::Trim(sp.cigar, true, &low_trimmed, &high_trimmed);

        get_transformed_sequence(dna, sp.cigar,
                                 dummy_iter, dummy_iter,
                                 transcript_sequence,
                                 zero_terminate);

        if (! sp.same_strand && revcomp_neg_strand_tx)
        {
            cis::ReverseComplement(transcript_sequence, transcript_length,
                                   transcript_sequence);
        }

        fprintf(output_fasta_fh, ">%s\n", sp.target_dna.c_str());
        if (line_length > 0 && line_length < transcript_length)
        {
            char * sub_sequence_buf = new char[line_length + 1];
            sub_sequence_buf[line_length] = '\0';
            for (size_t offset = 0; offset < transcript_length; offset += line_length)
            {
                strncpy(sub_sequence_buf, transcript_sequence + offset, line_length);
                fprintf(output_fasta_fh, "%s\n", sub_sequence_buf);
            }
            delete sub_sequence_buf;
        }
        else
        {
            fprintf(output_fasta_fh, "%s\n", transcript_sequence);
        }

        delete transcript_sequence;
    }
    fclose(output_fasta_fh);
    return 0;
           
}
