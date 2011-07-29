#include "readsim_aux.h"

#include <cstring>
#include <fstream>
#include <getopt.h>

#include "cisortho/memory.h"
#include "dep/tools.h"
#include "cisortho/dnacol.h"
#include "cisortho/dna.h"
#include "sam_helper.h"
#include "sam_buffer.h"
#include "fragment_generator.h"

/*
  Quality code definition from Wikipedia

  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)

 */

int sim_reads_usage(char const* ddef, size_t ldef, 
                    char const* sdef,
                    char cdef, size_t udef,
                    char const* mdef, size_t Mdef, size_t rdef)
{
    fprintf(stderr, 
            "\nUsage:\n\n"
            "sim reads [OPTIONS] genome.dnas transcripts.gtf species transcripts.expr sim.sam\n"
            "    sim.frag sim.read1.fastq [sim.read2.fastq]\n\n"
            "Options:\n\n"
            "-d  STRING   dna directory for finding pieces in genome_dna.fa file [%s]\n"
            "-l  INT      output read length [%Zu]\n"
            "-p  FLAG     produce paired-end reads (if absent, single-end produced)\n"
            "-b  FLAG     if absent, qname format id:chr:strand:start:CIGAR[,chr:strand:start:CIGAR]+\n"
            "-h  FLAG     sample the sense strand of the transcript\n"
            "-q  FLAG     create pairs on same strand (if absent, pairs are on opposite strands)\n"
            "-s  STRING   sampling scheme (see below) [%s]\n"
            "-f  STRING   fastq file 1, quality strings to randomly sample for error model\n"
            "-g  STRING   fastq file 2 (if paired), together with fastq file 1\n"
            "-i  FLAG     output only one of alignment-identical simulated reads\n"
            "-c  CHAR     quality code in fastq input defining quality score of zero [%c]\n"
            "-u  INT      minimum median Phred quality score to use quality string [%Zu]\n"
            "-m  STRING   somatic mutation file: species, contig, position, orig_base, mut_base [%s]\n"
            "-M  INT      number bytes maximum memory to use [%Zu]\n"
            "-r  INT      random seed for repeatability [%Zu]\n\n",
            ddef, ldef, sdef, cdef, udef, mdef, Mdef, rdef);

    fprintf(stderr,
            "genome.dnas         dna index file generated from make_dnas_file and fasta2cisfasta\n"
            "transcripts.gtf     gtf-formatted annotation file defining all transcripts\n"
            "transcripts.expr    expression levels of transcripts, produced from 'sim expression'\n\n"

            "Output Files\n\n"

            "sim.sam             SAM containing alignment definitions\n"
            "sim.frag            fields: contig transcript strand num_fragments num_transcript_mols\n"
            "sim.read1.fastq     fastq containing first-in-pair simulated reads\n"
            "sim.read2.fastq     fastq containing second-in-pair simulated reads\n\n"
            "Note: if 'sim.read2.fastq' is not provided, sim.read1.fastq will contain\n"
            "      intercalated read pairs as first,second,first,second, etc.\n\n"

            "Note 2: if any of the three files are '/dev/null', will not print to that file\n"
            "and execution will be faster.\n\n"

            );

    fprintf(stderr,
            "Possible sampling schemes:\n\n"
            "'cut:cutprob,fmin,fmax,totmol'\twhere:\n\n"
            "cutprob (FLOAT) is the probability [0,1] of a nucleotide bond being cut\n"
            "fmin, fmax (INT) size range of retained fragments\n"
            "totmol (INT): total num transcript mols to generate\n\n"

            "'uniform:stepsize,fmin,fstep,fmax'\twhere:\n\n"

            "stepsize (INT) num bases apart of each start position on sampled transcript\n"
            "fmin, fstep, fmax (INT) e.g. 100,10,150 generates sizes 100,110,120,130,140,150\n\n"
            );
           
    fprintf(stderr,
            "read id format is the following:\n\n"
            "id:read1:D:S:P:C[,D:S:P:C]+:read2:D:S:P:C[,D:S:P:C]+:fragment_size:F\n\n"
            "where\n"
            "D = dna name or contig name\n"
            "S = strand [+ or -]\n"
            "P = start position, (where first base on a contig is 1)\n"
            "C = CIGAR string for this read only\n"
            "F = positive integer fragment size\n\n"
            "Any remaining alignment blocks not on the same 'D' (contig) may be appended with commas.\n"
            "For example:\n\n"
            "1:read1:chr1:+:10592:50M:read2:chr1:-:10792:50M:fragment_size:150\n\n"
            );
    return 1;
}


/********************************************************
OVERVIEW 

In this scheme, there is actual dna sequence provided, and our task is
to simulate reads from subregions (transcripts) within that sequence.
The subregions are defined in terms of a CIGAR coordinate
transformation, which is an expansion from the subregion to the full
dna sequence. (contains M and I states)

Finally, the sampling scheme, which for now produces either single or
paired-end reads, are defined in terms of another CIGAR that goes from
the transcript coordinate system to the reads, and thus is a
contraction (contains M and D states)

These three coordiate systems will be called:

'genome', 'transcript', and 'reads'

Definitions:
Transcript   absolute_position  paired_identity  SAM_sequence  fastq_sequence
pos_genomic  left               first            pos           pos
pos_genomic  right              second           pos           neg
neg_genomic  left               second           pos           pos
neg_genomic  right              first            pos           neg

absolute_position: the position on the reference genome's sense strand.
paired_identity  : whether the read is simulated as the first-in-pair or second-in-pair
SAM_sequence     : which genomic strand is used as source sequence
fastq_sequence   : " " "
(if 'neg', the pos sequence is reverse complemented)

*********************************************************/

void purge_sam_buffer(PAIRED_READ_SET & buffer,
                      PAIRED_READ_SET::iterator start,
                      PAIRED_READ_SET::iterator end,
                      bool flip_query_strand_flag,
                      FILE * first_fastq_fh,
                      FILE * second_fastq_fh,
                      FILE * output_sam_fh)
{
    if (first_fastq_fh != NULL)
    {
        for (PAIRED_READ_SET::iterator it = start; it != end; ++it)
        {
            print_paired_fastq_entries(first_fastq_fh, second_fastq_fh,
                                       (*it).first, (*it).second);
        }
    }
    if (output_sam_fh != NULL)
    {
        for (PAIRED_READ_SET::iterator it = start; it != end; ++it)
        {
            (*it).first->print(output_sam_fh, flip_query_strand_flag);
            (*it).second->print(output_sam_fh, flip_query_strand_flag);
        }
    }
    for (PAIRED_READ_SET::iterator it = start; it != end; ++it)
    {
        delete (*it).first;
        delete (*it).second;
    }
    buffer.erase(start, end);
}


int main_sim_reads(int argc, char ** argv)
{
    char c;
    char dna_directory_def[] = ".";
    size_t read_length_def = 50;
    
    bool do_paired_end_def = false;
    bool pairs_same_strand_def = false;
    size_t random_seed_def = 28374979;
    char somatic_mutation_file_def[] = "";
    char read_sampling_scheme_def[] = "uniform:100,0,100";

    bool do_blind_read_names_def = false;

    char zero_quality_code_def = '@';
    size_t min_median_quality_score_def = 30;

    size_t read_length = read_length_def;
    bool do_paired_end = do_paired_end_def;
    bool pairs_same_strand = pairs_same_strand_def;
    bool do_blind_read_names = do_blind_read_names_def;
    bool ignore_read_ids_for_uniqueness = false;

    bool do_sample_sense_strand = false;

    char read_sampling_scheme[1000];
    strcpy(read_sampling_scheme, read_sampling_scheme_def);

    char dna_directory[1000];
    strcpy(dna_directory, dna_directory_def);

    char qual_file1[1000] = "";
    char qual_file2[1000] = "";

    char zero_quality_code = zero_quality_code_def;
    size_t min_median_quality_score = min_median_quality_score_def;

    char somatic_mutation_file[1000];
    strcpy(somatic_mutation_file, somatic_mutation_file_def);

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB
    size_t max_mem = max_mem_def;

    size_t random_seed = random_seed_def;

    while ((c = getopt(argc, argv, "d:l:pqbhis:f:g:c:u:m:M:r:")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'd': strcpy(dna_directory, optarg); break;
        case 'l': read_length = static_cast<size_t>(atof(optarg)); break;
        case 'p': do_paired_end = true; break;
        case 'q': pairs_same_strand = true; break;
        case 'b': do_blind_read_names = true; break;
        case 'h': do_sample_sense_strand = true; break;
        case 's': strcpy(read_sampling_scheme, optarg); break;
        case 'f': strcpy(qual_file1, optarg); break;
        case 'g': strcpy(qual_file2, optarg); break;
        case 'c': zero_quality_code = optarg[0]; break;
        case 'i': ignore_read_ids_for_uniqueness = true; break;
        case 'u': min_median_quality_score = static_cast<size_t>(atof(optarg)); break;
        case 'm': strcpy(somatic_mutation_file, optarg); break;
        case 'M': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'r': random_seed = static_cast<size_t>(atof(optarg)); break;

        default: return sim_reads_usage(dna_directory_def, 
                                        read_length_def, 
                                        read_sampling_scheme_def,
                                        zero_quality_code_def,
                                        min_median_quality_score_def,
                                        somatic_mutation_file_def,
                                        max_mem_def,
                                        random_seed_def); 
            break;
        }
    }

    int min_req_args = 7;
    int min_arg_count = optind + min_req_args;

    if ((argc != min_arg_count) && argc != (min_arg_count + 1))
    {
        return sim_reads_usage(dna_directory_def,
                               read_length_def, 
                               read_sampling_scheme_def,
                               zero_quality_code_def,
                               min_median_quality_score_def,
                               somatic_mutation_file_def,
                               max_mem_def,
                               random_seed_def);
    }
    
    char * dna_index_file = argv[optind];
    char * gtf_file = argv[optind + 1];
    char * species = argv[optind + 2];
    char * expression_file = argv[optind + 3];
    char * output_sam_file = argv[optind + 4];
    char * output_fragment_file = argv[optind + 5];
    char * output_first_fastq_file = argv[optind + 6];
    char output_second_fastq_file[1000] = "";

    if (argc == min_arg_count + 1)
    {
        strcpy(output_second_fastq_file, argv[optind + min_req_args]);
    }

    FILE * output_fragment_fh = open_if_present(output_fragment_file, "w");
    FILE * output_first_fastq_fh = open_if_present(output_first_fastq_file, "w");
    FILE * output_second_fastq_fh = open_if_present(output_second_fastq_file, "w");

    FILE * output_second_fastq_used_fh =
        output_second_fastq_fh == NULL ? output_first_fastq_fh :
        output_second_fastq_fh;
    
    std::set<SequenceProjection> genome_to_tx = gtf_to_sequence_projection(gtf_file, species);

    std::map<std::string, std::string> transcript_gene_map = gtf_to_transcript_gene_map(gtf_file);

    std::set<SequenceProjection, LessSequenceProjectionTarget> tx_to_genome;
    std::pair<std::set<SequenceProjection, LessSequenceProjectionTarget>::iterator, bool> tx_ins;
    std::set<SequenceProjection>::const_iterator sp_it;
    for (sp_it = genome_to_tx.begin(); sp_it != genome_to_tx.end(); ++sp_it)
    {
        tx_ins = tx_to_genome.insert(InvertProjection((*sp_it)));
        assert(tx_ins.second);
    }

    assert(tx_to_genome.size() == genome_to_tx.size());


    //parse dna file
    cis::dna_collection dnacollection;
    dnacollection.make_ready(dna_index_file, std::string(dna_directory));

    //parse a file of somatic point mutations
    LOCUS_SET somatic_mutations =
    parse_somatic_mutations(somatic_mutation_file);

    const bool do_reverse_second_quals = true;

    bool paired_reads_are_same_stranded = false;

    ReadSampler read_sampler(somatic_mutations, qual_file1, qual_file2,
                             read_length, zero_quality_code,
                             max_mem,
                             do_reverse_second_quals,
                             do_blind_read_names);

    read_sampler.Initialize();

    if (read_sampler.qual1_buf_file.valid != read_sampler.qual2_buf_file.valid)
    {
        fprintf(stderr, "Error: must provide either none or both qual file sources\n");
        exit(1);
    }

    bool do_simulate_errors = read_sampler.qual1_buf_file.valid;

    FragmentGenerator frag_generator;
    frag_generator.Initialize(read_sampling_scheme, expression_file, random_seed);

    if (frag_generator.fragment_length_min < read_length)
    {
        fprintf(stderr, "Error: Specified fragment generation minimum length of %Zu is "
                "shorter than the read length of %Zu.\n"
                "Please specify a higher minimum or shorter read length\n",
                frag_generator.fragment_length_min, read_length);
        exit(1);
    }
    //check that the expression file and set of projections match up

    if (tx_to_genome.size() != frag_generator.transcript_info.size())
    {
        fprintf(stderr, "Error: Number of distinct transcripts in GTF file (%Zu) "
                "does not match number in expression file (%Zu)\n",
                tx_to_genome.size(), frag_generator.transcript_info.size());
        exit(1);
    }
    
    for (std::set<SequenceProjection>::const_iterator sp_it = tx_to_genome.begin();
         sp_it != tx_to_genome.end(); ++sp_it)
    {
        SequenceProjection const& sp = *sp_it;
        unique_transcript tx_from_projection(sp.target_dna.c_str(),
                                             sp.source_dna.c_str(),
                                             sp.same_strand ? '+' : '-');
        
        if (frag_generator.transcript_info.find(tx_from_projection) ==
            frag_generator.transcript_info.end())
        {
            fprintf(stderr, "Error: provided GTF file contains one or more isoforms "
                    "not listed in expression file.  For example, %s on %s (%c).\n",
                    tx_from_projection.transcript_id.c_str(), 
                    tx_from_projection.contig_name.c_str(), 
                    tx_from_projection.strand);
            exit(1);
        }
    }
        
    SamOrder sam_order((ignore_read_ids_for_uniqueness ? SAM_POSITION : SAM_POSITION_RID),
                       "MIN_ALIGN_GUIDE");

    LessSAMLinePair sam_pair_less_fcn(&sam_order);

    FILE * output_sam_fh = open_if_present(output_sam_file, "w");

    //the coordinate transformation will be a transcript->genome input
    //and a transcript->reads.  When fused, this yields
    //genome->reads.  From this, we can retrieve the actual dna sequence

    bool first_read_on_transcript = true;

    //print SAM header
    if (output_sam_fh != NULL)
    {
        fprintf(output_sam_fh, "@HD\tVN:sim\tSO:coordinate\n");
        for (cis::DNAS::const_iterator dna_iter = dnacollection.begin();
             dna_iter != dnacollection.end(); ++dna_iter)
        {
            fprintf(output_sam_fh, "@SQ\tSN:%s\tLN:%Zu\n", (*dna_iter)->name.c_str(),
                    (*dna_iter)->length());
        }
    }

    if (do_paired_end)
    {        
        //CIGAR to describe the relationship between genome dna and matepair
        SamBuffer sam_buffer(&sam_order, paired_reads_are_same_stranded);

        for (std::set<SequenceProjection>::const_iterator 
             sp_it = tx_to_genome.begin(); sp_it != tx_to_genome.end(); ++sp_it)
        {
            SequenceProjection const& sp = *sp_it;

            Cigar::CIGAR_ITER low_trimmed, high_trimmed;
            Cigar::Trim(sp.cigar, true, &low_trimmed, &high_trimmed);
            size_t span_on_genome = Cigar::Length(low_trimmed, high_trimmed, false);

            int64_t start_on_genome = Cigar::LeftOffset(sp.cigar, false);
            int64_t end_on_genome = start_on_genome + span_on_genome;

            cis::dna_t const* dna = 
            GetDNAByName(dnacollection, sp.species, sp.target_dna);

            if (dna == NULL)
            {
                fprintf(stderr, "Skipping contig %s missing from sequence input\n",
                        sp.target_dna.c_str());
                continue;
            }


            LOCUS_SET::const_iterator 
            somatic_mutation_start = somatic_mutations.end(), 
            somatic_mutation_end = somatic_mutations.end();

            //find the (small handful) of possible somatic mutations
            //to apply to this source dna
            if (! somatic_mutations.empty())
            {
                Locus low_test(sp.species.c_str(), dna->name.c_str(), start_on_genome, 'A', 'A');
                somatic_mutation_start = somatic_mutations.lower_bound(low_test);

                Locus high_test(sp.species.c_str(), dna->name.c_str(), end_on_genome, 'A', 'A');
                somatic_mutation_end = somatic_mutations.upper_bound(high_test);
            }

            std::pair<SamLine const*, SamLine const*> read_pair;

            unique_transcript tx_from_projection(sp.target_dna.c_str(),
                                                 sp.source_dna.c_str(),
                                                 sp.same_strand ? '+' : '-');

            FragmentGenerator::BOUNDS bounds = 
                frag_generator.Sample(tx_from_projection);

            TranscriptInfo const& tinfo = (*frag_generator.transcript_info.find(tx_from_projection)).second;

            FragmentGenerator::BOUNDS::const_iterator bounds_iter;
            
            // fprintf(stderr, "transcript %s, contig %s at %Zu. Num fragments: %Zu\n",
            //         tx_from_projection.transcript_id.c_str(),
            //         tx_from_projection.contig_name.c_str(),
            //         Cigar::LeftOffset(sp.cigar, false),
            //         bounds.size()
            //         );

            fflush(stdout);

            first_read_on_transcript = true;
            size_t this_num_fragments = bounds.size();

            for (bounds_iter = bounds.begin(); bounds_iter != bounds.end(); ++bounds_iter)
            {
                size_t start_pos = (*bounds_iter).first;
                size_t end_pos = (*bounds_iter).second;

                read_pair = read_sampler.sample_pair(sp, dna, start_pos, end_pos,
                                                     min_median_quality_score,
                                                     do_simulate_errors,
                                                     do_sample_sense_strand);
                
                std::pair<SamLine const*, bool> first_inserted = sam_buffer.insert(read_pair.first);
                assert(first_inserted.second);

                std::pair<SamLine const*, bool> second_inserted = sam_buffer.insert(read_pair.second);
                if (! second_inserted.second)
                {
                    --this_num_fragments;
                }
                
                else
                {
                    if (first_read_on_transcript)
                    {
                        sam_buffer.safe_advance_lowbound(read_pair.first);
                        first_read_on_transcript = false;
                    }
                    
                }
                
            }

            std::string gene_id = transcript_gene_map[tx_from_projection.transcript_id];
            fprintf(output_fragment_fh, "%s\t%s\t%s\t%c\t%Zu\t%Zu\t%Zu\n",
                    tx_from_projection.contig_name.c_str(),
                    gene_id.c_str(),
                    tx_from_projection.transcript_id.c_str(),
                    tx_from_projection.strand,
                    tinfo.length,
                    tinfo.count,
                    this_num_fragments);

            sam_buffer.purge(output_sam_fh, 
                             output_first_fastq_fh,
                             output_second_fastq_used_fh,
                             false);

            // purge_sam_buffer(sam_buffer, sam_buffer.begin(), sam_buffer.end(),
            //                  flip_query_strand_flag,
            //                  output_first_fastq_fh,
            //                  output_second_fastq_used_fh,
            //                  output_sam_fh);

        }

        sam_buffer.purge(output_sam_fh, 
                         output_first_fastq_fh,
                         output_second_fastq_used_fh,
                         true); //purge remainder

        // purge_sam_buffer(sam_buffer, sam_buffer.begin(), sam_buffer.end(),
        //                  flip_query_strand_flag,
        //                  output_first_fastq_fh,
        //                  output_second_fastq_used_fh,
        //                  output_sam_fh);

    } // if (do_paired_end)

    //clean up
    close_if_present(output_sam_fh);
    close_if_present(output_first_fastq_fh);
    close_if_present(output_second_fastq_fh);
    close_if_present(output_fragment_fh);

    return 0;
}
