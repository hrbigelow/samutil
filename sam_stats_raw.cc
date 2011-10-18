#include "sam_stats.h"
#include "sam_stats_aux.h"
#include "matrix_tools.h"
#include "sam_helper.h"
#include "cigar_ops.h"
#include "dep/tools.h"

#include <cassert>

#include <map>
#include <string>
#include <set>


//generate the temporary raw stats output file
int main_raw(int argc, char ** argv)
{
    
    const size_t num_consensus_match = 2;
    const size_t num_strands = 2;

    size_t mcs_default = 40;
    size_t mqs_default = 40;
    size_t nas_default = 200;
    size_t mrl_default = 76;
    size_t cbs_default = 10000;
    bool smp_default = false;

    size_t max_consensus_quality = mcs_default; // c
    size_t max_quality_score = mqs_default; // q
    size_t num_alignment_scores = nas_default; // a
    size_t max_read_length = mrl_default; // r
    size_t consensus_buf_size = cbs_default;
    bool skip_max_prescan = smp_default;

    //const bool use_top_coord = true;

    //sam format files give ones-based positions
    const bool sam_is_ones_based = true;

    extern char *optarg;
    char c;

    
    while ((c = getopt(argc, argv, "sc:q:a:r:b:")) != -1) {
        switch(c) {
        case 's': skip_max_prescan = true; break;
        case 'c': max_consensus_quality = static_cast<size_t>(atof(optarg)); break;
        case 'q': max_quality_score = static_cast<size_t>(atof(optarg)); break;
        case 'a': num_alignment_scores = static_cast<size_t>(atof(optarg)); break;
        case 'r': max_read_length = static_cast<size_t>(atof(optarg)); break;
        case 'b': consensus_buf_size = static_cast<size_t>(atof(optarg)); break;
        default: return main_usage(); break;
        }
    }

    if (argc - optind != 2)
    {
        return raw_usage(smp_default, mcs_default, mqs_default, 
                         nas_default, mrl_default, cbs_default);
    }

    char * sam_file = argv[optind];
    char * pileup_file = argv[optind + 1];

    if (! skip_max_prescan)
    {
        max_consensus_quality = find_max_consensus_quality(pileup_file);
        samfile_limits(sam_file, 
                       & max_read_length,
                       & max_quality_score,
                       & num_alignment_scores);

        fprintf(stderr,
                "Max consensus_quality: %Zu\n"
                "Max read length      : %Zu\n"
                "Max base quality     : %Zu\n"
                "Num alignment scores : %Zu\n",
                max_consensus_quality,
                max_read_length,
                max_quality_score,
                num_alignment_scores);
    }

    size_t num_read_positions = max_read_length;
    size_t num_base_qualities = max_quality_score + 1;

    FILE * sam_fh = fopen(sam_file, "r");
    if (sam_fh == NULL)
    {
        fprintf(stderr, "Couldn't open SAM input file %s\n", sam_file);
        exit(1);
    }

    FILE * pileup_fh = fopen(pileup_file, "r");
    if (pileup_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup file %s\n", pileup_file);
        exit(1);
    }

    size_t num_consensus_qualities = max_consensus_quality + 1; // (must include 0 as a cons score)

    //  BaseQualStrandReader reader;
    //  reader.initialize(pileup_file);

    //assume source coordinates are given with ones_displacement.
    size_t const ones_offset = 1;

    //char cur_chrom[100];
    //size_t cur_read_start;

    size_t const num_features = 6;

    size_t dim_sizes[] =
        {
            max_read_length,
            num_consensus_match,
            num_consensus_qualities,
            num_base_qualities,
            num_strands,
            num_alignment_scores
        };

    size_t digits[num_features + 1];
    size_t * digit_chunks[num_features];

    digits[0] = 1;
    for (size_t d = 0; d != 6; ++d)
    {
        digits[d + 1] = digits[d] * dim_sizes[d];
        digit_chunks[d] = new size_t[dim_sizes[d]];
        for (size_t d2 = 0; d2 != dim_sizes[d]; ++d2)
        {
            digit_chunks[d][d2] = d2 * digits[d];
        }

    }

    //initialize / allocate counts to zero
    size_t num_categories = digits[6];

    int * alignment_stats = new int[num_categories];
    std::fill(alignment_stats, alignment_stats + num_categories, 0);

    int * consensus_score_buf = new int[consensus_buf_size];
    char * consensus_base_buf = new char[consensus_buf_size];

    char prev_refname[1000];
    strcpy(prev_refname, "");

    char prev_end_contig_name[1000];
    //char cur_end_contig_name[1000];

    strcpy(prev_end_contig_name, "");

    //size_t prev_end_genome_coord = 0;
    //size_t cur_end_genome_coord;

    bool has_alignment_score;

    SamLine *samline;

    //describes the set of contiguous positions of consensus held in the buffer
    Cigar::CIGAR_VEC metacontig_to_buf_cigar;

    std::map<int, int> alignment_scores;
    size_t alignment_score_max_index = 0;
    size_t alignment_score_index;

    //parse the header
    char sequence_name[1000];
    size_t sequence_length;

    int alignment_score;
    int strand;

    META_CONTIG contig_offsets;
    META_CONTIG contig_lengths;
    META_CONTIG::iterator contig_iter;

    size_t total_sequence_length = 0;

    while ((samline = new SamLine(sam_fh, sam_is_ones_based)) != NULL)
    {
        if (samline->parse_flag == HEADER)
        {
            if (2 == sscanf(samline->tag_string, "SQ\tSN:%[^\t]\tLN:%zu", sequence_name, &sequence_length))
            {
                contig_lengths[std::string(sequence_name)] = sequence_length;
                contig_offsets[std::string(sequence_name)] = total_sequence_length;
                total_sequence_length += sequence_length;
            }
            delete samline;
        }
        else
        {
            break;
        }
        
    }

    //int consensus_score;
    //char consensus_base;
    int consensus_match;
    size_t quality_score;

    //size_t begin_to_buf_offset = 0; // distance from begin to buffer
    size_t begin_to_contig_offset; // distance from 'begin' to start of contig
    //size_t buf_to_read_offset; // from start of buffer to start of read
    size_t contig_to_read_offset; // from begin of contig to begin of read (found in SAM file alignment)

    metacontig_to_buf_cigar =
        fill_buffer(pileup_fh, contig_offsets, Cigar::CIGAR_VEC(), ones_offset,
                    consensus_buf_size, max_consensus_quality,
                    consensus_score_buf, consensus_base_buf);

    //size_t metacontig_content_length = Cigar::Length(metacontig_to_buf_cigar, use_top_coord);
    //size_t buf_content_length = Cigar::Length(metacontig_to_buf_cigar, ! use_top_coord);

    size_t total_match_bases = 0;
    size_t num_reads_read = 1; // we start this loop after first read is read.

    while (samline->parse_flag != END_OF_FILE)
    {
        if (samline->query_unmapped())
        {
            //!!! do something here
            // fprintf(stderr, "Warning: %s unmapped\n",
            //         samline->qname);
            continue;
        }

        if (samline->raw_read_length() > max_read_length)
        {
            fprintf(stderr, "Error: parsed read length %Zu exceeds maximum length"
                    " specified of %Zu.\n"
                    "Please re-run with larger maximum specified length.\n",
                    samline->raw_read_length(),
                    max_read_length);
            exit(1);
        }

        if (num_reads_read % 10000 == 0)
        {
            fprintf(stderr, "Processed read %Zu\n", num_reads_read);
            fflush(stderr);
        }

        bool new_contig = strcmp(samline->rname, prev_refname) != 0;

        if (new_contig)
        {
            begin_to_contig_offset = contig_offsets[std::string(samline->rname)];
        }

        contig_to_read_offset = samline->zero_based_pos();
        size_t begin_to_read_offset = begin_to_contig_offset + contig_to_read_offset;


        Cigar::CIGAR_VEC metacontig_to_read_cigar = 
            Cigar::FromString(samline->cigar, begin_to_read_offset);

        // !!! how inefficient is this?
        Cigar::CIGAR_VEC buf_to_read_cigar =
            Cigar::TransitiveMerge(metacontig_to_buf_cigar, 
                                   metacontig_to_read_cigar, false);

        // size_t metalength1 = Cigar::Length(metacontig_to_buf_cigar, true);
        // size_t metalength2 = Cigar::Length(metacontig_to_read_cigar, true);

        // size_t readlen1 = Cigar::Length(metacontig_to_read_cigar, false);
        // size_t readlen2 = Cigar::Length(buf_to_read_cigar, false);

        // assert(readlen1 == readlen2);

        // size_t buflen1 = Cigar::Length(metacontig_to_buf_cigar, false);

        if ((*buf_to_read_cigar.rbegin()).op == Cigar::I)
        {
            //a section at the end of the read exceeds the buffer
            if (feof(pileup_fh))
            {
                fprintf(stderr, "error: read alignment exceeds pileup coordinates "
                        "for read %s at %s:%Zu\n", 
                        samline->qname, samline->rname, samline->ones_based_pos());
                exit(1);
            }

            //hopefully, this means that the beginning of the read
            //lies within the buffer, and so we can advance the buffer
            //to resolve it.
            assert(buf_to_read_cigar[0].op == Cigar::D);

            size_t leap_size = buf_to_read_cigar[0].length;

            metacontig_to_buf_cigar =
                advance_consensus_buffers(pileup_fh, 
                                          contig_offsets,
                                          metacontig_to_buf_cigar,
                                          leap_size,
                                          ones_offset,
                                          max_consensus_quality,
                                          consensus_score_buf,
                                          consensus_base_buf);

            buf_to_read_cigar =
                Cigar::TransitiveMerge(metacontig_to_buf_cigar, 
                                       metacontig_to_read_cigar, false);
        }

        //set the pointers of consensus_score and consensus_base to
        //the point in the buffer corresponding to the beginning of
        //the aligned read.

        strand = samline->query_on_pos_strand() ? 0 : 1;
        
        alignment_score = samline->alignment_score(& has_alignment_score);

        // assert(has_alignment_score);
        
        if (alignment_scores.find(alignment_score) == alignment_scores.end())
        {
            alignment_scores[alignment_score] = alignment_score_max_index++;
        }

        if (alignment_scores.size() > num_alignment_scores)
        {
            fprintf(stderr, "Error: number of distinct alignment scores exceeds"
                    " estimate for the maximum number of alignment scores %Zu."
                    " Please re-run and choose a higher estimate.",
                    num_alignment_scores);
            exit(1);
        }

        alignment_score_index = alignment_scores[alignment_score];

        //how to efficiently do a dual traversal of the buf_to_read_cigar?
        size_t read_position = 0;
        size_t buf_position = 0;
        size_t stats_index;

        for (size_t chunk = 0; chunk != buf_to_read_cigar.size(); ++chunk)
        {
            Cigar::Unit & block = buf_to_read_cigar[chunk];
            switch (block.op)
            {
            case Cigar::M:
                for (size_t p = 0; p != block.length; ++p, ++read_position, ++buf_position)
                {
                    //consensus_score = consensus_score_buf[buf_position];
                    // consensus_base = consensus_base_buf[buf_position];
                    char read_base = toupper(samline->seq[read_position]);
                    
                    consensus_match = 
                        consensus_base_buf[buf_position] == read_base ? 1 : 0;

                    quality_score = QualityCodeToQuality(samline->qual[read_position]);

                    if (quality_score > max_quality_score)
                    {
                        fprintf(stderr, "Error: found quality score %Zu exceeds maximum"
                                " chosen quality score %Zu.\n"
                                "Please rerun with a higher maximum quality score\n",
                                quality_score, max_quality_score);
                        exit(1);
                    }

                    ++total_match_bases;
                    stats_index = 
                        digit_chunks[0][read_position]
                        + digit_chunks[1][consensus_match]
                        + digit_chunks[2][consensus_score_buf[buf_position]]
                        + digit_chunks[3][quality_score]
                        + digit_chunks[4][strand]
                        + digit_chunks[5][alignment_score_index];
                    
                    alignment_stats[stats_index]++;
                    
                }
                break;

            case Cigar::I:
            case Cigar::S:
                //insertion, going from buffer-to-read.  (read has extra bases)
                read_position += block.length;
                break;
                
            case Cigar::D:
                //buffer has extra bases
                buf_position += block.length;
                break;

            case Cigar::N:
            case Cigar::H:
            case Cigar::P:
            case Cigar::None:
                fprintf(stderr, "Cigar ops N, H, None, and P not implemented yet\n");
                exit(1);
                break;
            }
        }

        if (new_contig)
        {
            strcpy(prev_refname, samline->rname);
        }

        ++num_reads_read;
        delete samline;
        samline = new SamLine(sam_fh, sam_is_ones_based);
                  
    }

    //printf("total match bases: %Zu\n", total_match_bases);

    //output raw memory, 
    printf("%Zu\t%Zu\t%Zu\t%Zu\t%Zu\t%Zu\n",
           num_read_positions,
           num_consensus_match,
           num_consensus_qualities,
           num_base_qualities,
           num_strands,
           num_alignment_scores);

    for (std::map<int, int>::iterator a = alignment_scores.begin();
         a != alignment_scores.end(); ++a)
    {
        printf("%i\t%i\n", (*a).second, (*a).first);
    }

    fflush(stdout);

    //     char sep;
    //     for (size_t s = 0; s != num_categories; ++s)
    //     {
    //         sep = (s % 10 == 0) ? '\n' : '\t';
    //         printf("%i%c", alignment_stats[s], sep);
    //     }
    //     printf("\n");

    fwrite(alignment_stats, sizeof(alignment_stats[0]),
           num_categories, stdout);

    delete alignment_stats;
    delete consensus_score_buf;
    delete consensus_base_buf;

    for (size_t d = 0; d != num_features; ++d)
    {
        delete digit_chunks[d];
    }

    return 0;
}
