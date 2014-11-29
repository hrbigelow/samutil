#include "sam_stats_aux.h"
#include "dep/pileup_tools.h"
#include "sam_line.h"
#include "dep/tools.h"

#include <cstdio>
#include <map>
#include <string>
#include <cassert>
#include <set>


size_t find_max_consensus_quality(char const* pileup_file)
{

    //prescan pileup file for max consensus score
    FILE * pileup_fh = fopen(pileup_file, "r");
    if (pileup_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup file %s\n", pileup_file);
        exit(1);
    }

    int indel_histo_size = 0;
    bool with_soapsnp = true;
    bool do_load_alignment = true;
    size_t ones_offset = 1;

    PileupSummary pileup(indel_histo_size, with_soapsnp, 
                         ! do_load_alignment, ones_offset);

    int max_consensus_quality = 0;

    while (pileup.load_line(pileup_fh))
    {
        if (pileup._is_indel_line)
        {
            continue;
        }
        max_consensus_quality = std::max(max_consensus_quality,
                                         pileup._consensus_quality);
    }
    fclose(pileup_fh);
    return static_cast<size_t>(max_consensus_quality);
}


//find maximum read length, quality, total number of alignment scores
//in a samfile
void samfile_limits(char const* sam_file, 
                    size_t * max_read_length,
                    size_t * max_base_quality,
                    size_t * num_distinct_align_scores)
{

    FILE * sam_fh = fopen(sam_file, "r");
    if (sam_fh == NULL)
    {
        fprintf(stderr, "Couldn't open SAM input file %s\n", sam_file);
        exit(1);
    }

    SamLine * samline;

    *max_read_length = 0;
    *max_base_quality = 0;
    *num_distinct_align_scores = 0;

    std::set<int> alignment_scores;
    bool has_alignment_score; // unused
    size_t quality_score;

    bool source_is_ones_based = true;

    while ((samline = new SamLine(sam_fh, source_is_ones_based)) != NULL)
    {
        if (samline->parse_flag == HEADER)
        {
            delete samline;
            continue;
        }
        *max_read_length = std::max(*max_read_length,
                                    samline->raw_read_length());
        
        alignment_scores.insert(samline->alignment_score(& has_alignment_score));
        char * qp; // qual pointer
        for (qp = samline->qual; *qp != '\0'; ++qp)
        {
            quality_score = QualityCodeToQuality(*qp);
            *max_base_quality = std::max(*max_base_quality, quality_score);
        }
        delete samline;
    }
    *num_distinct_align_scores = alignment_scores.size();
}


// fill 'fill_size' positions of a pileup into 'consensus_score_buf'
// and 'consensus_base_buf'
// the 
Cigar::CIGAR_VEC fill_buffer(FILE * pileup_fh,
                             std::map<std::string, size_t> const& contig_offsets,
                             Cigar::CIGAR_VEC const& prefix_metacontig_to_buf_cigar,
                             size_t ones_offset,
                             size_t fill_size,
                             size_t max_consensus_quality,
                             int * consensus_score_buf,
                             char * consensus_base_buf)
{

    size_t buf_pos;
    size_t cur_meta_contig_left_bound; //coordinate along concatenated set of alphabetized contigs
    size_t match_length = 0;

    size_t begin_to_contig_offset = 0;

    //initialize with a copy of the prefix, and start appending
    Cigar::CIGAR_VEC full_metacontig_to_buf_cigar(prefix_metacontig_to_buf_cigar);

    char prev_contig_name[1000];
    prev_contig_name[0] = '\0';

    size_t prec_meta_contig_right_bound = Cigar::Length(prefix_metacontig_to_buf_cigar, false);

    int indel_histo_size = 0;
    bool with_soapsnp = true;
    bool do_load_alignment = true;

    PileupSummary pileup(indel_histo_size, with_soapsnp, 
                         ! do_load_alignment, ones_offset);

    for (buf_pos = 0; buf_pos != fill_size; ++buf_pos)
    {
        if (! pileup.load_line(pileup_fh))
        {
            //we've reached the end of the file
            break;
        }

        if (pileup._is_indel_line)
        {
            //we are just ignoring indel lines for now
            --buf_pos;
            continue;
        }

        bool new_contig = strcmp(pileup._reference, prev_contig_name) != 0;

        if (new_contig)
        {
            std::map<std::string, size_t>::const_iterator cit =
                contig_offsets.find(std::string(pileup._reference));

            if (cit == contig_offsets.end())
            {
                fprintf(stderr, "Error: advance_consensus_buffers: "
                        "didn't find %s in list of contigs\n",
                        pileup._reference);
                exit(1);
            }
            else
            {
                begin_to_contig_offset = (*cit).second;
            }
        }

        cur_meta_contig_left_bound = begin_to_contig_offset + pileup._position;

        if (pileup._consensus_quality > static_cast<int>(max_consensus_quality))
        {
            fprintf(stderr, 
                    "Error: fill_buffer: pileup contained consensus quality %i\n"
                    "exceeding maximum consensus quality %Zu.\n"
                    "Please re-run with higher setting for max consensus quality\n",
                    pileup._consensus_quality,
                    max_consensus_quality);
            exit(1);
        }

        consensus_score_buf[buf_pos] = pileup._consensus_quality;
        consensus_base_buf[buf_pos] = toupper(pileup._consensus_base);

        assert(cur_meta_contig_left_bound >= prec_meta_contig_right_bound);

        if (cur_meta_contig_left_bound > prec_meta_contig_right_bound)
        {
            //this coord is non-contiguous with the preceding one.
            //the full_metacontig_to_buf_cigar should reflect this.
            if (buf_pos != 0)
            {
                //finished a legitimate match block
                full_metacontig_to_buf_cigar.push_back(Cigar::Unit(Cigar::M, match_length));
            }

            //also finished a delete block
            full_metacontig_to_buf_cigar.push_back
                (Cigar::Unit(Cigar::D, cur_meta_contig_left_bound - prec_meta_contig_right_bound));

            match_length = 0;
        }

        
        match_length++;

        prec_meta_contig_right_bound = cur_meta_contig_left_bound + 1;
        strcpy(prev_contig_name, pileup._reference);
    }

    //at the end, close out the current match state
    full_metacontig_to_buf_cigar.push_back(Cigar::Unit(Cigar::M, match_length));

    return full_metacontig_to_buf_cigar;
}


//advance the consensus buffer by 'leap_size'.  can this be
//generalized to do the initial buffer loading?
//calling with a leap size of zero initially fills the buffer
//assume 

Cigar::CIGAR_VEC 
advance_consensus_buffers(FILE * pileup_fh,
                          std::map<std::string, size_t> const& contig_offsets,
                          Cigar::CIGAR_VEC const& metacontig_to_buf_cigar,
                          size_t leap_size,
                          size_t ones_offset,
                          size_t max_consensus_quality,
                          int * consensus_score_buf,
                          char * consensus_base_buf)
{
    if (feof(pileup_fh))
    {
        fprintf(stderr, "advance_consensus_buffers: Cannot advance, end of file reached\n");
        exit(1);
    }

    bool const use_top_coord = true;

    //CIGAR(pileup_buffer, meta_contig_coord)
    size_t prev_buf_content_size = 
        Cigar::Length(metacontig_to_buf_cigar, ! use_top_coord);

    //sizes are given in buffer coordinates (top coord
    size_t overlap = prev_buf_content_size - leap_size;


    int * score_next_ptr = consensus_score_buf + overlap;
    char * base_next_ptr = consensus_base_buf + overlap;

    //copy the suffix of the buf to the beginning.
    std::copy(consensus_score_buf + leap_size, 
              consensus_score_buf + prev_buf_content_size,
              consensus_score_buf);

    std::copy(consensus_base_buf + leap_size,
              consensus_base_buf + prev_buf_content_size,
              consensus_base_buf);

    //initialize the new_cigar with the last <overlap> bases 
    Cigar::CIGAR_VEC prefix_cigar = 
        Cigar::Substring(metacontig_to_buf_cigar, prev_buf_content_size - overlap, 
                         prev_buf_content_size, ! use_top_coord);
    
    Cigar::CIGAR_VEC full_metacontig_to_buf_cigar =
        fill_buffer(pileup_fh, contig_offsets, 
                    prefix_cigar, ones_offset, leap_size, 
                    max_consensus_quality,
                    score_next_ptr, base_next_ptr);

    return full_metacontig_to_buf_cigar;
}
    
