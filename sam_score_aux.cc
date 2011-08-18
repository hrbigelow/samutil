#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "sam_score_aux.h"
#include "cigar_ops.h"
#include "dep/tools.h"

#include <algorithm>
#include <cstdio>
#include <cmath>

CollapseScoreTriplet::CollapseScoreTriplet(size_t _min, size_t _max) : 
    min_score(_min), max_score(_max) 
{ 
    if (max_score <= min_score)
    {
        fprintf(stderr, "Error: CollapseScoreTriplet: "
                "max_score must be strictly larger than min_score\n");
        exit(1);
    }
    this->shift_bits = static_cast<size_t>(ceilf(log2f(static_cast<float>(max_score - min_score))));
}

size_t CollapseScoreTriplet::operator()(size_t top, size_t sec, size_t given) const
{
    size_t result = top - this->min_score;
    result = result<<this->shift_bits;
    result += sec - this->min_score;
    result = result<<this->shift_bits;
    result += given - this->min_score;
    return result;
}


/*
  A packed score is equal to (fragment_score * 2) + 1 for genome alignments,
  or (fragment_score * 2) for transcriptome alignments,
  or default_missing for anything not aligned.

  The packed score is used to stratify by the combination of
  fragment_score X G/T space.

 */

//Calculate raw fragment score as a simple sum of raw read scores, but
//downgraded by one point if in the genome, and censored by a fragment
//length range.  Warning: depends on the isize field, which is
//incorrect for genome-projected transcript alignments.
//Returns the packed raw score, where the 'default missing' value is equal
//to the max valid score + 1
RAW_SCORE_T PairedEndPackedScore(SamLine const* a,
                                 SamLine const* b,
                                 char const* raw_score_tag,
                                 RAW_SCORE_T max_valid_fragment_score, /* numerically maximum */
                                 size_t min_allowed_fragment_length,
                                 size_t max_allowed_fragment_length)
{
    bool a_has_score;
    bool b_has_score;

    RAW_SCORE_T fragment_score;

    RAW_SCORE_T a_raw_score = 
        a->alignment_score(raw_score_tag, max_valid_fragment_score + 1, &a_has_score);

    RAW_SCORE_T b_raw_score = 
        b->alignment_score(raw_score_tag, max_valid_fragment_score + 1, &b_has_score);

    size_t fragment_length = static_cast<size_t>(a->isize);

    assert(a->isize >= 0);
    char xp_tag[2];
    char tag_type;

    if (fragment_length >= min_allowed_fragment_length
        && fragment_length <= max_allowed_fragment_length
        && (a_raw_score + b_raw_score) <= max_valid_fragment_score)
    {
        
        int add = (a->has_tag("XP", tag_type, xp_tag) && strcmp(xp_tag, "G") == 0)
            ? 1 : 0;
        fragment_score = (a_raw_score + b_raw_score) * 2 + add;
    }
    else
    {
        fragment_score = max_valid_fragment_score + 1;
    }

    return fragment_score;
}


std::vector<RAW_SCORE_T> 
GetPackedScores(SamBuffer const& sam_buffer, 
                size_t min_allowed_fragment_length,
                size_t max_allowed_fragment_length,
                char const* raw_score_tag,
                RAW_SCORE_T max_valid_fragment_score)
{
    std::vector<RAW_SCORE_T> packed_scores(sam_buffer.unique_entry_pairs.size());
                                          
    PAIRED_READ_SET::iterator pit;

    size_t element_position = 0;
    for (pit = sam_buffer.unique_entry_pairs.begin();
         pit != sam_buffer.unique_entry_pairs.end(); ++pit)
    {
        SamLine const* first = ((*pit).first);
        SamLine const* second = ((*pit).second);

        RAW_SCORE_T packed_score = 
            PairedEndPackedScore(first, second, raw_score_tag, 
                                 max_valid_fragment_score,
                                 min_allowed_fragment_length,
                                 max_allowed_fragment_length);
        
        packed_scores[element_position++] = packed_score;
    }
    return packed_scores;
}


RAW_SCORE_T UnpackScore(RAW_SCORE_T packed_score,
                        RAW_SCORE_T max_valid_fragment_score)
{
    return packed_score == max_valid_fragment_score + 1 ? packed_score : packed_score>>1;
}

//assume sam_buffer has alignments of only one physical
//fragment in its 'unique_entry_pairs' buffer.
//calculate all fragment_scores, and return the top two


/*
1. load all fragment scores into the map of raw_score -> rank

2. initialize map rank values in order.

3. look up and cache mapq calibration sub-table based on top two
ranking raw scores

4. set all rank, primary flag, and mapq calibration from lookups
*/
void SetScoreFields(SamBuffer const& sam_buffer, 
                    size_t min_fragment_length,
                    size_t max_fragment_length,
                    char const* raw_score_tag,
                    RAW_SCORE_T max_valid_fragment_score,
                    bool larger_raw_score_is_better,
                    int const* score_cpd,
                    CollapseScoreTriplet const& score_triplet)
{
                           
    PAIRED_READ_SET::iterator pit;

    RawScoreBetter score_order(larger_raw_score_is_better);
    //rank, counts
    std::map<RAW_SCORE_T, std::pair<size_t, size_t>, RawScoreBetter> strata(score_order);
    std::map<RAW_SCORE_T, std::pair<size_t, size_t>, RawScoreBetter>::iterator mit;

    //1. first traversal.  for each pair: load (if not exists)
    //fragment score into the map of raw_score -> rank
    std::vector<RAW_SCORE_T> packed_scores = 
        GetPackedScores(sam_buffer,
                        min_fragment_length,
                        max_fragment_length,
                        raw_score_tag,
                        max_valid_fragment_score);
    
    size_t elem = 0;
    for (elem = 0; elem != packed_scores.size(); ++elem)
    {
        mit = strata.find(packed_scores[elem]);
        if (mit == strata.end())
        {
            mit = strata.insert(strata.end(), std::make_pair(packed_scores[elem], std::make_pair(0, 0)));
        }
        //increment counts;
        ++((*mit).second).second;
    }

    //2. initialize map rank values in order
    size_t rank = 1;
    for (mit = strata.begin(); mit != strata.end(); ++mit)
    {
        ((*mit).second).first = rank++;
    }

    //3. look up and cache mapq calibration sub-table based on top two
    //ranking raw scores
    RAW_SCORE_T top_fragment_score = 
        strata.size() > 0 ? UnpackScore((*strata.begin()).first, max_valid_fragment_score) 
        : max_valid_fragment_score + 1;

    RAW_SCORE_T sec_fragment_score = 
        strata.size() > 1 ? UnpackScore((*++strata.begin()).first, max_valid_fragment_score) 
        : max_valid_fragment_score + 1;

    //4. Set all rank, primary flag, and mapq scores
    elem = 0;
    char rank_tag_value[10];
    char size_tag_value[10];

    size_t new_mapq;
    bool first_encountered_top_stratum = true;
    size_t stratum_size;
    
    for (pit = sam_buffer.unique_entry_pairs.begin();
         pit != sam_buffer.unique_entry_pairs.end(); ++pit)
    {
        SamLine * first = const_cast<SamLine *>((*pit).first);
        SamLine * second = const_cast<SamLine *>((*pit).second);

        rank = strata[packed_scores[elem]].first;
        sprintf(rank_tag_value, "%zu", rank);

        stratum_size = strata[packed_scores[elem]].second;
        sprintf(size_tag_value, "%zu", stratum_size);

        //set rank tag
        first->add_tag("XY", 'i', rank_tag_value);
        first->add_tag("XZ", 'i', size_tag_value);

        second->add_tag("XY", 'i', rank_tag_value);
        second->add_tag("XZ", 'i', size_tag_value);

        //set primary flag
        if (rank == 1 && first->mapped_in_proper_pair() && first_encountered_top_stratum)
        {
            //alignment is primary.
            first->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
            second->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
            first_encountered_top_stratum = false;
        }
        else if (first->query_unmapped())
        {
            //alignment is 'primary' per downstream weird interpretations
            first->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
            second->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
        }
        else
        {
            first->flag |= SamFlags::ALIGNMENT_NOT_PRIMARY;
            second->flag |= SamFlags::ALIGNMENT_NOT_PRIMARY;
        }

        RAW_SCORE_T fragment_score = UnpackScore(packed_scores[elem], max_valid_fragment_score);

        new_mapq = 
            score_cpd[score_triplet(top_fragment_score, sec_fragment_score, fragment_score)];

        bool mapq_correct = first->query_unmapped() ? (new_mapq == 0) : true;

        assert(mapq_correct);

        first->mapq = new_mapq;
        second->mapq = new_mapq;

        ++elem;
    }
}




//count the number of correctly aligned bases between this samline and
//its guide position.  also count and number of bases in the guide and test.
//if number correct is zero, set number of guide and test to zero as well
size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index,
                         size_t * num_bases_guide,
                         size_t * num_bases_test)
{

    size_t num_correct_bases = 0;
    *num_bases_guide = 0;
    *num_bases_test = 0;

    if (strcmp(samline->rname, guide_coords.contig) == 0
        && samline->query_on_pos_strand() == guide_coords.pos_stranded)
    {
        //proceed to measure actual base overlap
        Cigar::CIGAR_VEC test_cigar = Cigar::FromString(samline->cigar, samline->zero_based_pos());
        Cigar::CIGAR_VEC merge_cigar = 
            Cigar::TransitiveMerge(guide_cigar, guide_cigar_index, 
                                   test_cigar, true, false);
        
        num_correct_bases = 
            Cigar::CountAlignedPositions(merge_cigar, num_bases_guide, num_bases_test);
    }
    return num_correct_bases;
}


//parse the next SAM entry and load it into sam_buffer.  assume
//entries are sorted by SAM_RID_POSITION.  upon the occurrence of an
//entry whose fragment is new (or end of file), check that all
//pervious entries are properly paired.
void NextLine(FILE * unscored_sam_fh, 
              SamBuffer & sam_buffer,
              bool allow_absent_seq_qual,
              bool * new_fragment, 
              bool * seen_a_read,
              char * prev_qname,
              size_t * prev_fragment_id,
              SamLine ** low_bound)
{
    SamLine * samline = new SamLine(unscored_sam_fh, allow_absent_seq_qual);

    switch (samline->parse_flag)
    {
    case END_OF_FILE: 
        delete samline;
        if ((*low_bound) != NULL)
        {
            // last time around; purge everything.
            delete (*low_bound);
        }

        *new_fragment = true;
        if (! sam_buffer.yet_unpaired_entries.empty())
        {
            //all entries should be properly paired when we encounter a new fragment.  violation.
            fprintf(stderr, "Error: Reached end of file, yet this entry is paired in sequencing"
                    " but as yet have not found its mate:\n");
            (*sam_buffer.yet_unpaired_entries.begin())->print(stderr, true);
            exit(1);
        }
        break;

    case PARSE_ERROR:
        fprintf(stderr, "Parse error in input sam file");
        exit(1);
        break;

    case HEADER:
        delete samline;
        *new_fragment = false;
        *seen_a_read = false;
        break;

    case DATA_LINE:
        *new_fragment = *prev_fragment_id != samline->fragment_id && *seen_a_read;
        
        if (*new_fragment && ! sam_buffer.yet_unpaired_entries.empty())
        {
            //all entries should be properly paired when we encounter a new fragment.  violation.
            fprintf(stderr, "Error: input SAM buffer not sorted by read_id, fragment pair\n"
                    "This entry is paired in sequencing but as yet have not found its mate:\n");
            (*sam_buffer.yet_unpaired_entries.begin())->print(stderr, true);
            fprintf(stderr, "\nand this entry, on a different fragment, was found, violating the ordering\n");
            samline->print(stderr, true);
                
            exit(1);
        }

        //check sort order
        if (! (*prev_fragment_id <= samline->fragment_id))
        {
            fprintf(stderr,
                    "Error: input SAM file not sorted by fragment_id\n"
                    "Previous query name: %s\n"
                    "Current query name: %s\n",
                    prev_qname, samline->qname);
            exit(1);
        }

        strcpy(prev_qname, samline->qname);
        *prev_fragment_id = samline->fragment_id;

        CONTIG_OFFSETS::const_iterator contig_iter = sam_buffer.sam_order->contig_offsets.begin();

        samline->SetFlattenedPosition(sam_buffer.sam_order->contig_offsets, & contig_iter);
        std::pair<SamLine const*, bool> 
            insert_result = sam_buffer.insert(samline);

        *seen_a_read = true;

        if (*new_fragment)
        {
            if ((*low_bound) != NULL)
            {
                delete (*low_bound);
            }
            (*low_bound) = new SamLine(*insert_result.first);
        }
        break;
    }
}


typedef std::map<RAW_SCORE_T, size_t> LENGTH_HISTO;
typedef std::map<RAW_SCORE_T, size_t, RawScoreBetter> SCORE_HISTO;
typedef std::map<RAW_SCORE_T, LENGTH_HISTO, RawScoreBetter> LENGTH_HISTO_BY_SCORE; // key: score, value: length histo

//initialize these so that the best scores are at the beginning.


std::vector<size_t>
TallyFragmentLengths(FILE ** sam_fh, 
                     RawScoreBetter const& score_order,
                     SamBuffer & tally_buffer,
                     char const* raw_score_tag,
                     RAW_SCORE_T max_valid_fragment_score,
                     size_t num_frags_to_use)
{

    LENGTH_HISTO_BY_SCORE length_counts = LENGTH_HISTO_BY_SCORE(score_order);
    SCORE_HISTO score_totals(score_order);

    bool new_fragment;
    bool seen_a_read = false;
    size_t prev_fragment_id = 0;
    char prev_qname[1024] = "";
    bool allow_absent_seq_qual = true; // why not?
    size_t template_length;
    size_t num_frags_used = 0;

    std::vector<RAW_SCORE_T> packed_scores;
    std::vector<RAW_SCORE_T>::iterator fit;

    RAW_SCORE_T top_score;
    SamLine const* top_fragment_entry;
    SamLine * low_bound = NULL;

    while (! feof(*sam_fh))
    {

        NextLine(*sam_fh, tally_buffer, allow_absent_seq_qual,
                 &new_fragment, &seen_a_read, prev_qname, &prev_fragment_id, 
                 &low_bound);

        if (new_fragment)
        {
            packed_scores = 
                GetPackedScores(tally_buffer, 0, SIZE_MAX, raw_score_tag, max_valid_fragment_score);

            fit = std::min_element(packed_scores.begin(), packed_scores.end(), score_order);

            if (fit == packed_scores.end()
                || (*fit) == max_valid_fragment_score)
            {
                continue;
            }

            top_score = UnpackScore(*fit, max_valid_fragment_score);

            PAIRED_READ_SET::const_iterator mpit = tally_buffer.unique_entry_pairs.begin();
            std::advance(mpit, std::distance(packed_scores.begin(), fit));
            top_fragment_entry = (*mpit).first;
            
            assert(top_fragment_entry->isize > 0);
            template_length = static_cast<size_t>(top_fragment_entry->isize);
            length_counts[top_score][template_length]++;
            score_totals[top_score]++;
            num_frags_used++;
            tally_buffer.purge(NULL, NULL, NULL, low_bound);
            if (num_frags_to_use < num_frags_used)
            {
                break;
            }
        }
    }
    rewind(*sam_fh);

    //at this point, tally_buffer should be spent

    //select the top-scoring N% of fragments as a heuristic proxy for 'correctly aligned' fragments
    size_t partial_sum = 0;

    SCORE_HISTO::const_iterator sit_end = score_totals.begin();
    LENGTH_HISTO_BY_SCORE::const_iterator lit_end = length_counts.begin();

    while (partial_sum < num_frags_used)
    {
        partial_sum += (*sit_end++).second;
        ++lit_end;
    }

    //compute mean and sd of elite fragment alignments
    std::vector<size_t> top_frag_lengths;
    LENGTH_HISTO_BY_SCORE::const_iterator lit;
    LENGTH_HISTO::const_iterator hit;
    for (lit = length_counts.begin(); lit != lit_end; ++lit)
    {
        for (hit = (*lit).second.begin(); hit != (*lit).second.end(); ++hit)
        {
            std::fill_n(std::back_inserter(top_frag_lengths), (*hit).second, (*hit).first);
        }
    }

    return top_frag_lengths;
}


void QuantileFragmentEstimate(size_t min_allowed_fragment_length,
                              size_t max_allowed_fragment_length,
                              float frag_dist_low_quantile,
                              float frag_dist_high_quantile,
                              FILE ** sam_fh,
                              RawScoreBetter const& score_order,
                              SamBuffer & tally_buffer,
                              char const* raw_score_tag,
                              RAW_SCORE_T max_valid_fragment_score,
                              size_t num_top_scoring_frags_used,
                              size_t * min_est_fragment_length, 
                              size_t * max_est_fragment_length)
{

    if (min_allowed_fragment_length > max_allowed_fragment_length)
    {
        fprintf(stderr, "Minimum allowed fragment length (-n) (%Zu) is greater than Max length (-x) (%Zu)\n",
                min_allowed_fragment_length, max_allowed_fragment_length);
        exit(1);
    }

    if (frag_dist_low_quantile < 0.0
        || frag_dist_low_quantile > 1.0
        || frag_dist_high_quantile < 0.0
        || frag_dist_high_quantile > 1.0
        || frag_dist_low_quantile >= frag_dist_high_quantile)
    {
        fprintf(stderr, "Low quantile (-q) (%f) and high quantile (-Q) (%f) must be in (0,1) and low < high\n",
                frag_dist_low_quantile, frag_dist_high_quantile);
        exit(1);
    }


    std::vector<size_t> top_frag_lengths =
        TallyFragmentLengths(sam_fh, 
                             score_order, 
                             tally_buffer, 
                             raw_score_tag,
                             max_valid_fragment_score,
                             num_top_scoring_frags_used);

    if (top_frag_lengths.empty())
    {
        *min_est_fragment_length = min_allowed_fragment_length;
        *max_est_fragment_length = max_allowed_fragment_length;
    }
    else
    {
        std::vector<size_t>::iterator nth_iter;
        nth_iter = top_frag_lengths.begin() + frag_dist_low_quantile * top_frag_lengths.size();
        std::nth_element(top_frag_lengths.begin(), nth_iter, top_frag_lengths.end());
        *min_est_fragment_length = std::max(*nth_iter, min_allowed_fragment_length);
        
        nth_iter = top_frag_lengths.begin() + frag_dist_high_quantile * top_frag_lengths.size();
        if (nth_iter == top_frag_lengths.end())
        {
            --nth_iter;
        }
        std::nth_element(top_frag_lengths.begin(), nth_iter, top_frag_lengths.end());
        *max_est_fragment_length = std::min(*nth_iter, max_allowed_fragment_length);
    }

}


//Parse and load a score calibration file as produced from 'score_dist'
CollapseScoreTriplet 
ParseScoreCalibration(FILE * score_cal_fh,
                      char * raw_score_tag,
                      bool * larger_score_is_better,
                      int * & score_cpd)
{

    char larger_score_better_char;
    int max_valid_score;

    size_t num_parsed = fscanf(score_cal_fh, 
                               "score_tag: %s\n"
                               "max_valid_fragment_score: %i\n"
                               "larger_score_better: %c\n",
                               raw_score_tag,
                               &max_valid_score,
                               &larger_score_better_char);

    if (num_parsed != 3)
    {
        fprintf(stderr, "Error: score calibration file doesn't have score_tag, "
                "max_valid_fragment_score, or larger_is_better fields.\n"
                "Please produce with 'samutil score_dist'\n");
        exit(1);
    }

    *larger_score_is_better = larger_score_better_char == 'Y';

    RAW_SCORE_T top_raw_score, sec_raw_score, given_raw_score;
    int mapq;

    CollapseScoreTriplet score_triplet(0, max_valid_score + 1);

    size_t max_possible_joint = score_triplet(max_valid_score + 1,
                                              max_valid_score + 1,
                                              max_valid_score + 1);

    score_cpd = new int[max_possible_joint + 1];
    std::fill(score_cpd, score_cpd + max_possible_joint + 1, 0);

    while (! feof(score_cal_fh))
    {
        fscanf(score_cal_fh, "%i\t%i\t%i\t%i%*[^\n]\n",
               &top_raw_score,
               &sec_raw_score,
               &given_raw_score,
               &mapq);

        score_cpd[score_triplet(top_raw_score, sec_raw_score, given_raw_score)] = mapq;
    }

    return score_triplet;
}
