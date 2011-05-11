#include <algorithm>
#include <cstdio>
#include <cmath>

#include "cigar_ops.h"
#include "sam_score_aux.h"


AlignmentSpace ParseAlignmentSpace(char as)
{
    for (size_t i = 0; i != strlen(AlignmentSpaces); ++i)
    {
        if (as == AlignmentSpaces[i])
        {
            return static_cast<AlignmentSpace>(i);
        }
    }
    fprintf(stderr, "ParseAlignmentSpace: Error: Don't know this space: %c\n", as);
    exit(1);
}

char PrintAlignmentSpace(AlignmentSpace const& as)
{
    return AlignmentSpaces[as];
}

SpaceScore::SpaceScore(size_t _s, AlignmentSpace _as) : score(_s), space(_as) { }
SpaceScore::SpaceScore() : score(0), space(GENOME) { }

bool SpaceScore::operator==(SpaceScore const& b) const
{
    return this->score == b.score
        && this->space == b.space;
}


// this is only to be used for packing SpaceScore objects into containers
// the order does NOT have to do with 'better' or 'worse'.
// For that, see 'SpaceScore_better'
bool SpaceScore::operator<(SpaceScore const& b) const
{
    return this->score < b.score
        || (this->score == b.score
            && this->space < b.space);
}


SpaceScore_better::SpaceScore_better(bool _lscore, bool _lspace) :
    larger_score_better(_lscore), larger_space_better(_lspace) { }


bool SpaceScore_better::operator()(SpaceScore const& a, SpaceScore const& b) const
{
    int mulscore = this->larger_score_better ? 1 : -1;
    int mulspace = this->larger_space_better ? 1 : -1;

    int score_cmp = a.score < b.score ? mulscore * -1 : (a.score == b.score ? 0 : mulscore * 1);
    int space_cmp = a.space < b.space ? mulspace * -1 : (a.space == b.space ? 0 : mulspace * 1);

    return score_cmp > 0 || (score_cmp == 0 && space_cmp > 0);
}



ScorePair::ScorePair(SpaceScore _t, SpaceScore _s) : top_score(_t), sec_score(_s) { }
ScorePair::ScorePair() : top_score(SpaceScore(0, GENOME)), sec_score(SpaceScore(0, GENOME)) { }


bool ScorePair::operator<(ScorePair const& s) const
{
    SpaceScore_better space_score_better(false, false);
    return space_score_better(s.top_score, this->top_score)
        || (this->top_score == s.top_score
            && space_score_better(s.sec_score, this->sec_score));
}


SpaceScore FragmentScore(size_t first_raw_score, 
                         size_t second_raw_score,
                         AlignmentSpace alignment_space,
                         size_t missing_default_score,
                         size_t template_length,
                         size_t min_allowed_template_length,
                         size_t max_allowed_template_length)
{
    size_t fragment_score;
    if (template_length >= min_allowed_template_length
        && template_length <= max_allowed_template_length
        && first_raw_score != missing_default_score
        && second_raw_score != missing_default_score)
    {
        fragment_score = first_raw_score + second_raw_score;
    }
    else
    {
        fragment_score = missing_default_score;
    }
    return SpaceScore(fragment_score, alignment_space);
}



//assume sam_buffer has alignments of only one physical
//fragment in its 'unique_entry_pairs' buffer.
//calculate all fragment_scores, and return the top two
ScorePair FindTopTwoScores(SamBuffer const& sam_buffer, 
                           size_t min_allowed_fragment_length,
                           size_t max_allowed_fragment_length,
                           char const* tag,
                           size_t missing_default_score,
                           SpaceScore_better const& score_best_is_less,
                           ScoreVec * fragment_scores)
{
                           
    (*fragment_scores).clear();
    PAIRED_READ_SET::iterator pit;
    bool first_has_score;
    bool second_has_score;
    ScorePair score_pair;
    bool has_space;

    for (pit = sam_buffer.unique_entry_pairs.begin();
         pit != sam_buffer.unique_entry_pairs.end(); ++pit)
    {
        SamLine const* first = ((*pit).first);
        SamLine const* second = ((*pit).second);
        
        AlignmentSpace alignment_space = 
            ParseAlignmentSpace(first->alignment_space(AlignmentSpaceMissingDefault, &has_space));

        AlignmentSpace alignment_space2 = 
            ParseAlignmentSpace(second->alignment_space(AlignmentSpaceMissingDefault, &has_space));

        assert(alignment_space == alignment_space2);

        assert(first->isize >= 0);

        SpaceScore fragment_score = 
            FragmentScore(first->alignment_score(tag, missing_default_score, &first_has_score), 
                          second->alignment_score(tag, missing_default_score, &second_has_score),
                          alignment_space,
                          missing_default_score,
                          static_cast<size_t>(first->isize),
                          min_allowed_fragment_length,
                          max_allowed_fragment_length);

        if (fragment_score.score == 10
            && fragment_score.space == TRANSCRIPTOME)
        {
            first->print(stdout, false);
        }
        
        (*fragment_scores).push_back(fragment_score);
    }
    
    //3. retrieve relevant table
    assert(! (*fragment_scores).empty());

            
    if ((*fragment_scores).size() == 1)
    {
        score_pair = 
            ScorePair((*fragment_scores)[0], 
                      SpaceScore(missing_default_score, 
                                 ParseAlignmentSpace(AlignmentSpaceMissingDefault)));
    }
    else
    {
        std::vector<SpaceScore> fragment_scores_copy((*fragment_scores));

        std::nth_element(fragment_scores_copy.begin(), 
                         fragment_scores_copy.begin() + 2,
                         fragment_scores_copy.end(),
                         score_best_is_less);
        
        std::sort(fragment_scores_copy.begin(),
                  fragment_scores_copy.begin() + 2,
                  score_best_is_less);

        score_pair = ScorePair(fragment_scores_copy[0],
                               fragment_scores_copy[1]);
    }

    return score_pair;
}



//count the number of correctly aligned bases between this samline and its guide position
size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index)
{
    size_t num_correct_bases = 0;
    if (strcmp(samline->rname, guide_coords.contig) == 0
        && samline->query_on_pos_strand() == guide_coords.pos_stranded)
    {
        //proceed to measure actual base overlap
        Cigar::CIGAR_VEC test_cigar = Cigar::FromString(samline->cigar, samline->zero_based_pos());
        Cigar::CIGAR_VEC merge_cigar = 
            Cigar::TransitiveMerge(guide_cigar, guide_cigar_index, 
                                   test_cigar, true, false);
        
        num_correct_bases = Cigar::CountAlignedPositions(merge_cigar);
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
              bool * ignore_sambuffer_bound,
              bool * seen_a_read,
              char * prev_qname,
              size_t * prev_qid)
{
    SamLine const* samline = new SamLine(unscored_sam_fh, allow_absent_seq_qual);

    switch (samline->parse_flag)
    {
    case END_OF_FILE: 
        delete samline;
        *ignore_sambuffer_bound = true; //last time around; purge everything.
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
        break;

    case DATA_LINE:
        *new_fragment = ((SamLine::numeric_start_fragment_ids
                          ? *prev_qid != samline->qid
                          : strcmp(prev_qname, samline->qname) != 0)
                         && *seen_a_read);
        
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
        if (! (SamLine::numeric_start_fragment_ids
               ? *prev_qid <= samline->qid
               : strcmp(prev_qname, samline->qname) <= 0))
        {
            if (SamLine::numeric_start_fragment_ids)
            {
                fprintf(stderr, "Error: input SAM file not sorted by integer field of read_id\n"
                        "(and you are running in integer read id mode)\n");
            }
            else
            {
                fprintf(stderr, "Error: input SAM file not sorted by ascii-interpreted read_id.\n"
                        "Hint: if read ids are sorted by integer beginning, use option -n.\n");
            }
            fprintf(stderr, 
                    "Previous query name: %s\n"
                    "Current query name: %s\n",
                    prev_qname, samline->qname);
            exit(1);
        }

        strcpy(prev_qname, samline->qname);
        *prev_qid = samline->qid;

        sam_buffer.insert(samline);

        *seen_a_read = true;

        if (new_fragment)
        {
            sam_buffer.safe_advance_lowbound(samline);
        }
        break;
    }
}


typedef std::map<size_t, size_t> LENGTH_HISTO;
typedef std::map<SpaceScore, size_t, SpaceScore_better> SCORE_HISTO;
typedef std::map<SpaceScore, LENGTH_HISTO, SpaceScore_better> LENGTH_HISTO_BY_SCORE; // key: score, value: length histo

//initialize these so that the best scores are at the beginning.


std::vector<double>
TallyFragmentLengths(FILE ** sam_fh, 
                     SpaceScore_better const& score_best_is_less,
                     SamBuffer & tally_buffer,
                     char const* score_tag,
                     size_t missing_default_score,
                     float fraction_top_scoring_frags_used)
{

    LENGTH_HISTO_BY_SCORE length_counts = LENGTH_HISTO_BY_SCORE(score_best_is_less);
    SCORE_HISTO score_totals(score_best_is_less);

    bool new_fragment;
    bool ignore_sambuffer_bound;
    bool seen_a_read = false;
    size_t prev_qid = 0;
    char prev_qname[1024] = "";
    bool allow_absent_seq_qual = true; // why not?
    size_t template_length;
    size_t total_fragments;
    ScoreVec fragment_scores;
    ScoreVec::iterator fit;

    while (! feof(*sam_fh))
    {

        NextLine(*sam_fh, tally_buffer, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            ScorePair score_pair = 
                FindTopTwoScores(tally_buffer, 0, SIZE_MAX, score_tag, missing_default_score, 
                                 score_best_is_less, & fragment_scores);

            //finds the best score, (i.e. the 'min_element', since we
            //are using the best-is-less logic
            fit = min_element(fragment_scores.begin(), fragment_scores.end(),
                              score_best_is_less);

            SpaceScore max_score = *fit;
            PAIRED_READ_SET::const_iterator mpit = tally_buffer.unique_entry_pairs.begin();
            std::advance(mpit, std::distance(fragment_scores.begin(), fit));

            assert((*mpit).first->isize >= 0);
            template_length = static_cast<size_t>((*mpit).first->isize);
            length_counts[max_score][template_length]++;
            score_totals[max_score]++;
            total_fragments++;
            tally_buffer.purge(NULL, NULL, NULL, ignore_sambuffer_bound);
        }
    }
    rewind(*sam_fh);

    //at this point, tally_buffer should be spent

    //select the top-scoring N% of fragments as a heuristic proxy for 'correctly aligned' fragments
    size_t partial_sum = 0;
    size_t min_num_used_fragments = 
        static_cast<size_t>(floor(fraction_top_scoring_frags_used * static_cast<float>(total_fragments)));

    SCORE_HISTO::const_iterator sit_end = score_totals.begin();
    LENGTH_HISTO_BY_SCORE::const_iterator lit_end = length_counts.begin();

    while (partial_sum < min_num_used_fragments)
    {
        partial_sum += (*sit_end++).second;
        ++lit_end;
    }

    //compute mean and sd of elite fragment alignments
    std::vector<double> elite_fragment_lengths;
    LENGTH_HISTO_BY_SCORE::const_iterator lit;
    LENGTH_HISTO::const_iterator hit;
    for (lit = length_counts.begin(); lit != lit_end; ++lit)
    {
        for (hit = (*lit).second.begin(); hit != (*lit).second.end(); ++hit)
        {
            std::fill_n(std::back_inserter(elite_fragment_lengths), (*hit).second, 
                        static_cast<double>((*hit).first));
        }
    }

    return elite_fragment_lengths;
}
