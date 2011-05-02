#include <algorithm>
#include <cstdio>

#include "cigar_ops.h"
#include "sam_raw_score_aux.h"

ScorePair::ScorePair(size_t _t, size_t _s) : top_score(_t), sec_score(_s) { }
ScorePair::ScorePair() : top_score(0), sec_score(0) { }

bool ScorePair::operator<(ScorePair const& s) const
{
    return this->top_score < s.top_score
        || (this->top_score == s.top_score
            && this->sec_score < s.sec_score);
}


size_t FragmentScore(size_t first_raw_score, 
                     size_t second_raw_score,
                     size_t missing_score_default,
                     int template_length,
                     size_t min_allowed_template_length,
                     size_t max_allowed_template_length)
{
    size_t fragment_score;
    if (template_length >= min_allowed_template_length 
        && template_length <= max_allowed_template_length
        && first_raw_score != missing_score_default
        && second_raw_score != missing_score_default)
    {
        fragment_score = first_raw_score + second_raw_score;
    }
    else
    {
        fragment_score = missing_score_default;
    }
    return fragment_score;
}



//assume sam_buffer has alignments of only one physical
//fragment in its 'unique_entry_pairs' buffer.
//calculate all fragment_scores, and return the top two
ScorePair FindTopTwoScores(SamBuffer const& sam_buffer, 
                           size_t min_allowed_fragment_length,
                           size_t max_allowed_fragment_length,
                           char const* tag,
                           size_t missing_score_default,
                           bool larger_is_better,
                           std::vector<size_t> * fragment_scores)
{
                           
    (*fragment_scores).clear();
    PAIRED_READ_SET::iterator pit;
    bool first_has_score;
    bool second_has_score;
    ScorePair score_pair;

    for (pit = sam_buffer.unique_entry_pairs.begin();
         pit != sam_buffer.unique_entry_pairs.end(); ++pit)
    {
        SamLine const* first = ((*pit).first);
        SamLine const* second = ((*pit).second);
        size_t fragment_score = 
            FragmentScore(first->alignment_score(tag, missing_score_default, &first_has_score), 
                          second->alignment_score(tag, missing_score_default, &second_has_score),
                          missing_score_default,
                          first->isize,
                          min_allowed_fragment_length,
                          max_allowed_fragment_length);

        (*fragment_scores).push_back(fragment_score);
    }
            
    //3. retrieve relevant table
    assert(! (*fragment_scores).empty());

            
    if ((*fragment_scores).size() == 1)
    {
        score_pair = ScorePair((*fragment_scores)[0], missing_score_default);
    }
    else
    {
        // std::binary_function<size_t, size_t, bool> * comp = 
        //     larger_is_better ? &std::greater<size_t>() : &std::less<size_t>();

        std::vector<size_t> fragment_scores_copy((*fragment_scores));

        if (larger_is_better)
        {
            std::nth_element(fragment_scores_copy.begin(), 
                             fragment_scores_copy.begin() + 2,
                             fragment_scores_copy.end(),
                             std::greater<size_t>());

            std::sort(fragment_scores_copy.begin(),
                      fragment_scores_copy.begin() + 2,
                      std::greater<size_t>());
        }
        else
        {
            std::nth_element(fragment_scores_copy.begin(), 
                             fragment_scores_copy.begin() + 2,
                             fragment_scores_copy.end(),
                             std::less<size_t>());

            std::sort(fragment_scores_copy.begin(),
                      fragment_scores_copy.begin() + 2, 
                      std::less<size_t>());
        }
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
        Cigar::CIGAR_VEC test_cigar = Cigar::FromString(samline->cigar, samline->pos);
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
              bool ones_based_pos,
              bool allow_absent_seq_qual,
              bool * new_fragment, 
              bool * ignore_sambuffer_bound,
              bool * seen_a_read,
              char * prev_qname,
              size_t * prev_qid)
{
    SamLine const* samline = new SamLine(unscored_sam_fh, ones_based_pos, allow_absent_seq_qual);

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
            bool did_advance = sam_buffer.safe_advance_lowbound(samline);
        }
        break;
    }
}
