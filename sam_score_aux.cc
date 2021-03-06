#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "sam_score_aux.h"
#include "cigar_ops.h"
#include "dep/tools.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <climits>





FragmentScore::~FragmentScore()
{
    if (this->space_ordering != NULL)
    {
        delete this->space_ordering;
        this->space_ordering = NULL;
    }

    if (this->score_table != NULL)
    {
        delete this->score_table;
        this->score_table = NULL;
    }

    for (CONTIG_SPACE::iterator it = this->contig_space.begin();
         it != this->contig_space.end(); ++it)
    {
        delete (*it).first;
    }
}


void FragmentScore::init(char const* score_cal_file, char const* contig_space_file)
{
    FILE * score_cal_fh = 
        open_or_die(score_cal_file, "r", "Input score calibration file");

    size_t num_parsed = 
        fscanf(score_cal_fh, "score_tag: %s\n", this->raw_score_tag);

    if (num_parsed != 1)
    {
        fprintf(stderr, 
                "Error: score calibration file %s doesn't have score_tag field."
                "Please produce with 'samutil score_dist'\n", score_cal_file);
        exit(1);
    }


    int top_raw_score, sec_raw_score, given_raw_score;
    int mapq;

    this->max_fragment_score = -INT_MAX;
    this->min_fragment_score = INT_MAX;

    while (! feof(score_cal_fh))
    {
        int num_parsed =
            fscanf(score_cal_fh, "%i\t%i\t%i\t%*i\n", &top_raw_score, &sec_raw_score, &given_raw_score);

        if (num_parsed != 3)
        {
            fprintf(stderr, "Error: score calibration file %s has bad format.\n", score_cal_file);
            exit(1);
        }

        this->max_fragment_score = std::max(this->max_fragment_score, 
                                            std::max(top_raw_score, sec_raw_score));

        this->min_fragment_score = std::min(this->min_fragment_score, 
                                            std::min(top_raw_score, sec_raw_score));
    }

    this->table_shift_bits = 
        static_cast<size_t>(ceilf(log2f(static_cast<float>(this->max_fragment_score - this->min_fragment_score))));

    // score_table_index uses max_index, so must first set it to a permissive value
    this->max_index = SIZE_MAX;
    this->max_index = this->score_table_index(this->max_fragment_score,
                                              this->max_fragment_score,
                                              this->max_fragment_score) + 1;

    this->score_table = new int[this->max_index];
    std::fill(this->score_table, this->score_table + this->max_index, 0);

    rewind(score_cal_fh);
    fscanf(score_cal_fh, "score_tag: %*s\n");

    while (! feof(score_cal_fh))
    {
        int num_parsed = 
            fscanf(score_cal_fh, "%i\t%i\t%i\t%i\n",
                   &top_raw_score,
                   &sec_raw_score,
                   &given_raw_score,
                   &mapq);

        if (num_parsed != 4)
        {
            fprintf(stderr, "Error: score calibration file %s has bad format.\n", score_cal_file);
            exit(1);
        }

        size_t index = this->score_table_index(top_raw_score, sec_raw_score, given_raw_score);
        assert(index <= this->max_index);

        this->score_table[index] = mapq;
        this->larger_score_better = top_raw_score > sec_raw_score;
    }
    fclose(score_cal_fh);

    this->worst_fragment_score = 
        this->larger_score_better 
        ? this->min_fragment_score 
        : this->max_fragment_score;

    //want the best score at the beginning
    bool ascending_order = ! this->larger_score_better;

    this->raw_score_comp = RawScoreComp(ascending_order);

    FILE * contig_space_fh = 
        open_or_die(contig_space_file, "r", "Input contig alignment space file");

    char first_char;

    while (! feof(contig_space_fh))
    {
        if ((first_char = fgetc(contig_space_fh)) == '#')
        {
            //ignore comment line
            char * dummy = NULL;
            size_t dummy_size;
            getline(&dummy, &dummy_size, contig_space_fh);
            delete dummy;
        }
        else
        {
            //actual data line
            ungetc(first_char, contig_space_fh);
            this->space_ordering = new char[256];
            fscanf(contig_space_fh, "%s\n", this->space_ordering);
            break;
        }
    }

    char contig_name[256];
    char space_code;
    while (! feof(contig_space_fh))
    {
        fscanf(contig_space_fh, "%s\t%1c\n", contig_name, &space_code);
        char * this_contig_name = new char[strlen(contig_name) + 1];
        strcpy(this_contig_name, contig_name);
        this->contig_space.insert(std::make_pair(this_contig_name, space_code));
    }

    fclose(contig_space_fh);
}


int FragmentScore::raw_score(SamLine const* samline) const
{
    assert(samline->flag.is_rsam_format);

    int fragment_score;

    int unqual_score = samline->tags.raw_score_present 
        ? samline->tags.raw_score
        : this->worst_fragment_score;

    size_t fragment_length = Cigar::Length(Cigar::FromString(samline->cigar, 0), false);

    if (fragment_length >= this->min_fragment_length
        && fragment_length <= this->max_fragment_length
        && unqual_score >= this->min_fragment_score
        && unqual_score <= this->max_fragment_score)
    {
        //fragment score is valid
        fragment_score = unqual_score;
    }
    else
    {
        //outside the range.  give it the worst value
        fragment_score = this->worst_fragment_score;
    }

    return fragment_score;
}

char FragmentScore::alignment_space(SamLine const* a) const
{
    if (a->flag.this_fragment_unmapped
        || (a->flag.is_rsam_format && (! a->flag.all_fragments_mapped)))
    {
        return this->space_ordering[0];
    }

    CONTIG_SPACE::const_iterator cs_iter;
    cs_iter = this->contig_space.find(a->rname);
    
    if (cs_iter == this->contig_space.end())
    {
        fprintf(stderr, "Error: no Alignment Space definition "
                "for this reference contig: %s\n",
                a->rname);
        exit(1);
    }
    return (*cs_iter).second;
}


//retrieve the index in the score table corresponding to the entry
//packed in a single size_t by the offset
size_t FragmentScore::score_table_index(int top_score, int sec_score, int given_score) const
{
    size_t index = top_score - this->min_fragment_score;
    index = index<<this->table_shift_bits;
    index += sec_score - this->min_fragment_score;
    index = index<<this->table_shift_bits;
    index += given_score - this->min_fragment_score;

    // censor the index if it is out-of-bounds
    return std::min(index, this->max_index);
}


//order by [raw_score, space] with best score at the beginning
bool FragmentScore::operator()(ScoreSpace const& a, ScoreSpace const& b) const
{
    return this->raw_score_comp(a.raw_score, b.raw_score)
        || (a.raw_score == b.raw_score
            && (index(this->space_ordering, a.alignment_space) <
                index(this->space_ordering, b.alignment_space)));
}



FragmentScoreWrap::FragmentScoreWrap(FragmentScore const* _fs) : fragment_score(_fs) { }

bool FragmentScoreWrap::operator()(ScoreSpace const& a, ScoreSpace const& b) const
{
    return (*this->fragment_score)(a, b);
}


//Set mapq, primary alignment flag, XP, XY, and XZ tags
void set_score_fields(SINGLE_READ_SET::iterator beg,
                      SINGLE_READ_SET::iterator end,
                      FragmentScore const& fragment_scoring)
{

    SINGLE_READ_SET::iterator pit;

    size_t N = std::distance(beg, end);
    int * raw_scores = new int[N];
    int * raw_scores_copy = new int[N + 2];

    FragmentScoreWrap fs_wrap(&fragment_scoring);

    std::map<ScoreSpace, size_t, FragmentScoreWrap> stratum_hist(fs_wrap);
    std::map<ScoreSpace, size_t, FragmentScoreWrap>::iterator sit;

    size_t i = 0;
    for (pit = beg; pit != end; ++pit, ++i)
    {
        SamLine * samline = const_cast<SamLine *>((*pit));

        assert(samline->flag.is_rsam_format);

        raw_scores[i] = fragment_scoring.raw_score(samline);

        if (! samline->tags.alignment_space_present)
        {
            char as = fragment_scoring.alignment_space(samline);
            samline->tags.alignment_space = as;
            samline->tags.alignment_space_present = true;
        }

        
        ScoreSpace ss(raw_scores[i], samline->tags.alignment_space);
        sit = stratum_hist.find(ss);
        if (sit == stratum_hist.end())
        {
            sit = stratum_hist.insert(std::make_pair(ss, 0)).first;
        }
        (*sit).second++;
    }

    //extract the top two raw scores, or substitute with missing defaults
    std::partial_sort_copy(raw_scores, raw_scores + N,
                           raw_scores_copy, raw_scores_copy + N,
                           fragment_scoring.raw_score_comp);
    
    int * raw_unique_end = std::unique(raw_scores_copy, raw_scores_copy + N);
    
    *raw_unique_end = fragment_scoring.worst_fragment_score;
    *++raw_unique_end = fragment_scoring.worst_fragment_score;
    
    int top_fragment_score = raw_scores_copy[0];
    int sec_fragment_score = raw_scores_copy[1];

    int new_mapq;
    int fragment_score;
    char fragment_space;

    // size_t num_strata = stratum_hist.size();

    bool first_encountered_top_stratum = true;

    for (i = 0, pit = beg; pit != end; ++pit, ++i)
    {
        SamLine * samline = const_cast<SamLine *>((*pit));

        fragment_score = raw_scores[i];
        fragment_space = samline->tags.alignment_space;

        size_t table_index = 
            fragment_scoring.score_table_index(top_fragment_score, sec_fragment_score, fragment_score);

        new_mapq = fragment_scoring.score_table[table_index];

        bool mapq_correct = samline->flag.this_fragment_unmapped ? (new_mapq == 0) : true;

        assert(mapq_correct);

        sit = stratum_hist.find(ScoreSpace(fragment_score, fragment_space));
        assert(sit != stratum_hist.end());

        size_t index = std::distance(stratum_hist.begin(), sit);

        samline->tags.stratum_size = (*sit).second;
        samline->tags.stratum_size_present = true;

        samline->tags.stratum_rank = index + 1;
        samline->tags.stratum_rank_present = true;

        //set mapq
        samline->mapq = new_mapq;

        //set primary flag
        if (samline->tags.stratum_rank == 1 
            && samline->flag.all_fragments_mapped 
            && first_encountered_top_stratum)
        {
            //alignment is primary.
            samline->flag.alignment_not_primary = 0;
            first_encountered_top_stratum = false;
        }
        else if (! samline->flag.all_fragments_mapped)
        {
            //alignment is 'primary' per downstream weird interpretations
            samline->flag.alignment_not_primary = 0;
        }
        else
        {
            samline->flag.alignment_not_primary = 1;
        }

    }
    
    delete raw_scores;
    delete raw_scores_copy;
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
        && samline->flag.this_fragment_on_neg_strand != guide_coords.pos_stranded)
    {
        size_t expanded_samline_pos;
        //proceed to measure actual base overlap
        Cigar::CIGAR_VEC test_cigar = Cigar::FromString(samline->cigar, 0);

        Cigar::CIGAR_VEC merge_cigar = 
            Cigar::Expand(guide_coords.blocks, test_cigar, samline->zero_based_pos(), 
                          &expanded_samline_pos, true);

        merge_cigar.insert(merge_cigar.begin(), Cigar::Unit(Cigar::Ops[Cigar::D], expanded_samline_pos));

        num_correct_bases = 
            Cigar::CountAlignedPositions(merge_cigar, num_bases_guide, num_bases_test);
    }
    return num_correct_bases;
}


score_rsam_alloc_binary::score_rsam_alloc_binary(FragmentScore const* _fragment_scoring,
                                                 size_t _arll) : 
    fragment_scoring(_fragment_scoring),
    average_rsam_line_length(_arll) { }


std::vector<char> *
score_rsam_alloc_binary::operator()(std::pair<SAMIT, SAMIT> const& range,
                                    SamBuffer * sam_buffer)
{
    if (range.first == range.second)
    {
        return new std::vector<char>();
    }

    SAMIT cur;

    for (cur = range.first; cur != range.second; ++cur)
    {
        InsertResult ins = sam_buffer->insert(*cur);
        if (! ins.was_inserted)
        {
            delete ins.remaining_entry;
        }
    }
    // At this point, the buffer should have joined everything together
    if (! sam_buffer->incomplete_entries.empty())
    {
        fprintf(stderr, "Error: at a fragment boundary but there are still incomplete SAM entries\n");
        exit(1);
    }
    // score them.
    SINGLE_READ_SET::iterator beg, sit, end;
    beg = sam_buffer->unique_entries.begin();

    for (end = sam_buffer->unique_entries.begin(); 
         end != sam_buffer->unique_entries.end(); ++end)
    {
        if ((*beg)->fragment_id != (*end)->fragment_id)
        {
            set_score_fields(beg, end, *this->fragment_scoring);
            beg = end;
        }
    }
    set_score_fields(beg, end, *this->fragment_scoring);

    size_t num_entries = sam_buffer->unique_entries.size();
    std::vector<char> * result = new std::vector<char>;
    result->reserve(num_entries * this->average_rsam_line_length);
    
    char tmp_buffer[1024];

    for (sit = sam_buffer->unique_entries.begin(); 
         sit != sam_buffer->unique_entries.end(); ++sit)
    {
        (*sit)->sprint(tmp_buffer);
        delete (*sit);
        result->insert(result->end(), tmp_buffer, tmp_buffer + strlen(tmp_buffer));
    }
    sam_buffer->unique_entries.clear();

    return result;
}


parse_sam_unary::parse_sam_unary() { }

SamLine * parse_sam_unary::operator()(char * sam_string)
{
    SamLine * rec = new SamLine(sam_string);

    if (rec->parse_flag == PARSE_ERROR)
    {
        fprintf(stderr, "Encountered formatting error in 'parse_sam_unary'.\n");
        exit(1);
    }
    return rec;
}


set_flattened_pos_unary::set_flattened_pos_unary(SamOrder const* _sam_order) :
    sam_order(_sam_order)
{
}

void set_flattened_pos_unary::operator()(SamLine * rec)
{
    // hack.  should really be semi-sequential for max efficiency
    CONTIG_OFFSETS::const_iterator contig_iter = 
        this->sam_order->contig_offsets.begin();

    rec->SetFlattenedPosition(this->sam_order->contig_offsets, 
                              & contig_iter);
}

/*
std::vector<size_t>
TallyFragmentLengths(FILE ** sam_fh, 
                     RawScoreBetter const& score_order,
                     SamBuffer & tally_buffer,
                     char const* raw_score_tag,
                     int max_valid_fragment_score,
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

    std::vector<int> packed_scores;
    std::vector<int>::iterator fit;

    int top_score;
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

            PAIRED_READ_SET::const_iterator mpit = tally_buffer.unique_entries.begin();
            std::advance(mpit, std::distance(packed_scores.begin(), fit));
            top_fragment_entry = (*mpit).first;
            
            assert(top_fragment_entry->tlen > 0);
            template_length = static_cast<size_t>(top_fragment_entry->tlen);
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
                              int max_valid_fragment_score,
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
*/
