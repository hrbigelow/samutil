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
            fscanf(score_cal_fh, "%i\t%i\t%i\t%*i%*[^\n]\n", &top_raw_score, &sec_raw_score, &given_raw_score);

        if (num_parsed != 3)
        {
            fprintf(stderr, "Error: score calibration file %s has bad format.\n", score_cal_file);
            exit(1);
        }

        this->max_fragment_score = std::max(this->max_fragment_score, top_raw_score);
        this->min_fragment_score = std::min(this->min_fragment_score, given_raw_score);
    }

    this->table_shift_bits = 
        static_cast<size_t>(ceilf(log2f(static_cast<float>(this->max_fragment_score - this->min_fragment_score))));

    size_t max_index = this->score_table_index(this->max_fragment_score,
                                               this->max_fragment_score,
                                               this->max_fragment_score);

    this->score_table = new int[max_index + 1];
    std::fill(this->score_table, this->score_table + max_index, 0);

    rewind(score_cal_fh);
    fscanf(score_cal_fh, "score_tag: %*s\n");

    while (! feof(score_cal_fh))
    {
        int num_parsed = 
            fscanf(score_cal_fh, "%i\t%i\t%i\t%i%*[^\n]\n",
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
        assert(index <= max_index);

        this->score_table[index] = mapq;
        this->larger_score_better = top_raw_score > sec_raw_score;
    }
    fclose(score_cal_fh);

    this->worst_fragment_score = 
        this->larger_score_better 
        ? this->min_fragment_score 
        : this->max_fragment_score;

    this->raw_score_comp = RawScoreComp(this->larger_score_better);


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
    bool has_score;
    int fragment_score;

    int unqual_score = 
        samline->alignment_score(this->raw_score_tag, this->worst_fragment_score, &has_score);

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
    if (a->this_fragment_unmapped())
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
    return index;
}


//order by [raw_score, space]
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
void set_score_fields(SamBuffer const& sam_buffer, 
                      FragmentScore const& fragment_scoring)
{

    SINGLE_READ_SET::iterator pit;

    size_t N = sam_buffer.unique_entries.size();
    int * raw_scores = new int[N];
    int * raw_scores_copy = new int[N + 2];

    FragmentScoreWrap fs_wrap(&fragment_scoring);

    std::map<ScoreSpace, size_t, FragmentScoreWrap> stratum_hist(fs_wrap);
    std::map<ScoreSpace, size_t, FragmentScoreWrap>::iterator sit;

    size_t i = 0;
    for (pit = sam_buffer.unique_entries.begin();
         pit != sam_buffer.unique_entries.end(); ++pit, ++i)
    {
        SamLine * samline = const_cast<SamLine *>((*pit));

        raw_scores[i] = fragment_scoring.raw_score(samline);

        if (samline->alignment_space == AlignSpaceMissing)
        {
            char as = fragment_scoring.alignment_space(samline);
            samline->alignment_space = as;
            char space_string[] = "-";
            space_string[0] = as;

            samline->add_tag("XP", 'A', space_string);
        }

        
        ScoreSpace ss(raw_scores[i], samline->alignment_space);
        sit = stratum_hist.find(ss);
        if (sit == stratum_hist.end())
        {
            sit = stratum_hist.insert(std::make_pair(ss, 0)).first;
        }
        (*sit).second++;
    }

    std::copy(raw_scores, raw_scores + N, raw_scores_copy);
    raw_scores_copy[N] = fragment_scoring.worst_fragment_score;
    raw_scores_copy[N + 1] = fragment_scoring.worst_fragment_score;

    std::nth_element(raw_scores_copy, raw_scores_copy + N, 
                     raw_scores_copy + N + 2, 
                     fragment_scoring.raw_score_comp);

    std::sort(raw_scores_copy + N, raw_scores_copy + N + 2, fragment_scoring.raw_score_comp);

    int new_mapq;
    int top_fragment_score = raw_scores_copy[N + 1];
    int sec_fragment_score = raw_scores_copy[N];
    int fragment_score;
    char fragment_space;
    char rank_tag_value[10];
    char size_tag_value[10];

    size_t num_strata = stratum_hist.size();

    bool first_encountered_top_stratum = true;

    for (i = 0, pit = sam_buffer.unique_entries.begin();
         pit != sam_buffer.unique_entries.end(); ++pit, ++i)
    {
        SamLine * samline = const_cast<SamLine *>((*pit));

        fragment_score = raw_scores[i];
        fragment_space = samline->alignment_space;

        size_t table_index = 
            fragment_scoring.score_table_index(top_fragment_score, sec_fragment_score, fragment_score);

        new_mapq = fragment_scoring.score_table[table_index];

        bool mapq_correct = samline->this_fragment_unmapped() ? (new_mapq == 0) : true;

        assert(mapq_correct);

        sit = stratum_hist.find(ScoreSpace(fragment_score, fragment_space));
        assert(sit != stratum_hist.end());

        size_t index = std::distance(stratum_hist.begin(), sit);
        size_t stratum_size = (*sit).second;
        size_t stratum_rank = num_strata - index;

        //set mapq
        samline->mapq = new_mapq;

        sprintf(rank_tag_value, "%zu", stratum_rank);
        sprintf(size_tag_value, "%zu", stratum_size);

        //set rank tag
        samline->add_tag("XY", 'i', rank_tag_value);
        samline->add_tag("XZ", 'i', size_tag_value);

        
        //set primary flag
        if (stratum_rank == 1 && samline->all_fragments_mapped() && first_encountered_top_stratum)
        {
            //alignment is primary.
            samline->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
            first_encountered_top_stratum = false;
        }
        else if (samline->this_fragment_unmapped())
        {
            //alignment is 'primary' per downstream weird interpretations
            samline->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
        }
        else
        {
            samline->flag |= SamFlags::ALIGNMENT_NOT_PRIMARY;
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
        && samline->this_fragment_on_pos_strand() == guide_coords.pos_stranded)
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
        if (! sam_buffer.incomplete_entries.empty())
        {
            //all entries should be properly paired when we encounter a new fragment.  violation.
            fprintf(stderr, "Error: Reached end of file, yet this entry is paired in sequencing"
                    " but as yet have not found its mate:\n");
            (*sam_buffer.incomplete_entries.begin())->print_sam(stderr);
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
        
        if (*new_fragment && ! sam_buffer.incomplete_entries.empty())
        {
            //all entries should be properly paired when we encounter a new fragment.  violation.
            fprintf(stderr, "Error: input SAM buffer not sorted by read_id, fragment pair\n"
                    "This entry is incomplete:\n");
            (*sam_buffer.incomplete_entries.begin())->print_sam(stderr);
            fprintf(stderr, "\nand this entry, on a different fragment, was found, violating the ordering\n");
            samline->print_sam(stderr);
                
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


// typedef std::map<int, size_t> LENGTH_HISTO;
// typedef std::map<int, size_t, RawScoreBetter> SCORE_HISTO;
// typedef std::map<int, LENGTH_HISTO, RawScoreBetter> LENGTH_HISTO_BY_SCORE; // key: score, value: length histo

//initialize these so that the best scores are at the beginning.


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
