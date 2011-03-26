/*
  Evaluate alignment accuracy of provided SAM file.
  Assumptions:
  1.  SAM file read ids have the guide alignment format:
  id:read1:D:S:P:C[,D:S:P:C]+:read2:D:S:P:C[,D:S:P:C]+

  where
  D = dna name or contig name
  S = strand [+ or -]
  P = start position, (where first base on a contig is 1)
  C = CIGAR string for this read only

  2.  

  Produces:
  A tabular output with:
  chromosome start_position guide_depth correct_placed_depth error_placed_depth

  Policy for multiple alignments of the same read:
  1.  The guide depth will only be tallied upon the first encounter of the read.
  2.  The 'correct' and 'error' depths though, will be tallied on every encounter of the read.

 */

#include <map>
#include <string>
#include <numeric>

#include "cigar_ops.h"
#include "sam_helper.h"
#include "dep/tools.h"
#include "align_eval_raw.h"


int align_eval_raw_usage(size_t qdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval raw [OPTIONS] alignment_sorted.sam jumps.txt align_stats.{full,by_half,cumul}.txt\n\n"
            "Options:\n\n"
            "-q  INT      minimum quality to include alignment[%Zu]\n"
            "-p  FLAG     (primary alignment) if present, require primary alignment.\n"
            "             in this case, do not tally 'correct' and 'error' jumps\n"
            "             for non-primary alignments (SAM flag 0x0100, 256) [absent]\n",
            qdef);

    fprintf(stderr,
            "jumps.txt has fields:\n"
            "<contig> <position> <guide_jump> <correct_jump> <error_jump>\n"
            "It is used as input for other align_eval functions\n\n"
            "Output files:\n\n"
            "align_stats.full.txt:  first_op_length num_bases_correct num_reads_in_category total_num_reads_in_category total_num_guidefields:\n"
            "align_stats.by_halfs.txt:  first_op_length num_bases_correct num_reads_in_category total_num_reads_in_category total_num_guidefields:\n"
            "<guide_cigar> <num_correct_placed_bases> <num_reads>\n"
            "It tallies the number of each primary alignment having the number of \n"
            "correctly placed bases\n");
    return 1;
}

bool remember_read(size_t read_id, std::vector<char> * seen_ids)
{
    size_t read_id_byte = read_id / 8;
    char read_id_mask = 1<<(read_id % 8);

    if ((*seen_ids).size() <= read_id_byte)
    {
        (*seen_ids).resize(read_id_byte * 2, '\0');
    }
    bool seen_this_read = (*seen_ids)[read_id_byte] & read_id_mask;
    (*seen_ids)[read_id_byte] |= read_id_mask;
    return seen_this_read;
}


void CigarFromSimSAMLine(char const* readname, bool first_in_pair,
                         bool is_ones_based_position,
                         read_coords * guide_coords)
{

    int ones_adjust = is_ones_based_position ? 1 : 0;
    char strand;

    size_t num_fields = sscanf(readname, "%zu:", &guide_coords->fragment_id);
    if (num_fields != 1)
    {
        fprintf(stderr, "Error: simulated read has improper name:\n%s\n"
                "Please see 'readsim'\n", readname);
        exit(1);
    }

    if (first_in_pair)
    {
        char const* read1_str = strstr(readname, "read1");
        if (read1_str == NULL)
        {
            fprintf(stderr, "Error: bad format of read name.  See readsim documentation\n");
            exit(1);
        }
        size_t num_fields =
            sscanf(read1_str, "read1:%[^:]:%c:%zu:%[^:]", 
                   guide_coords->contig, 
                   &strand, 
                   &guide_coords->position,
                   guide_coords->cigar);
        assert(num_fields == 4);
    }
    else
    {
        //find string 'read2'
        char const* read2_str = strstr(readname, "read2");
        if (read2_str == NULL)
        {
            fprintf(stderr, "Error: bad format of read name.  See readsim documentation\n");
            exit(1);
        }

        size_t num_fields =
            sscanf(read2_str, "read2:%[^:]:%c:%zu:%[^:]", 
                   guide_coords->contig, 
                   &strand, 
                   &guide_coords->position,
                   guide_coords->cigar);
        assert(num_fields == 4);
    }
    guide_coords->read_id = (guide_coords->fragment_id<<1) + (first_in_pair ? 0 : 1);
    guide_coords->position -= ones_adjust;
    guide_coords->pos_stranded = (strand == '+');
}



int main_align_eval_raw(int argc, char ** argv)
{

    size_t min_mapping_quality_def = 0;
    size_t min_mapping_quality = min_mapping_quality_def;

    bool require_primary_alignment = false;

    char c;
    while ((c = getopt(argc, argv, "q:p")) >= 0)
    {
        switch(c)
        {
        case 'q': min_mapping_quality = static_cast<size_t>(atof(optarg)); break;
        case 'p': require_primary_alignment = true; break;
        default: return align_eval_raw_usage(min_mapping_quality_def); break;
        }
    }
    
    int arg_count = optind + 5;
    if (argc != arg_count)
    {
        return align_eval_raw_usage(min_mapping_quality_def);
    }

    char const* alignment_sam_file = argv[optind];
    char const* output_jumps_file = argv[optind + 1];
    char const* dist_full_file = argv[optind + 2];
    char const* dist_by_halfs_file = argv[optind + 3];
    char const* dist_summary_file = argv[optind + 4];

    FILE * alignment_sam_fh = open_if_present(alignment_sam_file, "r");
    FILE * output_jumps_fh = open_if_present(output_jumps_file, "w");
    FILE * dist_full_fh = open_or_die(dist_full_file, "w", "Output distribution file in full");
    FILE * dist_by_halfs_fh = open_or_die(dist_by_halfs_file, "w", "Output distribution file by halfs");
    FILE * dist_summary_fh = open_or_die(dist_summary_file, "w", "Output distribution summary file");

    bool sam_is_ones_based = true;

    //evaluation should not depend on the actual sequence, just the alignment blocks.
    bool allow_absent_seq_qual = true;

    SamLine * samline;
    std::map<std::string, size_t> contig_lengths = ContigLengths(alignment_sam_fh);
    SetToFirstDataLine(&alignment_sam_fh);

    CONTIG_OFFSETS contig_offsets = ContigOffsets(contig_lengths);

    //stores features ordered by flattened transcript coordinate
    std::map<size_t, Feature> ordered_features;
    std::map<size_t, Feature>::iterator feature_iter;
    std::map<size_t, Feature>::iterator fit, fit_bound;

    size_t primary_alignment_mapq[256];
    std::fill(primary_alignment_mapq, primary_alignment_mapq + 256, 0);

    read_coords guide_coords;

    size_t num_read_id_bytes = 1000000;

    //use for the purpose of recording the guide read id
    std::vector<char> seen_read_ids(num_read_id_bytes, '\0');

    //for recording whether we have used them.  this is difference since
    std::vector<char> used_read_ids(num_read_id_bytes, '\0');

    //boundary position of left-most feature
    size_t old_feature_bound = 0, new_feature_bound;
    size_t feature_bound;

    size_t num_guide_bases = 0; //number of bases that were generated in the simulation
    size_t num_total_bases = 0; //total number of aligned bases in the alignment file
    size_t num_skipped_bases = 0; //skipped due to low quality or non-primary or unaligned
    size_t num_used_bases = 0; //used for tallying correct and error
    size_t num_correct_bases = 0;
    size_t num_error_bases = 0;
    size_t num_trimmed_bases = 0; //bases in guide sequence but absent in alignment

    //the following should be true:
    //1. if primicity is a criterion for being used, then num_used_bases <= num_guide_bases
    //2. num_used_bases = num_correct_bases + num_error_bases;

    size_t this_num_guide_bases;

    OFFSETS_ITER guide_offset_iter;
    OFFSETS_ITER align_offset_iter;

    //some statistics on alignment types.  key is cigar string, value is histogram of 
    //bases correctly placed.
    std::map<size_t, std::vector<size_t> > correct_placed_by_type;
    std::map<size_t, std::vector<size_t> >::iterator correct_placed_iter;
    std::map<size_t, size_t> total_guide_by_type;
    std::map<size_t, size_t> total_used_by_type;

    samline = new SamLine(alignment_sam_fh, sam_is_ones_based, allow_absent_seq_qual);

    while (samline->parse_flag == DATA_LINE)
    {
        CigarFromSimSAMLine(samline->qname, 
                            samline->first_read_in_pair(),
                            sam_is_ones_based,
                            &guide_coords);

        assert(num_used_bases == num_correct_bases + num_error_bases + num_trimmed_bases);
        assert(num_total_bases == num_used_bases + num_skipped_bases);
        assert((! require_primary_alignment) || (num_used_bases <= num_guide_bases));
        assert(! (samline->mapq > 0 && samline->alignment_not_primary()));

        // printf("%s\tID: %Zu\tG: %Zu\tT: %Zu\tU:%Zu\tS: %Zu\tC: %Zu\tE: %Zu\tdiff1: %Zi\tdiff2: %Zi\tdiff3: %Zi\n",
        //        samline->qname, guide_coords.read_id,
        //        num_guide_bases, num_total_bases, num_used_bases, num_skipped_bases,
        //        num_correct_bases, num_error_bases,
        //        num_total_bases - (num_used_bases + num_skipped_bases),
        //        num_used_bases - (num_correct_bases + num_error_bases),
        //        num_guide_bases - num_used_bases);


        //determine low_bound coordinate for read (least coordinate
        //between its guide and aligned positions)
        //tally the guide state
        guide_offset_iter = contig_offsets.find(guide_coords.contig);

        size_t align_flat_low_bound = 
            flattened_position(samline->rname, samline->zero_based_pos(), 
                               samline->flag, contig_offsets,
                               &align_offset_iter);

        assert(guide_offset_iter != contig_offsets.end());
        
        size_t guide_flat_low_bound = (*guide_offset_iter).second + guide_coords.position;

        //any cumulative even before this low bound may be printed and purged
        new_feature_bound = std::min(guide_flat_low_bound, align_flat_low_bound);

        if (old_feature_bound == 0)
        {
            //this is first iteration of the loop
            old_feature_bound = new_feature_bound;
        }

        if (old_feature_bound < new_feature_bound)
        {
            //input is sorted by MIN_ALIGN_GUIDE.  There will be no
            //more updates to any jumps before this.  purge all
            //features before this.  they are complete
            fit_bound = ordered_features.upper_bound(old_feature_bound);

            for (fit = ordered_features.begin(); fit != fit_bound; ++fit)
            {
                Feature const& feature = (*fit).second;
                size_t feature_bound = (*fit).first - (*feature.offset_iter).second;
                char const* contig = (*feature.offset_iter).first;

                if (feature.guide_jumps != 0
                    || feature.correct_jumps != 0
                    || feature.error_jumps != 0)
                {
                    fprintf(output_jumps_fh, "%s\t%Zu\t%i\t%i\t%i\n",
                            contig, feature_bound,
                            feature.guide_jumps, 
                            feature.correct_jumps,
                            feature.error_jumps);
                }
            }

            ordered_features.erase(ordered_features.begin(), fit_bound);
        }
                     
        else if (old_feature_bound == new_feature_bound)
        {
            //do not purge.  update as usual
        }
        else
        {
            //old_feature_bound > new_feature_bound.  reads out of order
            fprintf(stderr, "Error: input reads out of order.  Should be sorted by "
                    "the minimum of guide position and alignment position\n");
            exit(1);
        }
        old_feature_bound = new_feature_bound;

        bool seen_this_read = remember_read(guide_coords.read_id, &seen_read_ids);

        Cigar::CIGAR_VEC guide_cigar = 
            Cigar::FromString(guide_coords.cigar, guide_coords.position);

        size_t first_guide_op_len = guide_cigar[1].length;

        this_num_guide_bases = Cigar::Length(guide_cigar, false);

        num_total_bases += this_num_guide_bases;

        if (! seen_this_read)
        {
            total_guide_by_type[first_guide_op_len]++;

            feature_bound = guide_offset_iter->second;
            for (Cigar::CIGAR_ITER gi = guide_cigar.begin(); gi != guide_cigar.end(); ++gi)
            {
                switch((*gi).op)
                {
                case Cigar::M:
                    {
                        assert(old_feature_bound <= feature_bound);

                        Feature & start_feature = ordered_features[feature_bound];
                        start_feature.guide_jumps++;
                        start_feature.offset_iter = guide_offset_iter;
                        
                        //this could be a bug, if we allow guide reads to
                        //be fusions across contigs. !!!
                        Feature & end_feature = ordered_features[feature_bound + (*gi).length];
                        end_feature.guide_jumps--;
                        end_feature.offset_iter = guide_offset_iter;
                        
                        feature_bound += (*gi).length;
                        num_guide_bases += (*gi).length;

                    }
                    break;
                case Cigar::I:
                    assert(false);
                    break;
                case Cigar::D:
                case Cigar::N:
                case Cigar::H:
                case Cigar::P:
                    feature_bound += (*gi).length;
                    break;
                case Cigar::S:
                case Cigar::None:
                    break;
                }
            }
        }

        if (! samline->alignment_not_primary())
        {
            ++primary_alignment_mapq[samline->mapq];
        }

        if (samline->mapq < min_mapping_quality
            || (require_primary_alignment && samline->alignment_not_primary()))
        {
            num_skipped_bases += this_num_guide_bases;
            delete samline;
            samline = new SamLine(alignment_sam_fh, sam_is_ones_based, allow_absent_seq_qual);
            continue;
        }
        else
        {
            //using this read
            num_used_bases += this_num_guide_bases;
            total_used_by_type[first_guide_op_len]++;
        }
        
        bool used_this_read = remember_read(guide_coords.read_id, &used_read_ids);

        if (require_primary_alignment && used_this_read)
        {
            fprintf(stderr, "Error: requiring primary alignment, but have used this read before:\n");
            samline->print(stderr, false);
            exit(1);
        }
                    
        Cigar::CIGAR_VEC test_cigar = Cigar::FromString(samline->cigar, samline->pos);

        std::multiset<std::pair<size_t, size_t> > guide_cigar_index =
            Cigar::ComputeOffsets(guide_cigar);


        size_t this_num_correct_bases = 0;
        //tally the correct and incorrect features
        if (guide_offset_iter == align_offset_iter &&
            guide_coords.pos_stranded == samline->query_on_pos_strand())
        {
            //alignment is at least on the right contig and strand.
            //appropriate to do the merge.
            Cigar::CIGAR_VEC merge_cigar = 
                Cigar::TransitiveMerge(guide_cigar, guide_cigar_index, test_cigar, true, false);

            size_t guide_read_pos = 0;
            size_t test_read_pos = 0;

            //position on the contig for the first alignment among the test and guide
            feature_bound = (*guide_offset_iter).second;

            //once the merge is complete, there will (in general) be an offset between the guide
            //and test start (represented by an 'I' or 'D' state).  
            for (Cigar::CIGAR_ITER mi = merge_cigar.begin(); mi != merge_cigar.end(); ++mi)
            {
                
                Cigar::Unit const& unit = *mi;

                Feature & start_feature = ordered_features[feature_bound];
                start_feature.offset_iter = guide_offset_iter;
                
                Feature & end_feature = ordered_features[feature_bound + unit.length];
                end_feature.offset_iter = guide_offset_iter;

                switch(unit.op)
                {
                case Cigar::M:
                    assert(old_feature_bound <= feature_bound);
                    
                    if (test_read_pos == guide_read_pos)
                    {
                        start_feature.correct_jumps++;
                        end_feature.correct_jumps--;
                        num_correct_bases += unit.length;
                        this_num_correct_bases += unit.length;
                    }                    
                    else
                    {
                        start_feature.error_jumps++;
                        end_feature.error_jumps--;
                        num_error_bases += unit.length;
                    }
                    
                    feature_bound += unit.length;
                    guide_read_pos += unit.length;
                    test_read_pos += unit.length;

                    break;

                case Cigar::I:
                    //present in test but not in guide.  test is erroneously placed
                    assert(old_feature_bound <= feature_bound);
                    start_feature.error_jumps++;
                    end_feature.error_jumps--;
                    
                    feature_bound += unit.length;
                    test_read_pos += unit.length;
                    num_error_bases += unit.length;

                    break;

                case Cigar::D:
                case Cigar::N:
                    feature_bound += unit.length;
                    guide_read_pos += unit.length;
                    break;
                    
                case Cigar::H:
                    assert(false);
                    break;

                case Cigar::P:
                    feature_bound += unit.length;
                    break;

                case Cigar::S:
                case Cigar::None:
                    //this shouldn't happen since we're dealing with a merged CIGAR
                    assert(false);
                    break;
                }
            }
            /*
              initially i thought we needed to traverse the entire
              guide read alignment (bottom length of guide_cigar, or
              top length of merge_cigar).  but, it turns out, only the
              test read alignment (bottom length of merge_cigar) need
              be fully traversed.  The guide length is handled in the
              previous code, and is parsed once per read id.

              This code is executed once per alignment (multiple
              alignments per read)
            */
            
            //assert(guide_read_pos == samline->raw_read_length());
            assert(test_read_pos == samline->raw_read_length());
            num_trimmed_bases += this_num_guide_bases - test_read_pos;

        }

        else if (align_offset_iter != contig_offsets.end())
        {
            //guide and test not on same contig or strand.  classify
            //as completely wrong
            feature_bound = align_offset_iter->second;

            if (guide_offset_iter == align_offset_iter &&
                guide_coords.pos_stranded != samline->query_on_pos_strand() &&
                guide_coords.position == samline->zero_based_pos())
            {
                //samline->print(stderr, false);
            }

            for (Cigar::CIGAR_ITER ti = test_cigar.begin(); ti != test_cigar.end(); ++ti)
            {
                switch((*ti).op)
                {
                case Cigar::M:
                    {
                        assert(old_feature_bound <= feature_bound);
                        Feature & start_feature = ordered_features[feature_bound];
                        start_feature.error_jumps++;
                        start_feature.offset_iter = align_offset_iter;
                        
                        Feature & end_feature = ordered_features[feature_bound + (*ti).length];
                        end_feature.error_jumps--;
                        end_feature.offset_iter = align_offset_iter;
                        
                        feature_bound += (*ti).length;
                        num_error_bases += (*ti).length;
                    }
                    break;
                case Cigar::I:
                    //doesn't get counted here at least
                    break;
                case Cigar::D:
                case Cigar::N:
                case Cigar::H:
                case Cigar::P:
                    feature_bound += (*ti).length;
                    break;
                case Cigar::S:
                case Cigar::None:
                    break;
                }
            }
            num_trimmed_bases += this_num_guide_bases - Cigar::Length(test_cigar, false);
        }

        else
        {
            //test not mapped at all
            num_skipped_bases += this_num_guide_bases;
        }

        if (correct_placed_by_type.find(first_guide_op_len) == correct_placed_by_type.end())
        {
            correct_placed_by_type[first_guide_op_len] = std::vector<size_t>(this_num_guide_bases + 1);
        }
        correct_placed_by_type[first_guide_op_len][this_num_correct_bases]++;

        assert(num_used_bases == num_correct_bases + num_error_bases + num_trimmed_bases);
        assert(num_total_bases == num_used_bases + num_skipped_bases);
        assert((! require_primary_alignment) || (num_used_bases <= num_guide_bases));

        delete samline;
        samline = new SamLine(alignment_sam_fh, sam_is_ones_based, allow_absent_seq_qual);
    }
    if (samline->parse_flag != END_OF_FILE)
    {
        fprintf(stderr, "Error: samline is not the end of the file:\n%s\n\n",
                samline->line);
        exit(1);
    }

    delete samline;
    close_if_present(alignment_sam_fh);

    //print out remaining jumps
    for (feature_iter = ordered_features.begin();
         feature_iter != ordered_features.end();
         ++feature_iter)
    {
        Feature & feature = (*feature_iter).second;
        size_t feature_bound = (*feature_iter).first 
            - (*feature.offset_iter).second;
            
        char const* contig = (*feature.offset_iter).first;
        
        if (feature.guide_jumps != 0
            || feature.correct_jumps != 0
            || feature.error_jumps != 0)
        {
            fprintf(output_jumps_fh, "%s\t%Zu\t%i\t%i\t%i\n",
                    contig, feature_bound,
                    feature.guide_jumps, 
                    feature.correct_jumps,
                    feature.error_jumps);
        }
    }
    
    close_if_present(output_jumps_fh);

    //delete contig_offsets names
    char const** keys = new char const*[contig_offsets.size()];
    size_t k;
    CONTIG_OFFSETS::iterator ci;
    for (ci = contig_offsets.begin(), k = 0; ci != contig_offsets.end(); ++ci, ++k)
    {
        keys[k] = (*ci).first;
    }
    contig_offsets.clear();

    for (size_t k = 0; k != contig_offsets.size(); ++k)
    {
        delete keys[k];
    }
    delete keys;

    
    fprintf(dist_summary_fh, 
            "G (guide bases)           :%14Zu\n"
            "T (total in alignments)   :%14Zu\n"
            "U (used)                  :%14Zu\n"
            "C (correct)               :%14Zu\n"
            "E (error)                 :%14Zu\n"
            "R (trimmed)               :%14Zu\n"
            "S (skipped)               :%14Zu\n"
            "U/G                       :%14f\n"
            "C/G                       :%14f\n"
            "E/G                       :%14f\n"
            "C/U                       :%14f\n"
            "E/U                       :%14f\n\n"
,
            num_guide_bases, num_total_bases, num_used_bases, 
            num_correct_bases, num_error_bases, num_trimmed_bases, num_skipped_bases,
            100.0 * static_cast<double>(num_used_bases) / static_cast<double>(num_guide_bases),
            100.0 * static_cast<double>(num_correct_bases) / static_cast<double>(num_guide_bases),
            100.0 * static_cast<double>(num_error_bases) / static_cast<double>(num_guide_bases),
            100.0 * static_cast<double>(num_correct_bases) / static_cast<double>(num_used_bases),
            100.0 * static_cast<double>(num_error_bases) / static_cast<double>(num_used_bases));

    size_t q_cumul = 0;
    for (size_t q = 0; q != 256; ++q)
    {
        q_cumul += primary_alignment_mapq[q];
        fprintf(dist_summary_fh, "%Zu\t%Zu\t%Zu\n", q, primary_alignment_mapq[q], q_cumul);
    }

    fprintf(dist_full_fh, "#%s\t%s\t%s\t%s\t%s\n", 
            "first_op_length",
            "num_bases_correctly_placed",
            "num_reads_in_category",
            "total_number_reads_in_category",
            "total_number_guide");
    
    fprintf(dist_by_halfs_fh, "#%s\t%s\t%s\t%s\t%s\n", 
            "first_op_length",
            "num_reads_under_half_correct",
            "num_reads_over_half_correct",
            "total_used_by_type",
            "total_guide_by_type");
    
    for (correct_placed_iter = correct_placed_by_type.begin();
         correct_placed_iter != correct_placed_by_type.end();
         ++correct_placed_iter)
    {
        size_t first_guide_op_len = (*correct_placed_iter).first;
        std::vector<size_t> & dist = (*correct_placed_iter).second;
        // ncb == num_correct_bases
        for (size_t ncb = 0; ncb != dist.size(); ++ncb)
        {
            fprintf(dist_full_fh, "%Zu\t%Zu\t%Zu\t%Zu\t%Zu\n", 
                    first_guide_op_len, ncb, dist[ncb], 
                    total_used_by_type[first_guide_op_len],
                    total_guide_by_type[first_guide_op_len]);
        }
        size_t half = dist.size() / 2;
        size_t under_half = std::accumulate(dist.begin(), dist.begin() + half, 0);
        size_t over_half = std::accumulate(dist.begin() + half, dist.end(), 0);

        fprintf(dist_by_halfs_fh, 
                "%Zu\t%Zu\t%Zu\t%Zu\t%Zu\n", 
                first_guide_op_len, under_half, over_half,
                total_used_by_type[first_guide_op_len],
                total_guide_by_type[first_guide_op_len]);
    }

    fclose(dist_full_fh);
    fclose(dist_by_halfs_fh);
    fclose(dist_summary_fh);

    return 0;
}
