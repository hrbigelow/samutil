#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <unistd.h>

#include <cstdlib>
#include <algorithm>
#include <getopt.h>

#include "sam_score_aux.h"
#include "sam_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"
#include "file_utils.h"

#include <omp.h>
#include <parallel/algorithm>

int score_usage(size_t ldef, size_t Ldef, size_t mdef)
{
    fprintf(stderr,
            "\n\nUsage:\n\n"
            "samutil score [OPTIONS] calibration.qcal contig_space.txt unscored.rsort.sam scored.rsort.sam\n\n"
            "Options:\n\n"
            "-l     INT     min allowed fragment length for paired alignment [%Zu]\n"
            "-L     INT     max allowed fragment length for paired alignment [%Zu]\n"
            "-y     STRING  expected read layout. If parsing traditional SAM, this is required. \n"
            "               If parsing rSAM, this must be absent or blank. []\n"
            "-m     INT     number bytes of memory to use [%Zu]\n"
            "-t     INT     number of threads to use [1]\n"
            "-C     STRING  work in the directory named here [.]\n"
            "\n\n"
            "calibration.qcal: a table defining what mapq to assign.\n"
            "\n"
            "Format:\n"
            "\n"
            "score_tag: NM\n"
            "0<tab>1<tab>0<tab>255\n"
            "0<tab>2<tab>0<tab>255\n"
            "1<tab>2<tab>1<tab>255\n"
            "1<tab>3<tab>1<tab>255\n"
            "2<tab>3<tab>2<tab>255\n"
            "2<tab>4<tab>2<tab>255\n"
            "3<tab>4<tab>3<tab>255\n"
            "3<tab>5<tab>3<tab>255\n"
            "...\n"
            "\n"
            "First line defines the score tag to consult for raw score calculation\n"
            "Remaining lines have fields:\n"
            "top score, 2nd score, given score, assigned mapq.\n"
            "\n"
            "1) combinations encountered in the data that are missing from this table\n"
            "   are given mapq of zero.\n"
            "\n"
            "2) The logic of lower-is-better or higher-is-better is determined simply by\n"
            "   whether the top score in any row is higher or lower than the 2nd score\n"
            "\n"
            "3) These scores are currently calculated as the sum of score_tag values of\n"
            "   each read in a fragment. Missing tags, out-of-bounds score_tag value, or\n"
            "   out-of-bounds fragment length all result in assignment of the worst raw\n"
            "   score that is seen in the file\n"
            "\n"
            "contig_space.txt: Defines contig space ordering and contig assignments\n"
            "\n"
            "first non-comment line defines preference for contig-space ordering\n"
            "in the 2-D stratification scheme. The remainder of lines assign each\n"
            "contig to one of those spaces\n"
            "\n"
            "Format:\n"
            "\n"
            "# comments ...\n"
            "TG\n"
            "chr1<tab>G\n"
            "chr2<tab>G\n"
            "...\n"
            "ENST000005299<tab>T\n"
            "ENST899998402<tab>T\n"
            "...\n"
            "\n"
            "Assigns mapq to each alignment record based on calibration.qcal\n"
            "\n"
            "Assigns XY:i (stratum rank) and XZ:i (stratum size) tags to all records\n"
            "based on a 2-D stratification scheme (raw score, alignment space).\n"
            "jointly defined by 'calibration_.qcal' and 'contig_space.txt'\n"
            "\n"
            "unscored_rsort.sam: alignment file sorted by FRAGMENT (see samutil sort) and\n"
            "having alignment score tags mentioned in first line of calibration.qcal file.\n"
            "\n"
            "scored.rsort.sam:    output alignment file sorted by (read id / pair flag)\n"
            "with mapq field updated.\n"
            "mapq will reflect Phred-scaled probability of correct alignment.\n\n",
            ldef, Ldef, mdef
            );
    return 1;
}



int main_score(int argc, char ** argv)
{
    char c;

    size_t ldef = 0;
    size_t Ldef = 1000;
    size_t mdef = 1024l * 1024l * 1024l * 4l; // 4 GB memory

    size_t min_fragment_length = ldef;
    size_t max_fragment_length = Ldef;
    size_t max_mem = mdef;

    size_t num_threads = 1;

    char const* expected_read_layout = "";

    char const* working_dir = ".";

    while ((c = getopt(argc, argv, "l:L:y:m:t:C:")) >= 0)
    {
        switch(c)
        {
        case 'l': min_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'L': max_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'y': expected_read_layout = optarg; break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        default: return score_usage(ldef, Ldef, mdef); break;
        }
    }

    if (argc != optind + 4)
    {
        return score_usage(ldef, Ldef, mdef);
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);


    char * score_calibration_file = argv[optind];
    char * contig_space_file = argv[optind + 1];
    char * unscored_sam_file = argv[optind + 2];
    char * scored_sam_file = argv[optind + 3];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * scored_sam_fh = open_or_die(scored_sam_file, "w", "Output scored sam file");

    FragmentScore fragment_scoring(min_fragment_length, max_fragment_length);
    fragment_scoring.init(score_calibration_file, contig_space_file);

    SamOrder sam_order(SAM_RID_POSITION, "NONE");

    char * header_buf = ReadAllocSAMHeader(unscored_sam_fh);
    size_t header_length = strlen(header_buf);
    fwrite(header_buf, 1, header_length, scored_sam_fh);
    fflush(scored_sam_fh);

    sam_order.AddHeaderContigStats(header_buf);
    
    delete header_buf;

    std::vector<char *> sam_lines;
    std::vector<char *>::iterator sit;

    size_t nbytes_fragment = 0;
    size_t nbytes_unused = 0;
    size_t nbytes_read;
    char * last_fragment;

    
    
    // at any one point, we will have:  a chunk_buffer, a buffer of converted chunks,
    size_t chunk_size = max_mem / 9;
    char * chunk_buffer_in = new char[chunk_size + 1];
    char * read_pointer = chunk_buffer_in;


    while (! feof(unscored_sam_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, unscored_sam_fh);
        read_pointer[nbytes_read] = '\0';
        sam_lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        nbytes_fragment = strlen(last_fragment);

        if (! sam_order.Initialized() && ! sam_lines.empty())
        {
            sam_order.InitFromID(sam_lines[0]);
            SamLine::SetGlobalFlags(QNAMEFormat(sam_lines[0]), 
                                    expected_read_layout,
                                    fragment_scoring.raw_score_tag,
                                    fragment_scoring.worst_fragment_score,
                                    false);
        }
        
        std::vector<SamLine *> sam_records(sam_lines.size());
        __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                  sam_records.begin(), 
                                  parse_sam_unary());

        __gnu_parallel::for_each(sam_records.begin(), sam_records.end(), 
                                 set_flattened_pos_unary(& sam_order));

        size_t min_range = 1000; // arbitrarily set the single work unit size to 1000
        size_t est_num_work_units = 1 + (sam_lines.size() / min_range);
        std::vector<std::pair<SAMIT, SAMIT> > ranges;
        ranges.reserve(est_num_work_units);

        std::vector<SamBuffer *> sam_buffers;
        sam_buffers.reserve(est_num_work_units);

        SAMIT pre = sam_records.begin();
        SAMIT cur;
        size_t pre_fragment_id = (*pre)->fragment_id;
        size_t skip = 0;

        for (cur = pre; cur != sam_records.end(); ++cur, ++skip)
        {
            if (skip >= min_range)
            {
                if ((*cur)->fragment_id != pre_fragment_id)
                {
                    ranges.push_back(std::make_pair(pre, cur));
                    pre = cur;
                    pre_fragment_id = (*pre)->fragment_id;
                    skip = 0;
                }
            }
            pre_fragment_id = (*cur)->fragment_id;
        }
        

        // put in one last range by scanning backwards.
        // remember: reverse_iterator(i).base() == i and &*ri == &*(ri.base() - 1).
        // the revit.base() will serve as the proper bound.
        SAMVEC::reverse_iterator revit = sam_records.rbegin();

        // if this is the last chunk, consume the entire range,
        // otherwise, look for a fragment boundary
        pre_fragment_id = feof(unscored_sam_fh) ? SIZE_MAX : (*revit)->fragment_id;
        for (revit = sam_records.rbegin(); revit.base() != pre; ++revit)
        {
            if ((*revit)->fragment_id != pre_fragment_id)
            {
                ranges.push_back(std::make_pair(pre, revit.base()));
                break;
            }
        }

        sam_buffers.resize(ranges.size());
        for (size_t i = 0; i != ranges.size(); ++i)
        {
            sam_buffers[i] = new SamBuffer(& sam_order, false);
        }
        
        // restore newlines into chunk_buffer_in
        // mark start of unused portion
        size_t unused_index = std::distance(sam_records.begin(), revit.base());
        char * unused_start = *(sam_lines.begin() + unused_index);
        for (sit = sam_lines.begin() + unused_index; sit != sam_lines.end(); ++sit)
        {
            assert(*((*sit) + strlen(*sit)) == '\0');
            *((*sit) + strlen(*sit)) = '\n';
        }
        nbytes_unused = strlen(unused_start);

        // join SAM to rSAM, allocate string buffer, print, delete rSAM record
        size_t n_sampled_records = 1000;
        size_t average_rsam_line_length = 
            average_line_length(sam_records.begin(), sam_records.end(), n_sampled_records);

        std::vector<std::vector<char> *> scored_sam_records(ranges.size());
        score_rsam_alloc_binary rsam_print_aux(& fragment_scoring, average_rsam_line_length);

        __gnu_parallel::transform(ranges.begin(), ranges.end(), sam_buffers.begin(),
                                  scored_sam_records.begin(), rsam_print_aux);

        for (size_t r = 0; r != ranges.size(); ++r)
        {
            fwrite(&(*scored_sam_records[r])[0], 1, 
                   (*scored_sam_records[r]).size(), 
                   scored_sam_fh);

            delete scored_sam_records[r];
            delete sam_buffers[r];
        }

        memmove(chunk_buffer_in, unused_start, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

    }

    delete [] chunk_buffer_in;

    fclose(unscored_sam_fh);
    fclose(scored_sam_fh);
    return 0;
}
