#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <unistd.h>

#include <cstdlib>
#include <algorithm>
#include <getopt.h>

#include "sam_score_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"
#include "file_utils.h"

#include <omp.h>
#include <parallel/algorithm>

int score_mapq_usage(size_t ldef, size_t Ldef, size_t mdef)
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
            "calibration.qcal: a histogram over the set of alignment categories\n"
            "(top score, 2nd score, given score)\n"
            "tallying number of alignments that are correct or incorrect.\n\n"

            "contig_space.txt: defines contig alignment space for each contig in 'unscored.rsort.sam' file\n"
            "Format:\n\n"

            "# comments ...\n"
            "TG\n"
            "chr1<tab>G\n"
            "chr2<tab>G\n"
            "...\n"
            "ENST000005299<tab>T\n"
            "ENST899998402<tab>T\n"
            "...\n\n"

            "unscored_rsort.sam: alignment file sorted by (read id / pair flag) and\n"
            "having alignment score tags given in option -s.\n\n"

            "scored.rsort.sam:    output alignment file sorted by (read id / pair flag)\n"
            "with mapq field updated.\n"
            "mapq will reflect Phred-scaled probability of correct alignment.\n\n",
            ldef, Ldef, mdef
            );
    return 1;
}



int main_score_mapq(int argc, char ** argv)
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
        default: return score_mapq_usage(ldef, Ldef, mdef); break;
        }
    }

    if (argc != optind + 4)
    {
        return score_mapq_usage(ldef, Ldef, mdef);
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
    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(unscored_sam_fh);
    sam_order.AddHeaderContigStats(unscored_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt, expected_read_layout,
                            fragment_scoring.raw_score_tag,
                            fragment_scoring.worst_fragment_score);

    PrintSAMHeader(&unscored_sam_fh, scored_sam_fh);
    fflush(scored_sam_fh);

    std::vector<char *> sam_lines;
    std::vector<char *>::iterator sit;

    size_t nbytes_unused;
    size_t nbytes_read;

    // at any one point, we will have:  a chunk_buffer, a buffer of converted chunks,
    size_t chunk_size = max_mem / 2;
    char * chunk_buffer_in = new char[chunk_size + 1];

    while (! feof(unscored_sam_fh))
    {
        nbytes_read = fread(chunk_buffer_in, 1, chunk_size, unscored_sam_fh);
        chunk_buffer_in[nbytes_read] = '\0';
        sam_lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & nbytes_unused);

        if (nbytes_unused > 0)
        {
            //only do this if we need to. this resets a flag that affects feof()
            fseek(unscored_sam_fh, -1 * static_cast<off_t>(nbytes_unused), SEEK_CUR);
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
        pre_fragment_id = nbytes_unused > 0 ? (*revit)->fragment_id : SIZE_MAX;
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

        // now recycle unused part of buffer by replacing all newlines,
        // and then rewinding the file by that length
        off_t unused_length = 0;
        size_t unused_index = std::distance(sam_records.begin(), revit.base());
        for (sit = sam_lines.begin() + unused_index; sit != sam_lines.end(); ++sit)
        {
            unused_length += strlen((*sit)) + 1;
        }
        if (unused_length > 0)
        {
            fseek(unscored_sam_fh, - unused_length, SEEK_CUR);
        }

        // join SAM to rSAM, allocate string buffer, print, delete rSAM record
        std::vector<char *> scored_sam_records(ranges.size());
        score_rsam_alloc_binary rsam_print_aux(& fragment_scoring);
        __gnu_parallel::transform(ranges.begin(), ranges.end(), sam_buffers.begin(),
                                  scored_sam_records.begin(), rsam_print_aux);

        for (size_t r = 0; r != ranges.size(); ++r)
        {
            fputs(scored_sam_records[r], scored_sam_fh);
            delete [] scored_sam_records[r];
            delete sam_buffers[r];
        }

    }

    delete [] chunk_buffer_in;

    fclose(unscored_sam_fh);
    fclose(scored_sam_fh);
    return 0;
}
