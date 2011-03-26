/*
  Parse a SAM file and its corresponding pileup file.
  Tally a table of per-base counts partitioned by:
  consensus_score
  quality_score
  consensus_match
  alignment_score
  read_position
  strand

  From this internal table, output marginal statistics to answer certain questions:

  1.  What is the overall alignment score distribution? (alignment_score)
  2.  What is the estimated accuracy of quality scores?
  (fraction of bases that match consensus at each quality score,
  filtered by alignment_score and consensus_score)
  3.  Quality score distribution at each read position?
 */


/*
  Pseudo code:
  Initialize longest_read from user input.
  Initialize a table of counts to zeros.  Table is a flat array of RP x CM x CS x QS x ST x AS.
  Parse a window* of pileup lines, retaining the consensus, consensus_score, position.
  For each SAM line:
     Cache alignment_score, consensus_score.
     Update pointers into consensus base and consensus_score buffers
     For each base position:
        Calculate consensus_match / mismatch
        Tally the count for this base.

  The SAM file must be sorted by read alignment position.  So, we will
  always know from the position in the SAM file what our lowest
  position is, and by inference to the longest read, the highest
  position.  So, implement an on-demand buffering, in which every time
  the requested position is beyond that of what's in the buffer, advance the buffer.

  We need windowed reading of the file though, since we need access to
  the window containing the current read at least.
 */

#include <cstdio>
#include <cstring>

#include "sam_stats_aux.h"
#include "sam_stats.h"

int main_usage()
{
    fprintf(stderr,
            "\nUsage:\n\nsam_stats raw [OPTIONS] input.sam input.pileup > raw_stats.out\n"
            "\nsam_stats out [OPTIONS] raw_stats out_prefix\n\n");
    return 1;
}


int raw_usage(int mdef, int cdef, int qdef, int adef, int rdef, int bdef)
{
    fprintf(stderr, 
            "\nUsage: sam_stats raw [OPTIONS] input.sam input.pileup > raw_stats.out\n"
            "Options:\n\n"
            "-s FLAG  Skip prescan of sam and pileup files to find parameter settings [%i]\n"
            "-c INT   Number of consensus scores (highest consensus score in pileup file [%i])\n"
            "-q INT   Number of quality scores (highest quality score in reads [%i])\n"
            "-a INT   Number of distinct alignment scores (alignment scores will be binned) [%i]\n"
            "-r INT   Number of distinct read positions [%i]\n"
            "-b INT   Number of bases for consensus buffer [%i]\n"
            "\n",
            mdef, cdef, qdef, adef, rdef, bdef
            );
    return 1;
}


int out_usage()
{
    fprintf(stderr,
            "\nUsage: sam_stats out raw_stats > stats.txt\n"
            );
    return 1;
}


int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        return main_usage();
    }
    else if (strcmp(argv[1], "raw") == 0)
    {
        return main_raw(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "out") == 0)
    {
        return main_out(argc - 1, argv + 1);
    }
}
