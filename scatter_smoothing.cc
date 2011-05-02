/*
expect a file with:
average_value value
average_value value
...

somewhere in the input.

*/

#include <map>
#include <cstring>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>
#include "file_utils.h"
#include "dep/tools.h"
#include "cisortho/string_tools.h"

int main(int argc, char ** argv)
{
    char const* infile = argv[1];
    char const* outfile = argv[2];
    float diff_threshold = atof(argv[3]);
    char const* scanf_line_raw = argv[4];
    size_t max_mem = static_cast<size_t>(atof(argv[5]));

    char scanf_line[256];

    convert_escapes(scanf_line, scanf_line_raw);

    size_t const max_line = 100000;

    gzFile in_zfh = gzopen(infile, "r");
    FILE * out_fh = open_or_die(outfile, "w", "Output file");

    std::multimap<float, std::pair<float, char *> > index;
    typedef std::multimap<float, std::pair<float, char *> >::const_iterator ITER;

    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(in_zfh, max_mem, max_line);
    gzclose(in_zfh);

    if (chunk_lengths.size() > 1)
    {
        fprintf(stderr, "Error: range_join must load entire file in memory.\n"
                "Please launch with larger max_mem argument\n");
        exit(1);
    }

    FILE * in_fh = open_or_die(infile, "r", "Input file");

    char * buffer = new char[chunk_lengths[0] + 1];

    size_t nbytes_read = fread(buffer, 1, chunk_lengths[0], in_fh);
    assert(nbytes_read == chunk_lengths[0]);
    buffer[nbytes_read] = '\0';

    std::vector<char *> lines = FileUtils::find_line_starts(buffer);

    //convert to 
    char * line_end = buffer;
    while ((line_end = strchr(line_end, '\n')) != NULL)
    {
        *line_end++ = '\0';
    }

    float x, y;
    for (size_t i = 0; i != lines.size(); ++i)
    {
        sscanf(lines[i], scanf_line, &x, &y);
        index.insert(std::make_pair(x, std::make_pair(y, lines[i])));
    }

    ITER first_iter, second_iter;
    for (first_iter = index.begin(); first_iter != index.end(); ++first_iter)
    {
        float x_center = (*first_iter).first;
        ITER inner_low = index.lower_bound(x_center - diff_threshold);
        ITER inner_high = index.upper_bound(x_center + diff_threshold);

        size_t n = std::distance(inner_low, inner_high);
        double * val_buffer = new double[n];
        double * wgt_buffer = new double[n];

        //summarize the data here.
        size_t vbi = 0;
        for (second_iter = inner_low; second_iter != inner_high; ++second_iter)
        {
            float x = (*second_iter).first;
            float y = (*second_iter).second.first;
            val_buffer[vbi] = y - x;
            wgt_buffer[vbi] = gsl_ran_gaussian_pdf(x - x_center, diff_threshold);

            ++vbi;
        }
        double wsd = gsl_stats_wsd(wgt_buffer, 1, val_buffer, 1, n);
        fprintf(out_fh, "%f\t%s\n", wsd, (*first_iter).second.second);
        delete val_buffer;
        delete wgt_buffer;
        
    }
    delete buffer;

    fclose(in_fh);
    fclose(out_fh);

    return 0;
}
