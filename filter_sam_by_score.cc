#include <cstdio>
#include <cstring>
#include <algorithm>


/* Filter a SAM file, dropping all but the best alignment(s) for each read
   If there are two or more equal-scoring alignments, retain them.

   Assume the alignments of each read are grouped together.

   With the filtered file, one can compute the pileup and run the composition estimation
*/

int main(int argc, char ** argv)
{

    if (argc != 2 && argc != 3)
    {
        printf("\nUsage:\n\n"
               "filter_sam_by_score best.sam [rest.sam] < input.sam\n\n");
        return 0;
    }

    FILE * best_sam = fopen(argv[1], "w");
    FILE * rest_sam = NULL;
    if (argc == 3)
    {
        rest_sam = fopen(argv[2], "w");
    }
        
    size_t const MAX_LINE_LENGTH = 10000;
    size_t const MAX_READ_ALIGNMENTS = 1000;

    char * line_buffer_flat = new char[MAX_LINE_LENGTH * MAX_READ_ALIGNMENTS];
    char ** line = new char *[MAX_READ_ALIGNMENTS];
    int * alignment_score = new int[MAX_READ_ALIGNMENTS];

    char read_id[1000];
    char group_read_id[1000];
    group_read_id[0] = '\0';

    size_t next_li = 0;

    for (size_t li = 0; li != MAX_READ_ALIGNMENTS; ++li)
    {
        line[li] = line_buffer_flat + (li * MAX_LINE_LENGTH);
    }
    
    while(fgets(line[next_li++], MAX_LINE_LENGTH, stdin) != NULL)
    {
        if (line[next_li - 1][0] == '@')
        {
            fprintf(best_sam, "%s", line[next_li - 1]);
            if (rest_sam != NULL)
            {
                fprintf(rest_sam, "%s", line[next_li - 1]);
            }

            next_li = 0;
            continue;
        }

        sscanf(line[next_li - 1], "%s\t", read_id);
        if (group_read_id[0] == '\0')
        {
            strcpy(group_read_id, read_id);
        }

        if (strcmp(read_id, group_read_id) == 0)
        {
            //do nothing, already stored the next line
        }
        else
        {
            //process the set of lines, finding max alignment score
            --next_li; //we read one past.  back up

            for (size_t li = 0; li != next_li; ++li)
            {
                char * score_string = strstr(line[li], "AS:i:");
                sscanf(score_string, "AS:i:%i", &alignment_score[li]);
            }
            int max_alignment_score = 
                *std::max_element(alignment_score, alignment_score + next_li);

            for (size_t li = 0; li != next_li; ++li)
            {
                if (alignment_score[li] == max_alignment_score)
                {
                    fprintf(best_sam, "%s", line[li]);
                }
                else
                {
                    if (rest_sam != NULL)
                    {
                        fprintf(rest_sam, "%s", line[li]);
                    }
                }
            }
            strcpy(line[0], line[next_li]);
            next_li = 1;
            strcpy(group_read_id, read_id);
        }
    }
    
    delete line_buffer_flat;
    delete line;
    delete alignment_score;

    fclose(best_sam);
    if (rest_sam != NULL)
    {
        fclose(rest_sam);
    }

    return 0;
}
