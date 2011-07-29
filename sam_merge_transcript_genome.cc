
#include "sam_buffer.h"
#include "sam_order.h"
#include "dep/tools.h"

int merge_tg_usage()
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil merge_tg [OPTIONS] hybrid_alignments.asort.sam hybrid_alignments.merged.asort.sam\n\n"
            "Options:\n\n"
            "\n"
            "hybrid_alignments.asort.sam is the union of a set of unique genomic\n"
            "alignment records (no tag), and unique genome-projected transcript\n"
            "alignments (tag XP:A:T), all sorted by alignment position\n"
            "\n"
            "samutil merge_tg goes through, finding records that are duplicate alignments\n"
            "and outputs a single record with an updated tag XP:A:M\n\n"
            );
    return 1;
}


int main_merge_tg(int argc, char ** argv)
{
    char c;

    while ((c = getopt(argc, argv, "")) >= 0)
    {
        switch(c)
        {
        default: return merge_tg_usage(); break;
        }
    }

    //printf("optind: %i\n", optind);
    if (argc != 2 + optind)
    {
        return merge_tg_usage();
    }

    char * input_sam_file = argv[optind];
    char * output_sam_file = argv[optind + 1];

    FILE * input_sam_fh = open_or_die(input_sam_file, "r", "Input hybrid alignment file");
    FILE * output_sam_fh = open_or_die(output_sam_file, "w", "Output merged hybrid alignment file");

    bool paired_reads_are_same_stranded = false; // !!! shoudln't this be an option?

    // input SAM file is allowed to have or not have seq / qual.
    // basically, they are just payload, and unaffected by the transformation
    bool allow_absent_seq_qual = true;

    bool numeric_start_fragment_ids = false; //do not assume we have numeric starting IDS
    SamLine::SetGlobalFlags(numeric_start_fragment_ids);

    SamOrder input_sam_order(SAM_RID_POSITION, "READ_ID_FLAG");
    input_sam_order.AddHeaderContigStats(input_sam_fh);

    SamBuffer input_buffer(&input_sam_order, paired_reads_are_same_stranded);
    SamLine * samline;

    char fake_samline[1024];

    size_t prev_pos_index = 0;
    size_t cur_pos_index = 0;

    size_t num_lines_printed = 0;
    size_t num_lines_omitted = 0;

    std::pair<SamLine const*, bool> insert_result;
    
    PrintSAMHeader(&input_sam_fh, output_sam_fh);

    while (! feof(input_sam_fh))
    {
        // parse samline
        samline = new SamLine(input_sam_fh, allow_absent_seq_qual);

        switch (samline->parse_flag)
        {
        case END_OF_FILE: 
            delete samline;
            break;

        case HEADER:
            delete samline;
            break;

        case PARSE_ERROR:
            fprintf(stderr, "Parse error in input hybrid alignment file %s", input_sam_file);
            exit(1);
            break;

        case DATA_LINE:

            //check well-orderedness of input
            sprintf(fake_samline, "%s\t%i\t%s\t%Zu", samline->qname, samline->flag, samline->rname, samline->pos);
            cur_pos_index = (input_buffer.sam_order->*(input_buffer.sam_order->sam_index))(fake_samline);
            if (cur_pos_index < prev_pos_index)
            {
                fprintf(stderr, "Error: input is out of order.\n"
                        "Please sort by 'READ_ID_FLAG' using align_eval sort\n");
                exit(1);
            }
            prev_pos_index = cur_pos_index;

            insert_result = input_buffer.insert(samline);

            if (insert_result.second)
            {
                //entry was successfully inserted. do nothing for the time being.
                ++num_lines_printed;
                char tag_type;
                char tag_value[10];
                if (! insert_result.first->has_tag("XP", tag_type, tag_value))
                {
                    SamLine * variable_samline = const_cast<SamLine *>(insert_result.first);
                    variable_samline->add_tag("XP", 'A', "G");
                }
                else
                {
                    assert(strcmp(tag_value, "T") == 0);
                }
            }
            else
            {
                //entry is a duplicate. update the tag of the pre-existing entry
                //this entry represents the transcriptome and genome entries merged (M)
                SamLine * variable_samline = const_cast<SamLine *>(insert_result.first);
                variable_samline->add_tag("XP", 'A', "M");
                ++num_lines_omitted;
            }

            assert(insert_result.first != NULL);
            input_buffer.safe_advance_lowbound(insert_result.first);
            input_buffer.purge(output_sam_fh, NULL, NULL, false);
        }
    }
    input_buffer.purge(output_sam_fh, NULL, NULL, true);

    fclose(input_sam_fh);
    fclose(output_sam_fh);

    fprintf(stderr, "%zu lines printed, %zu duplicate lines omitted.\n", num_lines_printed, num_lines_omitted);
    return 0;
}
