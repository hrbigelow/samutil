/*
  output a subset of lines of <query.txt> whose selected <fields>
  match <select.txt>.  Assume query.txt and select.txt are sorted
  according to these fields.

  Does this also assume that <select.txt> has lines that are unique
  w.r.t the fields?
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "cisortho/string_tools.h"
#include "line_tools.h"


const size_t QUERY_MAX_LINE = 1000000;
const size_t GUIDE_MAX_LINE = 1000000;


bool parse_check_line(FILE * fh, char * line, const size_t max_length)
{
    char * ret = fgets(line, max_length, fh);
    if (ret == NULL)
    {
        if (feof(fh))
        {
            return false;
        }
        else
        {
            fprintf(stderr, "error parsing line\n");
            exit(1);
        }
    }
    else
    {
        return true;
    }
}


int main(int argc, char **argv)
{
    if (argc != 6)
    {
        fprintf(stderr, 
                "Usage: filter_matching_lines <guide> <guide_fmt> <query> <query_fmt> <field_types>\n\n"
                "<guide> provides lines that are unique w.r.t the fields defined by <guide_fmt>\n"
                "<query> may contain duplicate lines w.r.t the fields defined by <quiery_fmt>\n"
                "Both <guide> and <query> are assumed sorted.\n"
                "<field_types> are a string of [sif]+, i.e. si, or sisff\n"
                );
        exit(1);
    }

    char * guide_fname = argv[1];
    char * guide_fmt_raw = argv[2];

    char * query_fname = argv[3];
    char * query_fmt_raw = argv[4];

    char * field_types = argv[5];

    FILE * query_fh = fopen(query_fname, "r");
    FILE * guide_fh = fopen(guide_fname, "r");

    char query_fmt[100];
    char guide_fmt[100];

    convert_escapes(query_fmt, query_fmt_raw);
    convert_escapes(guide_fmt, guide_fmt_raw);

    size_t num_fields = strlen(field_types);

    char * buffer = new char[num_fields * 200];

    void ** query_fields = new void *[num_fields];
    void ** guide_fields = new void *[num_fields];

    for (size_t f = 0; f != num_fields; ++f)
    {
        query_fields[f] = buffer + (100 * f);
        guide_fields[f] = buffer + (100 * (f + num_fields));
    }

    char query_line[QUERY_MAX_LINE];
    char guide_line[GUIDE_MAX_LINE];


    //initialization
    fgets(query_line, QUERY_MAX_LINE, query_fh);
    LineTools::parse_fields(query_fmt, query_line, field_types, query_fields);

    fgets(guide_line, GUIDE_MAX_LINE, guide_fh);
    LineTools::parse_fields(guide_fmt, guide_line, field_types, guide_fields);

    while (1)
    {
        int g_vs_q = LineTools::compare_all(guide_fields, query_fields, 
                                            field_types, num_fields);
        
        if (g_vs_q < 0)
        {
            //parse guide line to catch up.
            if (! parse_check_line(guide_fh, guide_line, GUIDE_MAX_LINE))
            {
                break;
            }
            LineTools::parse_fields(guide_fmt, guide_line, field_types, guide_fields);
        }
        else if (g_vs_q == 0)
        {
            //print query
            printf("%s", query_line);
            if (! parse_check_line(query_fh, query_line, QUERY_MAX_LINE))
            {
                break;
            }
            LineTools::parse_fields(query_fmt, query_line, field_types, query_fields);
        }
        else
        {
            //no matching guide.  catch up with query
            if (! parse_check_line(query_fh, query_line, QUERY_MAX_LINE))
            {
                break;
            }
            LineTools::parse_fields(query_fmt, query_line, field_types, query_fields);
        }
    }

    delete buffer;
    delete query_fields;
    delete guide_fields;
}
