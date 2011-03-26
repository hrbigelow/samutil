#include "line_tools.h"
#include "cisortho/string_tools.h"

/*
  join two sorted files on multiple fields as specified by a list of join fields


 */
size_t const MAX_LINE = 1000000;
size_t const MAX_FIELD_STRING = 100;
size_t const MAX_FILES = 9;
size_t const MAX_JOIN_FIELDS = 10;
size_t const MAX_FORMAT_STRING = 1000;
size_t const MAX_FIELD_BYTES = 1000;

int multijoin_usage()
{
    fprintf(stderr,
            "\nUsage: multijoin [options] file1 file2 ... > outfile\n"
            "Options: \n\n"
            "-1 INT_LIST    list of field numbers from \E[1mfile1\E[0m to join on\n"
            "-2 INT_LIST    list of field numbers from \E[1mfile2\E[0m to join on\n"
            "...\n"
            "-n INT_LIST    list of field numbers from \E[1mfilen\E[0m to join on\n"

            "-o STRING      output field specifier string (e.g. 2a3b1a, see note)\n\n"
            "Note: output specifier string, e.g. 2a3b1a is a string of [0-9]+[a-z][0-9]+[a-z]\n");
    exit(1);
}

int main(int argc, char ** argv)
{

    //i.e. -1 "2,4,5" -2 "3,4,2" -3 "7,2,3"
    char field_string[MAX_FILES][MAX_FIELD_STRING];

    char c;
    while ((c = getopt(argc, argv, "1:2:3:4:5:6:7:8:9:o:")) >= 0)
    {
        switch(c)
        {
        case '1': 
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
            char cc[2]; 
            cc[0] = c; 
            cc[1] = '\0';
            size_t fn = static_cast<size_t>(atof(cc) - 1);
            strcpy(fields_string[fn], optarg); 
            break;
        case 'o': strcpy(output_spec, optarg); break;
        default: return multijoin_usage(); break;
        }
    }

    size_t num_input_files = argc - optind;
    if (num_input_files < 2)
    {
        return multijoin_usage();
    }

    //tabular input file names
    char * file[MAX_FILES];

    //holds one input line per file
    char line[MAX_FILES][MAX_LINE];

    //defines permutation ordering for requested join fields for each file
    size_t join_field_index[MAX_FILES][MAX_FIELDS];

    size_t field_max[MAX_FILES];
    std::fill(field_max, field_max + MAX_FILES, 0);

    char join_format_string[MAX_FILES][MAX_FORMAT_STRING];
    char output_format_string[MAX_FILES][MAX_FORMAT_STRING];

    size_t num_join_fields;

    //join fields for each
    void ** join_fields[MAX_FILES];
    void * join_fields_ptr_buf[MAX_FILES][MAX_JOIN_FIELDS];
    void join_fields_buf[MAX_FILES][MAX_JOIN_FIELDS][MAX_FIELD_BYTES];

    //holds the currently parsed output fields for each file
    void ** output_fields[MAX_FILES];
    void * output_fields_ptr_buf[MAX_FILES][MAX_JOIN_FIELDS];
    void output_fields_buf[MAX_FILES][MAX_JOIN_FIELDS][MAX_FIELD_BYTES];

    void * output_fields_ordered[MAX_FILES * MAX_FIELDS];
    
    std::vector<std::string> fields_tmp;
    
    char file_codes[] = "abcdefghi";

    size_t ofi = 0;
    fields_tmp = split_token(output_field_string, strlen(output_field_string), ",");
    for (size_t f = 0; f != fields_tmp.size(); ++f)
    {
        char file_code;
        size_t field_num;
        sscanf(fields_tmp[f].c_str(), "%zu%c", &field_num, file_code);
        size_t file_num = std::distance(file_codes, index(file_codes, file_code));
        output_fields_ordered[ofi++] = output_fields_ptr_buf[file_num][field_num];
    }

    size_t num_output_fields = ofi;

    //one-time initialization of ordering for output
    //the 

    for (size_t fn = 0; fn != num_input_files; ++fn)
    {
        file[fn] = argv[optind + fn];
        fields_tmp = split_token(fields_string[fn], strlen(fields_string[fn]), ",");
        if (fields_tmp.size() > MAX_JOIN_FIELDS)
        {
            fprintf("Error: number of requested join fields, %Zu"
                    ", is greater than max allowed, %Zu\n", fields_tmp.size(), MAX_JOIN_FIELDS);
            exit(1);
        }

        for (size_t f = 0; f != fields_tmp.size(); ++f)
        {
            field_index[fn][f] = static_cast<size_t>(atoi(fields_tmp[f]));
            field_max[fn] = std::max(field_max[fn], field_index[fn][f]);
        }

        //initialize join_fields and output_fields
        for (size_t f = 0; f != MAX_JOIN_FIELDS; ++f)
        {
            join_fields_ptr_buf[fn][f] = join_fields_buf[fn][f];
            output_fields_ptr_buf[fn][f] = output_fields_buf[fn][f];
        }

        join_fields[fn] = &join_fields_ptr_buf[fn][0];
        output_fields[fn] = &output_fields_ptr_buf[fn][0];


        //construct printf specifiers, i.e. for [2,7,5], construct "%*s\t%*s\t%s\t%*s\t%*s\t%s\t%*s\t%s"
        size_t f = 0;
        for (size_t p = 0; p != field_max[fn]; ++p)
        {
            if (p != 0)
            {
                strcat(format_string[fn], " ");
            }
            if (field_index[fn][f] == p)
            {
                strcat(format_string[fn], "%s");
                ++f;
                if (f == field_index[fn][f])
                {
                    break;
                }
            }
            else
            {
                strcat(format_string[fn], "%*s");
            }
        }
    }

    //check all same join fields.
    size_t num_join_fields = field_max[0] + 1;
    for (size_t fn = 0; fn != num_input_files; ++fn)
    {
        if (num_join_fields != field_max[fn])
        {
            fprintf(stderr, "Error: number of join fields given for each file must be the same"
                    "First file had %Zu fields, file %Zu given %Zu fields\n",
                    num_join_fields, fn, field_max[fn]);
            exit(1);
        }
    }

    //open all FILE handles
    FILE * fh[MAX_FILES];
    for (size_t fn = 0; fn != num_input_files; ++fn)
    {
        fh[fn] = fopen(file[fn], "r");
        if (fh[fn] == NULL)
        {
            fprintf(stderr, "Error: Couldn't open input file %s\n",
                    file[fn]);
            exit(1);
        }
    }

    //read one of each file
    for (size_t fn = 0; fn != num_input_files; ++fn)
    {
        fgets(line[fn], MAX_LINE, fh[fn]);
        LineTools::parse_fields(line[fn], join_format_string[fn], 
                                join_fields[fn]);

        LineTools::parse_fields(line[fn], output_format_string[fn], 
                                output_fields[fn]);

    }

    //parse fields
    //construct the joined line from the join string specification
    //i.e. a join string of 2a7b3c2b
    //from this string we need to know the maximum number of desired reported fields
    //in each file.  The fields are stored in raw (void *) format, and to join them,
    //we go through

    
    
    while (1)
    {
        void ** least_join_field = 
            std::min_element(join_fields, join_fields = num_input_files, less_join_fields);


        // !!! left off here.  Need to retool parse_fields.  See notebook for flowchart.

        for (size_t fn = 0; fn != num_input_files; ++fn)
        {
            
        std::pair<void **, void **> range = 
            std::equal_range(join_fields, join_fields + num_input_files,
                             join_fields[0], less_join_fields);
        
        

    }

    for (size_t fn = 0; fn != num_input_files; ++fn)
    {
        fclose(file[fn]);
    }


    delete file;
    delete line;
    delete line_buf;
    delete field_index;
    delete join_fields_buf;

    return 0;
}
