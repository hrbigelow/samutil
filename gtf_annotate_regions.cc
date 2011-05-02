/*
  Annotate a set of user-defined regions (dna, start, end) with the
  GTF entries they overlap.
 */

#include "nclist.h"
#include "gtf.h"

#include <vector>

int gtf_annotate_regions_usage(char const* rdef, char const* sdef)
{
    fprintf(stderr,
            "\nAuthor: Henry Bigelow,  hbigelow@amgen.com\n\n"
            "gtf_annotate_regions [OPTIONS] annot.gtf regions.txt annot_regions.gtf\n"
            "Augment regions.txt with overlapping annotations in annot.gtf, outputting to annot_regions.gtf.\n\n"
            "Options:\n\n"
            "-e FLAG    (orphans-are-errors) if present, a user-region that does not\n"
            "           overlap any gtf regions is considered an error [false]\n\n"
            "-1 FLAG    (ones-based) if present, regions.txt file assumed to be ones-based\n"
            "           (1,10) denotes the first 10 bases.  if absent, (0,10) denotes them\n\n"
            "-r STRING  (report type) one of {region,annotation,intersection}\n"
            "           For each intersecting region report either the coordinates in region.txt,\n"
            "           or annot.gtf, or the intersection of the two. [%s]\n\n"
            "-s STRING  source field for output in annot_regions.gtf [%s]\n"
            "-t FLAG    (with strand information) if present, expect regions.txt to have four columns:\n"
            "           contig, start_bound, end_bound, strand\n"
            "           If absent, expect regions.txt to have columns (contig, start_bound, end_bound)\n"
            "           In this case, strand will be reported as '.' [false]\n\n"
            ,
            rdef, sdef
            );
    return 1;
}


struct GTFValues
{
    char frame;
    char feature[1000];
    char attribute_string[1000];
    GTFValues(char _f, char const* _fe, char const* _as) : frame(_f)
    {
        strcpy(feature, _fe);
        strcpy(attribute_string, _as);
    }
};

enum ReportType {
    REPORT_REGION,
    REPORT_INTERSECTION,
    REPORT_ANNOTATION
};



int main(int argc, char **argv)
{
    char c;
    bool disallow_orphan_regions = false;
    bool ones_based_regions_file = false;
    bool with_strand = false;

    char const* report_type_def = "intersection";
    char const* report_type_string = report_type_def;

    char const* source_field_def = "gtf_annotate_regions";
    char const* source_field = source_field_def;


    while ((c = getopt(argc, argv, "e1r:s:t")) >= 0)
    {
        //fprintf(stderr, "c = %c, optind = %i\n", c, optind);
        switch(c)
        {
        case 'e': disallow_orphan_regions = true; break;
        case '1': ones_based_regions_file = true; break;
        case 'r': report_type_string = optarg; break;
        case 's': source_field = optarg; break;
        case 't': with_strand = true; break;
        default: return gtf_annotate_regions_usage(report_type_def, source_field_def); break;
        }
    }

    int min_req_args = 3;
    int min_arg_count = optind + min_req_args;

    if (argc != min_arg_count)
    {
        return gtf_annotate_regions_usage(report_type_def, source_field_def);
    }
    
    char * annot_gtf_file = argv[optind];
    char * user_regions_file = argv[optind + 1];
    char * annot_regions_file = argv[optind + 2];

    FILE * annot_gtf_fh = open_or_die(annot_gtf_file, "r", "Annotation file");
    FILE * user_regions_fh = open_or_die(user_regions_file, "r", "User regions file");
    FILE * annot_regions_fh = open_or_die(annot_regions_file, "w", "Output annotated user regions file");

    ReportType report_type;

    if (strcmp(report_type_string, "region") == 0)
    {
        report_type = REPORT_REGION;
    }
    else if (strcmp(report_type_string, "annotation") == 0)
    {
        report_type = REPORT_ANNOTATION;
    }
    else if (strcmp(report_type_string, "intersection") == 0)
    {
        report_type = REPORT_INTERSECTION;
    }
    else
    {
        fprintf(stderr, "Error: optino '-r' must be followed by 'region', 'annotation', or 'intersection'\n");
        exit(1);
    }

    GTFEntry gtf_entry;

    typedef std::vector<GTFValues> GTF_VALUES_VEC;
    typedef std::pair<BY_START::iterator, bool> GTF_INS;

    GTF_VALUES_VEC * gtf_values;
    GTF_INS ins_result;

    IntervalTree * gtf_tree = new IntervalTree();

    while (gtf_entry.get_next_record(annot_gtf_fh))
    {
        if (! gtf_entry.is_data_line)
        {
            continue;
        }

        gtf_values = new GTF_VALUES_VEC();
        
        ins_result = 
            gtf_tree->insert(gtf_entry.seqname, gtf_entry.start - 1, gtf_entry.end, gtf_values);
        if (! ins_result.second)
        {
            delete gtf_values;
            gtf_values = static_cast<GTF_VALUES_VEC *>((*ins_result.first).second->payload);
        }
        (*gtf_values).push_back(GTFValues(gtf_entry.frame, gtf_entry.feature, gtf_entry.attribute_string));
    }
    //print_tree(gtf_tree, 0);
    assert(check_tree(gtf_tree));

    int ones_adjust = ones_based_regions_file ? 1 : 0;

    std::vector<IntervalTree const*> overlapping_annots;
    std::vector<IntervalTree const*>::iterator ovit;

    int nscanned;
    char reg_contig[1024];
    size_t reg_start_bound, reg_end_bound, reported_start_bound, reported_end_bound;

    char reg_strand = '.';

    char score[] = ".";

    //check if empty file
    char next_char;
    if ((next_char = fgetc(user_regions_fh)) == EOF)
    {
        //do nothing
    }
    else
    {
        ungetc(next_char, user_regions_fh);
    }

    while (! feof(user_regions_fh))
    {
        if (with_strand)
        {
            nscanned = fscanf(user_regions_fh, "%s\t%zu\t%zu\t%c\n",
                              reg_contig, &reg_start_bound, &reg_end_bound, &reg_strand);
            if (nscanned != 4)
            {
                fprintf(stderr, 
                        "Error: bad format in user regions file %s\n"
                        "Must be: 'contig\\tstart\\tend\\tstrand\\n'", 
                        user_regions_file);
                exit(1);
            }
        }
        else
        {
            nscanned = fscanf(user_regions_fh, "%s\t%zu\t%zu\n",
                              reg_contig, &reg_start_bound, &reg_end_bound);
            if (nscanned != 3)
            {
                fprintf(stderr, 
                        "Error: bad format in user regions file %s\n"
                        "Must be: 'contig\\tstart\\tend\\n'", 
                        user_regions_file);
                exit(1);
            }
        }
        reg_start_bound -= ones_adjust;

        //the set of overlapping annotated regions
        overlapping_annots = gtf_tree->overlap(reg_contig, reg_start_bound, reg_end_bound);

        if (disallow_orphan_regions && overlapping_annots.empty())
        {
            fprintf(stderr, "Error: region %s:%Zu\t%Zu does not overlap with any annotation "
                    "and this is not allowed by flag -e\n", 
                    reg_contig, reg_start_bound, reg_end_bound);
            exit(1);
        }

        for (ovit = overlapping_annots.begin(); ovit != overlapping_annots.end(); ++ovit)
        {
            IntervalTree const& node = *(*ovit);
            GTF_VALUES_VEC const& gtf_values = *(static_cast<GTF_VALUES_VEC *>(node.payload));

            assert(gtf_values.size() == 1);
            switch(report_type)
            {
            case REPORT_REGION:
                reported_start_bound = reg_start_bound;
                reported_end_bound = reg_end_bound;
                break;
            case REPORT_INTERSECTION:
                reported_start_bound = std::max(reg_start_bound, node.start);
                reported_end_bound = std::min(reg_end_bound, node.end);
                break;
            case REPORT_ANNOTATION:
                reported_start_bound = node.start;
                reported_end_bound = node.end;
                break;
            default:
                fprintf(stderr, "Unknown region report type");
                exit(1);
                break;
            }

            for (GTF_VALUES_VEC::const_iterator git = gtf_values.begin(); 
                 git != gtf_values.end(); ++git)
            {
                fprintf(annot_regions_fh, "%s\t%s\t%s\t%Zu\t%Zu\t%s\t%c\t%c\t%s\n",
                        reg_contig, source_field, (*git).feature, 
                        reported_start_bound + 1, 
                        reported_end_bound,
                        score, reg_strand, (*git).frame, (*git).attribute_string);
            }
        }
    }

    free_payload<GTF_VALUES_VEC>(gtf_tree);
    
    fclose(annot_gtf_fh);
    fclose(user_regions_fh);
    fclose(annot_regions_fh);
    
}
