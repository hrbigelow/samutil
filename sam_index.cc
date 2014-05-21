#include "sam_index.h"
#include "seq_projection.h"
#include "file_utils.h"
#include "gtf.h"

//#include "sam_helper.h"
#include "sam_flag.h"

#include <cassert>
#include <map>
#include <set>
#include <endian.h>
#include <byteswap.h>

/*
  functions for indexing the lines in a SAM or BAM file.
 */

sam_index::sam_index(size_t *_key, size_t _so, int _ll) : 
    start_offset(_so), line_length(_ll) { 
    key.raw[0] = _key[0];
    key.raw[1] = _key[1];
}


sam_index::sam_index() :
    start_offset(0), line_length(0) { 
    key.raw[0] = 0;
    key.raw[1] = 0;
}

//partial ordering
bool samidx_less_offset(sam_index const& a, sam_index const& b)
{
    return a.start_offset < b.start_offset;
}


//partial ordering
bool samidx_less_key(sam_index const& a, sam_index const& b)
{
    return 
        a.key.raw[0] < b.key.raw[0] ||
        (a.key.raw[0] == b.key.raw[0] &&
         (a.key.raw[1] < b.key.raw[1] ||
          (a.key.raw[1] == b.key.raw[1] &&
           (a.start_offset < b.start_offset))));
}

bool samidx_equal_key(sam_index const& a, sam_index const& b)
{
    return 
        a.key.raw[0] == b.key.raw[0] &&
        a.key.raw[1] == b.key.raw[1];
}




contig_dict::contig_dict() : length(NULL), offset(NULL), space(NULL), projections(NULL) { }

contig_dict::~contig_dict() {
    if (length != NULL)
    {
        delete length;
    }
    if (offset != NULL)
    {
        delete offset;
    }
    if (space != NULL)
    {
        delete space;
    }
    if (projections != NULL)
    {
        delete [] projections;
    }
    for (HMAP::iterator hit = name_map.begin(); hit != name_map.end(); ++hit)
    {
        delete (*hit).first;
    }
    for (HMAP::iterator hit = projected_name_map.begin(); hit != projected_name_map.end(); ++hit)
    {
        delete (*hit).first;
    }
}


SAM_INDEX_TYPE sam_index_type_from_string(char const* index_type)
{
    static char const* index_type_strings[] = {
        "FRAGMENT",
        "FRAGMENT_ALIGN",
        "FRAGMENT_MAGALN",
        "FRAGMENT_PROJALN",
        "ALIGN",
        "MIN_ALIGN_GUIDE",
        "PROJALIGN"
    };


    SAM_INDEX_TYPE found_type = SAM_INDEX_UNDEFINED;

    for (size_t i = 0; i != 5; ++i)
    {
        if (strcmp(index_type, index_type_strings[i]) == 0)
        {
            found_type = SAM_INDEX_TYPE(i);
            break;
        }
    }
    return found_type;
}


// 
void init_contig_length_offset (char const* sam_header, contig_dict * dict)
{
    char * hdr_copy = new char[strlen(sam_header) + 1];
    strcpy(hdr_copy, sam_header);

    char * last_fragment;
    std::vector<char *> header_lines =
        FileUtils::find_complete_lines_nullify(hdr_copy, & last_fragment);
    
    assert(strlen(last_fragment) == 0);

    std::vector<char *>::const_iterator iter;
    char const* line;

    // pre-scan header_lines to count number of SQ lines
    size_t sq_count = 0;
    for (iter = header_lines.begin(); iter != header_lines.end(); ++iter)
    {
        line = *iter;
        if (strncmp(line, "@SQ", 3) == 0)
        {
            sq_count++;
        }
    }
    if (sq_count == 0)
    {
        fprintf(stderr, "No SQ lines in header, or no SN / LN fields defining contigs\n");
        exit(1);
    }

    dict->length = new size_t[sq_count + 1];
    dict->offset = new size_t[sq_count + 1];

    size_t ord = 0;
    size_t contig_offset = 0;
    for (iter = header_lines.begin(); iter != header_lines.end(); ++iter)
    {
        //if there is a SQ tag, parse and record
        line = *iter;
        if (strncmp(line, "@SQ", 3) == 0)
        {
            //this is a SQ line
            char contig_name[1000];
            
            size_t contig_length;
            // need to convert this to something tag-agnostic
            sscanf(line, "@SQ\tSN:%[^\t]\tLN:%zu", contig_name, &contig_length);
            char * contig_name_copy = new char[strlen(contig_name) + 1];
            strcpy(contig_name_copy, contig_name);
            dict->name_map.insert(std::make_pair(contig_name_copy, ord));
            dict->length[ord] = contig_length;
            dict->offset[ord] = contig_offset;

            contig_offset += contig_length;
            ord++;
        }
    }
    // now add one for '*' virtual contig (the 'unmapped' symbolic contig)
    dict->length[ord] = 0;
    dict->offset[ord] = contig_offset;
    char * contig_name_copy = new char[2];
    strcpy(contig_name_copy, "*");
    dict->name_map.insert(std::make_pair(contig_name_copy, ord));

    delete hdr_copy;

    // now rehash the name map
    dict->name_map.rehash(dict->name_map.size() * 10);

    std::map<size_t, size_t> bucket_hist;
    for (size_t i = 0; i != dict->name_map.bucket_count(); ++i)
    {
        bucket_hist[dict->name_map.bucket_size(i)]++;
    }
    // for (std::map<size_t, size_t>::const_iterator it = bucket_hist.begin();
    //      it != bucket_hist.end(); ++it)
    // {
    //     fprintf(stderr, "size %zu buckets: %zu\n", (*it).first, (*it).second);
    // }

}


void init_sequence_projection(char const* gtf_file, contig_dict * dict)
{
    char const* species = "dummy";
    std::set<SequenceProjection> tx_to_genome = 
        gtf_to_sequence_projection(gtf_file, species);
    
    std::set<SequenceProjection>::const_iterator tgi;
    
    assert(dict->projections == NULL);
    dict->projections = new SequenceProjection[tx_to_genome.size()];

    size_t ord = 0;
    for (tgi = tx_to_genome.begin(); tgi != tx_to_genome.end(); ++tgi)
    {
        char * transcript = new char[(*tgi).source_dna.size() + 1];
        strcpy(transcript, (*tgi).source_dna.c_str());
        dict->projections[ord] = SequenceProjection(*tgi);
        dict->projected_name_map.insert(std::make_pair(transcript, ord));
        ord++;
    }
}





size_t samidx_flattened_position(char const* contig,
                                 size_t position, 
                                 contig_dict const* dict)
{
    
    HMAP::const_iterator contig_iter = dict->name_map.find(contig);
        
    if (contig_iter == dict->name_map.end())
    {
        fprintf(stderr, "flattened_position: error: rname %s (at %Zu) "
                "does not exist in provided contig index\n",
                contig, position);
        exit(1);
    }
    
    return dict->offset[(*contig_iter).second] + position;
}


//produce index consistent with ordering [rname, pos]
size_t samline_position_align(char const* samline, contig_dict const* dict)
{

    char contig[1024];
    size_t ones_based_pos;
    sscanf(samline, "%*[^\t]\t%*u\t%[^\t]\t%zu", contig, &ones_based_pos);

    size_t zero_based_pos = ones_based_pos == 0 ? 0 : ones_based_pos - 1;

    return samidx_flattened_position(contig, zero_based_pos, dict);
}


// For all records, ordering is [rname, pos]. Unmapped reads come at
// the end.

// Before, for unmapped reads, I had fragment_id field sorted as well,
// but am taking it out, unless there is some hidden necessity I discover
// later...
   
// By convention, the pseudo-contig '*' represents 'unmapped', and is
// sorted at the end
size_t samline_projected_position_align(char const* samline, contig_dict const* dict)
{

    char * contig;
    char * cigar;

    size_t ones_based_pos;
    size_t flat_pos;

    sscanf(samline, "%*[^\t]\t%*u\t%a[^\t]\t%zu\t%*i\t%as", 
           &contig, &ones_based_pos, &cigar);

    size_t zero_based_pos = ones_based_pos == 0 ? 0 : ones_based_pos - 1;

    HMAP::const_iterator proj_iter = dict->projected_name_map.find(contig);
    if (proj_iter == dict->projected_name_map.end())
    {
        // don't do any projection
        flat_pos = samidx_flattened_position(contig, zero_based_pos, dict);
    }
    else
    {
        SequenceProjection const& proj = dict->projections[(*proj_iter).second];
        flat_pos = samidx_flattened_position(proj.target_dna.c_str(), 
                                          ExpandedStartPos(proj, zero_based_pos, cigar),
                                          dict);
    }

    // this doesn't seem that it should be necessary...
    // if (strcmp(contig, "*") == 0)
    // {
    //     flat_pos += (this->*(this->parse_fragment_id))(samline);
    // }
    
    delete contig;
    delete cigar;
    return flat_pos;
}


//compute the minimum collapsed start position between the alignment
//and the guide.
size_t samline_position_min_align_guide(char const* samline, contig_dict const* dict)
{
    char guide_chrom_left[32];
    char guide_chrom_right[32];
    int read_num_left, read_num_right;
    SamFlag flag;
    size_t guide_pos_left, guide_pos_right;
    char align_chrom[32];
    size_t align_pos;

    size_t flag_raw;

    int nfields = 
        sscanf(samline, 
               "%*[^:]:" // id field.  ignored
               "read%i:%[^:]:%*c:%zu:%*[^:]:%*u:"  //read chunk
               "read%i:%[^:]:%*c:%zu:%*[^:]:%*u:" //second read chunk
               "fragment_size:%*i\t"
               "%zu\t%[^\t]\t%zu", //part of alignment
               &read_num_left, guide_chrom_left, &guide_pos_left, 
               &read_num_right, guide_chrom_right, &guide_pos_right,
               &flag_raw, align_chrom, &align_pos);
    
    if (nfields != 9)
    {
        fprintf(stderr, "SamOrder::samline_position_min_align_guide:\n"
                "Bad ID format for this sort type:\n%s\n", samline);
        exit(1);
    }

    flag.set_raw(flag_raw);

    int read_num = flag.first_fragment_in_template ? 1 : 2;

    char * guide_chrom = read_num == read_num_left ? guide_chrom_left : guide_chrom_right;
    size_t guide_pos = read_num == read_num_left ? guide_pos_left : guide_pos_right;

    HMAP::const_iterator guide_iter = dict->name_map.find(guide_chrom);

    if (guide_iter == dict->name_map.end())
    {
        fprintf(stderr, "This samline has its guide contig (%s) "
                "that is not in SAM contig index\n\n%s\n", 
                guide_chrom, samline);
        exit(1);
    }
    size_t guide_index = dict->offset[(*guide_iter).second] + guide_pos;

    size_t align_index = 
        samidx_flattened_position(align_chrom, align_pos, dict);

    return std::min(guide_index, align_index);

}


/* ************************************************************ */
// Fragment ID parsers
//compute the flattened coordinate start position.  
//read id that starts out with illumina format:
//@LANE:TILE:XPOS:YPOS
// The layout, with high bits to the left, is the same
// llllllll tttttttt ttttxxxx xxxxxxxx xxxxxxxx xxyyyyyy yyyyyyyy yyyyyyyy
// 256 lanes, 4096 tiles, 4194304 x and y positions
/*
struct illumina_id
{
    unsigned int ypos : 22;
    unsigned int xpos : 22;
    unsigned int tile : 12;
    unsigned int lane : 8;
    size_t get_raw() const
    {
        return
            static_cast<size_t>
            (this->ypos % (1L<<22)
             | (this->xpos % (1L<<22))<<22
             | (this->tile % (1L<<12))<<44
             | (this->lane % (1L<<8))<<56);
    }
};


// convert all 
inline void convert16_aux(char const* raw, size_t * converted)
{
    strncpy(reinterpret_cast<char *>(converted), raw, 16);
    switch (BYTE_ORDER)
    {
    case LITTLE_ENDIAN:
        converted[0] = bswap_64(converted[0]);
        converted[1] = bswap_64(converted[1]);
        break;
    case BIG_ENDIAN:
        // do nothing
        break;
    case PDP_ENDIAN:
        fprintf(stderr, "Error: PDP_ENDIAN machines not supported\n");
        exit(56);
        break;
    }

}
*/


bool parse_fragment_id_numeric(char const* qname, bool frag_first, size_t position, sam_index *idx)
{
    size_t id;
    int nfields_read = sscanf(qname, "%zu", &id);
    if (nfields_read != 1)
    {
        return false;
    }
    if (frag_first)
    {
        idx->key.raw[0] = id;
        idx->key.raw[1] = position;
    }
    else
    {
        idx->key.raw[0] = position;
        idx->key.raw[1] = id;
    }
    return true;
}





// Illumina format
bool parse_fragment_id_illumina(const char *qname,
                                bool frag_first,
                                size_t position,
                                sam_index *idx,
                                index_dict_t *dict)
{

    char flowcell[256];
    int end_pos;

    unsigned int flowcell_id, lane, tile, xpos, ypos;

    int nfields_read = sscanf(qname, "%[^:]:%u:%u:%u:%u%n", 
                              flowcell, &lane, &tile, &xpos, &ypos,
                              &end_pos);

    // if there is extra stuff on qname after ypos,
    // it must be one of (<spaces>, '#', '/', or '\0')
    // anything else generates an error
    if (nfields_read != 5 
        || (
            ! isspace(qname[end_pos])
            && qname[end_pos] != '#' 
            && qname[end_pos] != '/'
            && qname[end_pos] != '\0')
        )
    {
        return false;
    }
    index_dict_t::iterator fit = (*dict).find(flowcell);
    if (fit == (*dict).end())
    {
        char * flowcell_copy = new char[strlen(flowcell) + 1];
        strcpy(flowcell_copy, flowcell);
        fit = (*dict).insert((*dict).end(), std::make_pair(flowcell_copy, (*dict).size()));
    }
    // now, fit should be valid
    flowcell_id = (*fit).second;

    if (frag_first)
    {
        idx->key.fp.F = flowcell_id;
        idx->key.fp.L = lane;
        idx->key.fp.T = tile;
        idx->key.fp.X = xpos;
        idx->key.fp.Y = ypos;
        idx->key.fp.P = position;
    }
    else
    {
        idx->key.pf.F = flowcell_id;
        idx->key.pf.L = lane;
        idx->key.pf.T = tile;
        idx->key.pf.X = xpos;
        idx->key.pf.Y = ypos;
        idx->key.pf.P = position;
    }

    return true;
}


/*
  Parses Casava 1.8 read id format.  This routine ignores
  instrument-name, run ID, flowcell ID, and is therefore unsuitable for
  running data that comes from different flowcells etc.

  I will remedy this in the future.

  @<instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos> \
  <read number>:<is filtered>:<control number>:<barcode sequence>

  The spec is vague w.r.t header format.  But, I will assume that the statement:

  The first line is prefixed by the “@” symbol and contains the read
  name. These names are parsed until the first encountered
  whitespace. Due to this behavior, adding additional tags to the header
  line is not problematic for extant FASTQ parsers.

  means that a valid first fastq line may or may not contain a space
  followed by extra characters.

  So, another valid format could be:


  @<instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos>

*/
bool parse_fragment_id_casava_1_8(const char *qname,
                                  bool frag_first,
                                  size_t position,
                                  sam_index *idx,
                                  index_dict_t *dict)
{
    int end_pos;
    int hash_end_pos;

    unsigned int flowcell_id, lane, tile, xpos, ypos;

    int nfields_read = sscanf(qname, "%*[^:]:%*[^:]:%*[^:]%n:%u:%u:%u:%u%n", 
                              &hash_end_pos, &lane, &tile, &xpos, &ypos, &end_pos);
    if (nfields_read != 4 ||
        (! isspace(qname[end_pos]) && qname[end_pos] != '\0'))
    {
        return false;
    }

    char flowcell[256];
    strncpy(flowcell, qname, hash_end_pos);
    flowcell[hash_end_pos] = '\0';

    index_dict_t::iterator fit = (*dict).find(flowcell);
    if (fit == (*dict).end())
    {
        char * flowcell_copy = new char[strlen(flowcell) + 1];
        strcpy(flowcell_copy, flowcell);
        fit = (*dict).insert((*dict).end(), std::make_pair(flowcell_copy, (*dict).size()));
    }

    // now, fit should be valid
    flowcell_id = (*fit).second;

    if (frag_first)
    {
        idx->key.fp.F = flowcell_id;
        idx->key.fp.L = lane;
        idx->key.fp.T = tile;
        idx->key.fp.X = xpos;
        idx->key.fp.Y = ypos;
        idx->key.fp.P = position;
    }
    else
    {
        idx->key.pf.F = flowcell_id;
        idx->key.pf.L = lane;
        idx->key.pf.T = tile;
        idx->key.pf.X = xpos;
        idx->key.pf.Y = ypos;
        idx->key.pf.P = position;
    }

    return true;
}


bool parse_fragment_id_zero(char const* /* qname */, sam_index *idx)
{
    idx->key.raw[0] = 0;
    idx->key.raw[1] = 0;

    return true;
}




SAM_QNAME_FORMAT qname_format(char const* sam_dataline)
{
    SAM_QNAME_FORMAT qname_format;
    sam_index dummy;
    size_t p = 0;
    index_dict_t dict;

    if (strlen(sam_dataline) == 0)
    {
        qname_format = SAM_FORMAT_UNDEFINED;
    }
    else if (parse_fragment_id_numeric(sam_dataline, false, p, &dummy))
    {
        qname_format = SAM_NUMERIC;
    }
    else if (parse_fragment_id_illumina(sam_dataline, false, p, &dummy, &dict))
    {
        qname_format = SAM_ILLUMINA;
    }
    else if (parse_fragment_id_casava_1_8(sam_dataline, false, p, &dummy, &dict))
    {
        qname_format = SAM_CASAVA18;
    }
    else
    {
        fprintf(stderr, "Error: SamOrder: Don't have a parser for this qname format: %s\n",
                sam_dataline);
        exit(1);
    }

    return qname_format;
}


// set the index fields of idx for this sam_line according to the chosen type and format
void set_sam_index(char const* samline, SAM_INDEX_TYPE itype, SAM_QNAME_FORMAT qfmt, 
                   contig_dict const* dict, 
                   index_dict_t *flowcell_dict,
                   sam_index * idx)
{
    size_t pos = 0;

    // calculate position
    switch(itype)
    {
    case SAM_INDEX_FID_ALN:
    case SAM_INDEX_ALN_FID:
        pos = samline_position_align(samline, dict);
        break;
    case SAM_INDEX_FID_MAGALN:
    case SAM_INDEX_MAGALN_FID:
        pos = samline_position_min_align_guide(samline, dict);
        break;
    case SAM_INDEX_FID_PROJALN:
    case SAM_INDEX_PROJALN_FID:
        pos = samline_projected_position_align(samline, dict);
        break;
    case SAM_INDEX_FID:
    case SAM_INDEX_UNDEFINED:
        break;
    }

    bool frag_first = false;

    // initialize position field
    switch(itype)
    {
    case SAM_INDEX_FID:
    case SAM_INDEX_FID_MAGALN:
    case SAM_INDEX_FID_PROJALN:
    case SAM_INDEX_FID_ALN:
        frag_first = true;
        break;

    case SAM_INDEX_ALN_FID:
    case SAM_INDEX_MAGALN_FID:
    case SAM_INDEX_PROJALN_FID:
        frag_first = false;
        break;
    case SAM_INDEX_UNDEFINED:
        break;
    }
    
    // initialize all other fields
    switch(qfmt)
    {
    case SAM_NUMERIC: parse_fragment_id_numeric(samline, frag_first, pos, idx); break;
    case SAM_ILLUMINA: parse_fragment_id_illumina(samline, frag_first, pos, idx, flowcell_dict); break;
    case SAM_CASAVA18: parse_fragment_id_casava_1_8(samline, frag_first, pos, idx, flowcell_dict); break;
    case SAM_FORMAT_UNDEFINED: 
    default:
        fprintf(stderr, "Error: Unknown QNAME format\n");
        exit(23);
        break;
    }
}


// update this index's flowcell id in place
void sam_update_flowcell_id(const unsigned int * remap, SAM_INDEX_TYPE itype, sam_index * idx)
{
    switch(itype)
    {
    case SAM_INDEX_FID:
    case SAM_INDEX_FID_ALN:
    case SAM_INDEX_FID_MAGALN:
    case SAM_INDEX_FID_PROJALN:
        idx->key.fp.F = remap[idx->key.fp.F];
        break;
    case SAM_INDEX_ALN_FID:
    case SAM_INDEX_MAGALN_FID:
    case SAM_INDEX_PROJALN_FID:
        idx->key.pf.F = remap[idx->key.pf.F];
        break;
    case SAM_INDEX_UNDEFINED:
        break;
    }

}


sam_index samidx_make_max()
{
    size_t sentinel_ary[2] = {SIZE_MAX, SIZE_MAX};
    return sam_index(sentinel_ary, SIZE_MAX, 0);
}



bool eqstr::operator()(const char* s1, const char* s2) const
{
    return strcmp(s1, s2) == 0;
}


//a hashing function.  Expects that k will either start or end with a string of digits.
//if it does, the integer value of those digits is used.
//Otherwise, if it is '*', zero is used.,
//Otherwise, the default hash of its string value is used.
size_t to_integer::operator()(char const* k) const
{
    size_t val = 0;
    char dummy; //used as a test whether sscanf reached the end of the string
    size_t dummy_int;

    if (sscanf(k, "%*[^0-9]%zu%c", &val, &dummy) == 1 
        || sscanf(k, "%zu%*[^0-9]%zu", &val, &dummy_int) == 1)
    {
        return val;
    }
    else if (k[0] == '*')
    {
        return 0;
    }
    else
    {
        return static_cast<std::hash<std::string> >(*this)(std::string(k, strlen(k)));
    }
}
