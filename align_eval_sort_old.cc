//sort an evaluation read by earliest boundary.

/*
  overview:
  1.  load SAM file into memory in contiguous chunks, separated at newlines.
  2.  make a vector of pointers to each line
  3.  sort the vector in memory
  4.  print out each sorted chunk to temp files (or a single file if one chunk)
  5.  merge all sorted chunks, (see below)

  The merging strategy shall be to read smaller chunks from each temp file,
  construct a set for each, and merge them.

  The difficulty is to print out only what is safe, and keep the rest.
  
  1.  Read a fixed chunk into a block of memory for each file (file_blocks[m])
  2.  Construct ordered sets of pointers into each block (line[m])
  3.  Merge the sets (merged)
  4.  Once you get to the 
*/


//Maintain a set of 'next' pointers, one for each block.
//Initialize std::vector<char const*> merged to
//size_t m = 0

#include "align_eval_sort.h"

#include <functional>
#include <vector>
#include <map>
// #include <sys/mman.h>
#include <cstdlib>

#include "align_eval_raw.h"

#include "file_utils.h"
#include "dep/tools.h"
#include "sam_helper.h"

bool less_coord(char const* ref1, size_t pos1,
                char const* ref2, size_t pos2)
{
    int refcmp = strcmp(ref1, ref2);
    return refcmp == -1 ||
        (refcmp == 0 &&
         pos1 < pos2);
}


bool less_samline_min_ag(char * s1, char * s2)
{
    SamLine * sa1 = new SamLine(s1, true);
    SamLine * sa2 = new SamLine(s2, true);
    read_coords ga1;
    read_coords ga2;
    CigarFromSimSAMLine(sa1->qname, sa1->first_read_in_pair(), true, &ga1);
    CigarFromSimSAMLine(sa2->qname, sa2->first_read_in_pair(), true, &ga2);

    char * min_contig1;
    char * min_contig2;
    size_t min_pos1;
    size_t min_pos2;

    if (less_coord(sa1->rname, sa1->zero_based_pos(),
                   ga1.contig, ga1.position))
    {
        min_contig1 = sa1->rname;
        min_pos1 = sa1->zero_based_pos();
    }
    else
    {
        min_contig1 = ga1.contig;
        min_pos1 = ga1.position;
    }
 
    if (less_coord(sa2->rname, sa2->zero_based_pos(),
                   ga2.contig, ga2.position))
    {
        min_contig2 = sa2->rname;
        min_pos2 = sa2->zero_based_pos();
    }
    else
    {
        min_contig2 = ga2.contig;
        min_pos2 = ga2.position;
    }
        
    bool is_less = less_coord(min_contig1, min_pos1, min_contig2, min_pos2);
    delete sa1;
    delete sa2;

    return is_less;
}
    
    
std::vector<char *> find_line_starts(char * lines)
{
    char * lines_tmp = lines;
    std::vector<char *> line_starts;

    line_starts.push_back(lines_tmp);
    while ((lines_tmp = strchr(lines_tmp, '\n')))
    {
        line_starts.push_back(++lines_tmp);
    }
    line_starts.pop_back();
    return line_starts;
}


bool less_first(std::pair<size_t, char *> const& a,
                std::pair<size_t, char *> const& b)
{
    return a.first < b.first;
}



std::vector<std::pair<size_t, char *> > 
build_index(std::vector<char *> const& lines, 
            size_t (* samline_pos)(char const*, CONTIG_OFFSETS const&), 
            CONTIG_OFFSETS const& contig_offsets)
{
    std::vector<char *>::const_iterator lsi;
    std::vector<std::pair<size_t, char *> > lines_index;
    lines_index.reserve(lines.size());
    
    for (lsi = lines.begin(); lsi != lines.end(); ++lsi)
    {
        size_t index = samline_pos(*lsi, contig_offsets);
        lines_index.push_back(std::make_pair(index, *lsi));
    }
    return lines_index;
}


//sort lines in a buffer using a comparison function,
//outputting them to a file
void sort_lines(std::vector<char *> const& lines, 
                size_t (* samline_pos)(char const*, CONTIG_OFFSETS const&), 
                CONTIG_OFFSETS const& contig_offsets,
                FILE * out_fh)
{
    
    
    std::vector<char *>::const_iterator lsi;
    std::vector<std::pair<size_t, char *> > lines_index = 
        build_index(lines, samline_pos, contig_offsets);

    std::sort(lines_index.begin(), lines_index.end(), &less_first);
    std::vector<std::pair<size_t, char *> >::const_iterator lsii;

    for (lsii = lines_index.begin(); lsii != lines_index.end(); ++lsii)
    {
        fprintf(out_fh, "%s\n", (*lsii).second);
        //FileUtils::print_until_delim((*lsii).second, '\n', out_fh);
    }
}





//merge a set of sorted tmp files.  
void merge_files(FILE * sorted_tmp_fhs[], size_t num_tmp_files, 
                 FILE * merged_fh, size_t max_mem)
{

    if (num_tmp_files == 1)
    {
        fseek(sorted_tmp_fhs[0], 0, std::ios::end);
        size_t num_file_bytes = ftell(sorted_tmp_fhs[0]);
        rewind(sorted_tmp_fhs[0]);
        FileUtils::print_chunk_by_line(sorted_tmp_fhs[0], 0, num_file_bytes, merged_fh);
        return;
    }

    size_t chunk_size = max_mem / num_tmp_files;

    std::map<char *, size_t> block_nums;
    std::vector<char **> merged;
    char ** file_blocks = new char *[num_tmp_files];
    char * file_blocks_buffer = new char[num_tmp_files * (chunk_size)];
    for (size_t t = 0; t != num_tmp_files; ++t)
    {
        file_blocks[t] = file_blocks_buffer + (t * (chunk_size));
        block_nums[file_blocks[t]] = t;
    }

    std::vector<std::vector<char *> > line_starts;

    size_t m = 0;

    while(1)
    {
        
        if (merged.empty() || (*(merged[m])++) == NULL)
        {
            //we've run out, or reached the first end in one of the sublists
            char ** current_elem = merged[m];

            char * this_block_cur = *current_elem;
            char * this_block_start = *(current_elem + 2);
            char * this_block_end = *(current_elem + 3);

            merged.clear();
        
            size_t remain = std::distance(this_block_cur, this_block_end);
        
            size_t block_num = block_nums[this_block_start];

            //Refresh and reload the block.

            //1. do a memcpy from ('next', end) to beginning,
            memcpy(this_block_start, this_block_cur, remain);
        
            //2. read file into mem-write pointer to end of block.
            size_t bytes_read = 
                fread(this_block_start + remain, 1, 
                      chunk_size - remain - 1, sorted_tmp_fhs[block_num]);

            size_t buffer_bytes = bytes_read + remain;
            this_block_start[buffer_bytes] = '\0';
        
            //3. Construct a vector of pointers into each block at line starts.
            line_starts[m] = find_line_starts(this_block_start);
            size_t num_line_starts = line_starts[m].size();

            line_starts[m].push_back(NULL);
            line_starts[m].push_back(this_block_start);
            line_starts[m].push_back(this_block_start + buffer_bytes);

            //4. Merge the vectors into one sorted vector of char const**
            for (size_t m = 0; m != num_tmp_files; ++m)
            {
                std::vector<char **> line_start_ptrs(num_line_starts);
                std::vector<char **>::iterator lp;
                std::vector<char *>::iterator l;
                for (lp = line_start_ptrs.begin(), l = line_starts[m].begin(); 
                     lp != line_start_ptrs.end(); ++lp, ++l)
                {
                    *lp = &(*l);
                }

                std::merge(merged.begin(), merged.end(), 
                           line_start_ptrs.begin(), line_start_ptrs.end(),
                           merged.begin());
            }
        }
        else
        {
            FileUtils::print_until_delim(*merged[m++], '\n', merged_fh);
        }
    }
    delete file_blocks;
    delete file_blocks_buffer;

}


int align_eval_sort_usage(char const* sdef, size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval sort [OPTIONS] alignment.sam alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-s  STRING    type of sorting to use {READ_ID, ALIGN, GUIDE, MIN_ALIGN_GUIDE}[%s]\n"
            "-m  INT       number bytes of memory to use [%Zu]\n",
            sdef, mdef);

    fprintf(stderr,
            "Sort orders are:\n"
            "READ_ID: sort by read id\n"
            "ALIGN: sort by alignment position\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_align_eval_sort(int argc, char ** argv)
{

    char * sort_type_def = "MIN_ALIGN_GUIDE";
    char * sort_type = sort_type_def;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    char c;
    while ((c = getopt(argc, argv, "s:m:")) >= 0)
    {
        switch(c)
        {
        case 's': sort_type = optarg; break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        default: return align_eval_sort_usage(sort_type_def, max_mem_def); break;
        }
    }

    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return align_eval_sort_usage(sort_type_def, max_mem_def);
    }

    char const* alignment_sam_file = argv[optind];
    char const* sorted_sam_file = argv[optind + 1];

    FILE * alignment_sam_fh = open_if_present(alignment_sam_file, "r");
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "w");

    size_t const max_line = 10000;

    size_t chunk_size = (max_mem - max_line);

    std::map<std::string, size_t> contig_lengths = ContigLengths(alignment_sam_fh);
    size_t header_length = ftell(alignment_sam_fh);

    CONTIG_OFFSETS contig_offsets = ContigOffsets(contig_lengths);

    std::vector<size_t> chunk_lengths = FileUtils::ChunkLengths(alignment_sam_fh, chunk_size);

    chunk_lengths[0] -= header_length;

    size_t num_chunks = chunk_lengths.size();

    char * chunk_buffer = new char[max_mem];

    size_t (* samline_index)(char const*, CONTIG_OFFSETS const&);
    if (strcmp(sort_type, "MIN_ALIGN_GUIDE") == 0)
    {
        samline_index = &samline_position_min_align_guide;
    }
    else
    {
        fprintf(stderr, "Sorry, %s sort type not implemented\n", sort_type);
        exit(1);
    }


    int file_des;
    char * tmp_file_template = "align_eval_sort.XXXXXX";
    size_t template_length = strlen(tmp_file_template);

    char ** tmp_files = new char *[num_chunks];
    FILE ** tmp_fhs = new FILE *[num_chunks];
    char * tmp_file_buf = new char[num_chunks * (template_length + 1)];
    for (size_t c = 0; c != num_chunks; ++c)
    {
        tmp_files[c] = tmp_file_buf + (c * (template_length + 1));
        strcpy(tmp_files[c], tmp_file_template);
    }

    //print out header.  currently not handled
    FileUtils::cat(chunk_buffer, max_mem, 
                   header_length, alignment_sam_fh, sorted_sam_fh);

    for (size_t chunk = 0; chunk != num_chunks; ++chunk)
    {
        file_des = mkstemp(tmp_files[chunk]);
        tmp_fhs[chunk] = fdopen(file_des, "w");

        //print a chunk, then null terminate it, since this is a mmap
        size_t nbytes_read = 
            fread(chunk_buffer, 1, chunk_lengths[chunk], alignment_sam_fh);

        assert(nbytes_read == chunk_lengths[chunk]);

        chunk_buffer[nbytes_read] = '\0';

        std::vector<char *> samlines = find_line_starts(chunk_buffer);

        //replace all newlines with '\0'
        char * line_end = chunk_buffer;
        while ((line_end = strchr(line_end, '\n')) != NULL)
        {
            *line_end++ = '\0';
        }

        sort_lines(samlines, samline_index, contig_offsets, tmp_fhs[chunk]);
        fflush(tmp_fhs[chunk]);
    }
    
    merge_files(tmp_fhs, num_chunks, sorted_sam_fh, max_mem);

    //remove temporary files
    for (size_t chunk = 0; chunk != num_chunks; ++chunk)
    {
        fclose(tmp_fhs[chunk]);
        remove(tmp_files[chunk]);
    }

    close_if_present(alignment_sam_fh);
    close_if_present(sorted_sam_fh);

    delete chunk_buffer;
    delete tmp_files;
    delete tmp_file_buf;
    delete tmp_fhs;

    return 0;
}
