#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "file_utils.h"

#include <cstring>
#include <cassert>
#include <cstdarg>



//optionally prints to a string or a file.  to make a string behave
//the same as a file, though, we need to increment it after printing.
int FileUtils::switch_printf(bool is_file, void ** buf, char const* format, ...)
{
    va_list arg;
    int done;
    
    va_start (arg, format);
    
    if (is_file)
    {
        done = vfprintf(static_cast<FILE *>(*buf), format, arg);
    }
    else
    {
        char * stringbuf = static_cast<char *>(*buf);
        done = vsprintf(stringbuf, format, arg);
        (*buf) += done;
    }

    va_end (arg);

    return done;
}



//safely reads and writes 'bytes_to_print' bytes from input_fh to
//output_fh or stops when input_fh reaches the end.
//internally uses buf (of size bufsize)
void FileUtils::cat(char * buf, 
                    size_t bufsize, 
                    size_t bytes_to_print, 
                    FILE * input_fh, 
                    FILE * output_fh)
{
    size_t bytes_left = bytes_to_print;

    while (! feof(input_fh))
    {
        size_t bytes_chunk = std::min(bytes_left, bufsize);
        fread(buf, 1, bytes_chunk, input_fh);
        fwrite(buf, 1, bytes_chunk, output_fh);
        bytes_left -= bytes_chunk;
        if (bytes_left == 0)
        {
            break;
        }
    }
    fflush(output_fh);
}


//virtually splits the contents  of input_fh into chunks, separated at
//newlines,  of  approximately  'chunk_approx_size', and  outputs  the
//'chunk'th chunk to output_fh.  

//the chunk size must be approximate because there is no guarantee
//that the newline will fall exactly on a chunk boundary.
void FileUtils::print_chunk_by_line(FILE * input_fh, size_t chunk, 
                                    size_t chunk_approx_size,
                                    FILE * output_fh)
{

    size_t const io_blocksize = 32 * 1024;
    char * buf = new char[io_blocksize];

    //go to beginning
    rewind(input_fh);
    
    size_t cur_chunk = 0;
    while (! feof(input_fh))
    {
        if (cur_chunk == chunk)
        {
            cat(buf, io_blocksize, chunk_approx_size, input_fh, output_fh);

            //now get remainder of line
            buf[0] = '\0';
            fgets(buf, io_blocksize - 1, input_fh);
            fprintf(output_fh, "%s", buf);
            break;
        }
        else
        {
            fseek(input_fh, chunk_approx_size, std::ios::cur);
            fgets(buf, io_blocksize - 1, input_fh);
        }
        ++cur_chunk;
    }
    delete [] buf;
}


//return the list of chunk lengths, where the position is the
//closest newline after chunk_approx_size
std::vector<size_t> FileUtils::ChunkLengths(gzFile input_fh, size_t chunk_approx_size,
                                            size_t max_line)
{
    std::vector<size_t> lengths;
    char * dummy_buf = new char[max_line]; //allow for an extremely long line
    gzrewind(input_fh);
    z_off_t pos, prev_pos = 0;

    while (! gzeof(input_fh))
    {
        gzseek(input_fh, chunk_approx_size - max_line, SEEK_CUR);
        gzgets(input_fh, dummy_buf, max_line - 1);
        pos = gztell(input_fh);
        if (gzeof(input_fh))
        {
            gzclearerr(input_fh);
            gzseek(input_fh, prev_pos, SEEK_SET);
            unsigned int len = UINT32_MAX - 1;
            char * tmp_buf = new char[len + 1];
            while (! gzeof(input_fh))
            {
                gzread(input_fh, tmp_buf, len);
            }
            delete [] tmp_buf;
            pos = gzseek(input_fh, 0L, SEEK_CUR);
        }
        lengths.push_back(static_cast<size_t>(pos - prev_pos));
        prev_pos = pos;
    }
    delete [] dummy_buf;
    gzrewind(input_fh);

    return lengths;
}


//
std::vector<char *> 
find_complete_lines_aux(char * lines, bool nullify_newlines, char ** last_fragment)
{
    char * lines_tmp = lines;
    std::vector<char *> line_starts;
    size_t total_length = strlen(lines);
    if (total_length == 0)
    {
        (*last_fragment) = lines;
        return line_starts;
    }

    size_t first_line_length = 
        std::max(100UL, 
                 static_cast<size_t>(std::distance(lines_tmp, strchrnul(lines_tmp, '\n'))));

    size_t est_num_lines = total_length / first_line_length;

    line_starts.reserve(est_num_lines);
        
    line_starts.push_back(lines_tmp);

    if (nullify_newlines)
    {
        while ((lines_tmp = strchr(lines_tmp, '\n')))
        {
            *lines_tmp = '\0';
            line_starts.push_back(++lines_tmp);
        }
    }
    else
    {
        while ((lines_tmp = strchr(lines_tmp, '\n')))
        {
            line_starts.push_back(++lines_tmp);
        }
    }
    (*last_fragment) = line_starts.back();
    line_starts.pop_back();

    return line_starts;
}

std::vector<char *> 
FileUtils::find_complete_lines(char * lines, char ** last_fragment)
{
    return find_complete_lines_aux(lines, false, last_fragment);
}

std::vector<char *> 
FileUtils::find_complete_lines_nullify(char * lines, char ** last_fragment)
{
    return find_complete_lines_aux(lines, true, last_fragment);
}




void FileUtils::print_until_delim(char const* data, char delim, FILE * out_fh)
{
    size_t bytes = strchr(data, delim) - data;
    fwrite(data, 1, bytes, out_fh);
}



BufferedFile::BufferedFile(char const* _file, size_t _chunk_size, size_t _max_line)
    : chunk_size(_chunk_size), max_line(_max_line)
{
    this->file = new char[strlen(_file) + 1];
    strcpy(this->file, _file);
    this->chunk_buffer = new char[chunk_size + _max_line + 1];

    this->line_iter = this->line_starts.end();
}


BufferedFile::~BufferedFile()
{
    delete [] this->file;
    delete [] this->chunk_buffer;
    if (this->fh != NULL)
    {
        gzclose(this->fh);
    }
}


bool BufferedFile::initialize()
{

    if (strcmp(this->file, "/dev/null") == 0 ||
        strcmp(this->file, "") == 0)
    {
        this->valid = false;
        return this->valid;
    }

    this->fh = gzopen(this->file, "r");
    if (this->fh == NULL)
    {
        fprintf(stderr, "BufferedFile::initialize(): "
                "could not open file %s for reading\n", this->file);
        return false;
    }

    this->chunk_lengths = 
        FileUtils::ChunkLengths(this->fh, this->chunk_size, this->max_line);
    this->chunk_iter = this->chunk_lengths.begin();
    gzrewind(this->fh);

    this->valid = true;
    return this->valid;
}



//return a null-terminated pointer to the next num_lines lines in the
//buffer, advancing the chunk or rewinding/reloading the file if there
//are fewer than that remaining.
char * BufferedFile::next_n_lines(size_t num_lines, bool * advanced_chunk, bool * reloaded_file)
{
    //is our current line_iter valid?  if not, we need to read another chunk
    size_t lines_left = std::distance(this->line_iter, this->line_starts.end());
    size_t nbytes_unused;
    char * last_fragment;

    *advanced_chunk = lines_left < num_lines;

    if (*advanced_chunk)
    {
        //reading another chunk.  is our chunk_iter valid?
        *reloaded_file = (this->chunk_iter == this->chunk_lengths.end());

        if (*reloaded_file)
        {
            //if not, we need to reload the file
            //assume that if the user is requesting num_lines at a time,
            //the file contains a multiple of num_lines
            if (lines_left != 0)
            {
                fprintf(stderr, "Error: BufferedFile::next_n_lines: %zu lines requested at a time, \n"
                        "but there are %zu lines left in this file\n", num_lines, lines_left);
                exit(1);
            }
            gzrewind(this->fh);
            this->chunk_iter = this->chunk_lengths.begin();
        }

        //now, our chunk_iter is valid.  reload the buffer
        //carry over the last 'lines_left' into the next buffer
        //these will be terminated by newlines
        size_t carry_over = 0;
        if (lines_left != 0)
        {
            carry_over = strlen(*this->line_iter);
            strcpy(this->chunk_buffer, *this->line_iter);
        }

        size_t num_bytes_read = 
            gzread(this->fh, this->chunk_buffer + carry_over, *this->chunk_iter);

        assert(num_bytes_read == *this->chunk_iter);
        ++this->chunk_iter;

        this->chunk_buffer[num_bytes_read + carry_over] = '\0';

        this->line_starts = FileUtils::find_complete_lines(this->chunk_buffer, & last_fragment);

        nbytes_unused = strlen(last_fragment);

        this->line_iter = this->line_starts.begin();

        //replace newlines with nulls
        char * line_end = chunk_buffer;
        size_t nth = 0;
        while ((line_end = strchr(line_end, '\n')) != NULL)
        {
            ++nth;
            if (nth % num_lines == 0)
            {
                *line_end = '\0';
            }
            ++line_end;
        }

    }
    else
    {
        *reloaded_file = false;
    }

    char * lines = *line_iter;
    line_iter += num_lines;
    return lines;
    
}
