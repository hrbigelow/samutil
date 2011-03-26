#include "litestream.h"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <cassert>
#include <cstring>

litestream::litestream(std::istream & base, 
                       litestream::postype size) :
    stream(&base), 
    buffer_start(NULL), 
    buffer_current(NULL), 
    buffer_end(NULL),
    buffer_offset(0),
    stream_size(0),
    request_size(size) { }

litestream::litestream(std::istream & base) :
    stream(&base), 
    buffer_start(NULL), 
    buffer_current(NULL), 
    buffer_end(NULL),
    buffer_offset(0),
    stream_size(0),
    request_size(-1) { }

void litestream::open(const char * filename){
    static_cast<std::ifstream *>(stream)->open(filename, 
                                               std::ios::in | std::ios::binary);
  
    stream->seekg(0, std::ios::end);
    stream_size = stream->tellg();
    stream->seekg(0, std::ios::beg);
  
    buffer_size = request_size == -1 ? stream_size : request_size;
    while (1){
        buffer_start = new char[buffer_size + 1];
        if (buffer_start != NULL) break;
        buffer_size /= 2;
    }
    buffer_offset = 0;
    refill_buffer();
}




/*
  Suppose

  file_size = 1500
  buffer_size = 300
  buffer_offset = 420
  buffer_current-buffer_start = 38

  A request for offset 50 from POS_BEGIN:
  The new buffer_current should now point to a location 370 behind its current value, i.e.
  buffer_current + 50 - 420

  A request for offset 50 from POS_CURRENT:
  The new buffer_current should now point to a location 50 ahead of its current value

  A request for offset 50 from POS_END:
  The new buffer_current should now point to a location that is
  1500 - (


*/


//                                    
//                  buffer_current  offset
//                         |          |
// fs                bs        be                              fe
// |                 |         |                               |
// |<-buffer_offset->|<------->|
//                      ^
//                      |
//                   buffer_size
   

//precondition: buffer_current is within [buffer_start, buffer_end]
//postcondition: buffer_current is set to point to the part of the
//old or possibly newly rewritten buffer corresponding to the requested
//position

/* the semantics of this seekg wrapper should be as a filter to the file based
   seekg.  it should only update buffer_current if the desired new position
   is within the current buffer range.
   if it is not, then this should read from the stream, filling the buffer with
   as much data as possible, either buffer_size or until the end of the file is
   reached.
   in that case, though, the file is not supposed to be read to the end.  it should
   still have good_bit set after the read.

   So, it simply calculates the new buffer_offset if one is needed, and passes
   it on to istream::seekg.  if it is out of bounds, the error will be signaled
   by the istream flags.  

*/
std::istream& 
litestream::seekg(litestream::postype offset, litestream::dirtype way){
  
    //ultimately, calculate a displacement for buffer_current,
    //and decide whether it is within bounds of the start.
    //also, do a check for out-of-bounds, either from the begin or end.
    //the out of bounds can be handled by the seekg of the underlying stream.

    postype read_delta = 0;

    postype read_pos = buffer_current - buffer_start;

    //buffer_offset is just shorthand for tellg() of the underlying
    //istream
    postype stream_pos = buffer_offset + read_pos;

    switch(way){
    case POS_BEGIN: read_delta = offset - stream_pos; break;
    case POS_CURRENT: read_delta = offset; break;
    case POS_END: read_delta = stream_size - stream_pos - offset; break;
    }

    //check whether we need to update the buffer
    //need to update buffer if new read_pos is negative

    //
    postype new_read_pos = read_pos + read_delta;
    postype buffer_delta = 0;
    if (new_read_pos < 0 || new_read_pos > buffer_size){ 
        buffer_delta = new_read_pos;
    }

    //length to read is the min of the buffer_size or
    //whatever's remaining.
    if (buffer_delta != 0){
        buffer_offset = stream_pos + read_delta;

        if (! this->refill_buffer()){
            buffer_current = buffer_start;
            buffer_end = buffer_start;
            return *stream;
        }

    } else {
        buffer_current += read_delta;
        assert(buffer_current >= buffer_start);
        assert(buffer_current <= buffer_end);
    }

    return *stream;
}


//reread the buffer based on a newly updated buffer_offset
//sync buffer_current and buffer_end
bool litestream::refill_buffer(){

    postype refill_length = std::min(buffer_size, stream_size - buffer_offset);

    if (refill_length > 0)
    {
        stream->seekg(buffer_offset, std::ios::beg);
        stream->read(buffer_start, refill_length);
        
        buffer_current = buffer_start;
        buffer_end = buffer_start + refill_length;

        assert(stream->good());
        return true;
    }
    else
    {
        return false;
    }
    //printf("Refilling buffer at offset %"PRId64"\n", buffer_offset);
 
}

//reads from the litestream object, reading as much from the buffer 
//as possible, and repositioning if necessary
std::istream & litestream::read(char * target, postype size){

    postype remain = buffer_end - buffer_current;
    if (size <= remain)
    {
        memcpy(target, buffer_current, size);
        buffer_current += size;
    } 
    else
    {
        memcpy(target, buffer_current, remain);
        buffer_offset += remain + (buffer_current - buffer_start);
        if (! refill_buffer()) return *stream;

        read(target + remain, size - remain);
    }

    return *stream;

}


//get a single character
std::istream & litestream::get(char & c){

    postype remain = buffer_end - buffer_current;
    if (remain > 0){
        c = *buffer_current;
        buffer_current++;
    } else {
        buffer_offset += buffer_current - buffer_start;
        if (! refill_buffer()) return *stream;
        c = *buffer_current;
        buffer_current++;
    }
    return *stream;
}
