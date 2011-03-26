#ifndef _LITESTREAM_H
#define _LITESTREAM_H

#include <istream>
#include <fstream>
#include <stdint.h>

/*
A class for fast forward-seek access of a file using a buffer.
This is a wrapper for istream.  istream uses a buffer to avoid frequent
disk access.  However, the buffer is only used for successive reads from
the same position.  Once seekg() is called, the buffer is automatically
reloaded from the file, regardless of whether the new position falls within
the buffer or not.

litestream wraps seekg() and read() with its own separate buffer.  With litestream,
each call to read() first checks whether the region of the stream requested
is covered by the buffer or whether a buffer refill is necessary.  seekg(), similarly
checks whether the requested get pointer position falls within the buffer and
only calls refill_buffer if the requested position falls outside.

If a buffer size is not provided, litestream will attempt to buffer the
entire file.  Upon failed allocation, it will attempt half that much,
and so on, until allocation succeeds
 */


class litestream {

 public:

  enum dirtype {
    POS_BEGIN,
    POS_CURRENT,
    POS_END
  };

  typedef int64_t postype;

  std::istream * stream;

  char * buffer_start;
  char * buffer_current;
  char * buffer_end;

  postype buffer_size;
  postype buffer_offset; //offset of the beginning of the buffer
  postype stream_size;
  postype request_size;

  //request a buffer of size.  attempts to allocate that size,
  //trying half the size at each unsuccessful attempt
  litestream(std::istream & base, postype size);
  litestream(std::istream & base);

  inline ~litestream() { if (this->is_open()) { this->close(); } }

  void open(const char * filename);

  inline void close() { 
    delete buffer_start;
    buffer_start = NULL;
    return static_cast<std::ifstream *>(stream)->close(); 
  }

  inline bool is_open() { 
    return static_cast<std::ifstream *>(stream)->is_open(); 
  }

  inline bool good() { return stream->good(); }

  bool refill_buffer();

  std::istream & read(char * dest, postype size);
  std::istream & seekg(postype off, dirtype way);

  inline postype tellg() { 
    return buffer_offset + buffer_current - buffer_start;
  }
  
  std::istream & get(char & c);

  inline std::ios::iostate rdstate() const { return stream->rdstate(); }

};


#endif // _LITESTREAM_H
