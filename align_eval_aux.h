#ifndef _ALIGN_EVAL_AUX_H
#define _ALIGN_EVAL_AUX_H

#include <vector>
#include <zlib.h>
#include <time.h>
#include <cstdio>

#include "sam_index.h"

typedef std::vector<sam_index>::iterator INDEX_ITER;

sam_index **
create_load_ranges(sam_index *line_index, size_t num_threads, size_t num_lines);


//transform a set of <num_threads> loads of index, remapping the
//flowcell ids according to work_remap
void apply_remap(size_t num_threads, sam_index **line_index_loads, 
                 SAM_INDEX_TYPE itype,
                 unsigned int **work_remap);

// initialize line_index_chunk, and augment flowcell_dict
void samlines_to_index(size_t num_threads, 
                       char **samlines,
                       size_t num_lines,
                       SAM_INDEX_TYPE itype,
                       const contig_dict *cdict,
                       sam_index *line_index_chunk,
                       index_dict_t *flowcell_dict);

std::pair<size_t, size_t> 
process_chunk(std::vector<char *> & samlines,
              char * chunk_buffer_in,
              char * chunk_buffer_out,
              size_t num_threads,
              SAM_INDEX_TYPE itype,
              const contig_dict *contig_dict,
              index_dict_t * index_dict, // this will accumulate each time process_chunk is called
              std::vector<sam_index> * line_index);

void
get_key_quantiles(std::vector<sam_index> const& line_index,
                  size_t num_chunks,
                  size_t * key_quantile_sizes,
                  size_t * key_quantile_nlines,
                  std::vector<sam_index> * key_quantile_sentinels);


size_t
set_start_offsets(sam_index *beg, sam_index *end, size_t initial_offset);


void 
write_final_merge(std::vector<sam_index> const& ok_index,
                  std::vector<INDEX_ITER> const& offset_quantiles,
                  std::vector<gzFile> const& tmp_fhs,
                  // std::vector<FILE *> const& tmp_fhs,
                  /* std::vector<z_stream *> const& zstreams, */
                  /* size_t zstream_bufsize, */
                  bool check_uniqueness,
                  FILE * out_dat_fh,
                  FILE * out_ind_fh);


/* std::vector<INDEX_ITER>  */
/* get_quantiles(std::vector<sam_index> * line_index, */
/*               bool (less_fcn)(sam_index const&, sam_index const&), */
/*               size_t num_chunks); */


sam_index const***
get_quantiles(sam_index const** line_index,
              size_t index_size,
              bool (less_fcn)(sam_index const&, sam_index const&),
              size_t num_chunks);



#endif // _ALIGN_EVAL_AUX_H
