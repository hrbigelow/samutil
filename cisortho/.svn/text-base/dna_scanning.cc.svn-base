#include "dna_scanning.h"
#include "dna.h"

#include "litestream.h"

using namespace cis;

//policy:
/* scan the dna in chunks, but ensure that all windows in the entire length of dna
   are presented to the scan function.  The windows shall be at least minwinsize long
   but as long as is permitted by the dna up to maxwinsize.
   In particular, if the chunk is not at the end of the dna, then every window presented
   will be maxwinsize.  Otherwise, every window up until the last <maxwinsize> nucleotides
   will be maxwinsize.  Then, each successive window will be one nucleotide shorter until
   the final window, which will be <minwinsize>. 

   One caveat:  Since the scan function is in fact comparing a prefix to 


 */

//iter stores 
void IterateDNAWindows(scanStore * scan_fcn,
					   void * store, 
					   int minwinsize, //the smallest window presented to the scan function.  only in effect at ends of DNA
					   int maxwinsize, //the largest window presented to the scan function.
					   int chunk_size, 
					   cis::dna_t const* dna){
	
	char * dnaseq = new char[chunk_size];
	NUC * idnaseq = new NUC[chunk_size];

	litestream & seq = dna->source();
	
	int64_t fh_start = dna->seek_start_pos;
	int64_t end = dna->length();

	seq.seekg(fh_start, litestream::POS_BEGIN);

	int jump = -(maxwinsize - 1);
	uint64_t fh_cur;
	
	while (1){ //break when we go out of bounds
		
		//only back up if this is not the first chunk
		if (seq.tellg() > fh_start) seq.seekg(jump, litestream::POS_CURRENT); //back up to overlap
		fh_cur = seq.tellg();
		
		int64_t chunk_start = fh_cur - fh_start;
		int64_t chunk_end = std::min((chunk_start + (int64_t)chunk_size), end);
		int64_t cur_chunk_size = chunk_end - chunk_start;
		
		seq.read(dnaseq, cur_chunk_size); //get next chunk
		
		TranslateSeq(idnaseq, dnaseq, cur_chunk_size);
		//idnaseq[cur_chunk_size] = cis::X;


		//this test works correctly even when the dna length is a multiple of chunk_size
		bool last_chunk = cur_chunk_size == (end - chunk_start);

		if (last_chunk){
			int lim = cur_chunk_size - minwinsize + 1;
			for (int winpos = 0; winpos != lim; ++winpos){
				int last_winsize = std::min(maxwinsize, (int)(cur_chunk_size - winpos));
				(*scan_fcn)(store, idnaseq + winpos, last_winsize, chunk_start + winpos, dna);
			}


			break; //finished
		}

		else {
			int lim = cur_chunk_size - maxwinsize + 1;
			for (int winpos = 0; winpos != lim; ++winpos)
				(*scan_fcn)(store, idnaseq + winpos, maxwinsize, chunk_start + winpos, dna);
		}

	}

	delete dnaseq;
	delete idnaseq; 

}
