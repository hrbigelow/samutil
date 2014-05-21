/*
  Work like samtools view, allowing conversion from sam to bam and
  back, but using parallel compression / decompression

  Overall, this would work as follows:

  For BAM file input:
  1. Base layer:  Slurp up raw chunks of the file without doing any processing.
  2. Traverse the binary, finding start positions of the BGZF blocks by reading the EXTRA fields.

  3. Do parallel decompression into a buffer (using ULEN, the
  uncompressed length field to determine output buffer size

  4. Serially, traverse the uncompressed output, finding line starts based on the BAM format.
  5. In parallel, convert to SAM structure, copying to memory.

  6. In serial, write out the converted structure.

  For SAM file input:
  1. Base layer: Slurp up raw chunk of the file.
  2. Find all line starts, and nullify (find_complete_lines_nullify())
  3. Convert to BAM format and copy to buffer.
  4. In parallel, compress BGZF_MAX_BLOCK_SIZE chunks, wrapping in BGZF headers.
  5. Write to output stream.

 */

int sam_view_usage(size_t mdef, size_t zdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil view [OPTIONS] reads.bam reads.sam\n"
            " -- or --\n"
            "samutil view [OPTIONS] reads.sam reads.bam\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-S  FLAG      if present, assume SAM input and produce BAM output (otherwise the reverse)\n"
            "-z  INT       zlib compression level, if producing BAM [0-8] [%Zu]\n",
            mdef, zdef);

    return 1;
}


int main_sam_view(int argc, char** argv)
{
    size_t gz_level_def = 5;
    size_t gz_level = gz_level_def;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;
    bool sam_input = false;

    char c;
    while ((c = getopt(argc, argv, "m:t:z:S")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'S': sam_input = true; break;
        case 'z': gz_level = static_cast<size_t>(atof(optarg)); break;
        default: return sam_view_usage(max_mem_def, gz_level_def);
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    __gnu_parallel::_Settings psettings;
    psettings.for_each_minimal_n = 2;
    __gnu_parallel::_Settings::set(psettings);

    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return sam_sort_usage(max_mem_def, gz_level_def);
    }

    if (gz_level > 8)
    {
        fprintf(stderr, "Error: zlib compression level > 8 is not supported\n");
        return sam_sort_usage(max_mem_def, gz_level_def);
    }

    char const* in_file = argv[optind];
    char const* out_file = argv[optind + 1];

    FILE * in_fh = open_if_present(in_file, "r");
    FILE * out_fh = open_if_present(out_fh, "w");

    // memory
    char * chunk_buffer_in;
    char * chunk_buffer_out;

    if (sam_input)
    {
        //   For SAM file input:
        //   1. Base layer: Slurp up raw chunk of the file.
        //   2. Find all line starts, and nullify (find_complete_lines_nullify())
        //   3. Convert to BAM format and copy to buffer.
        //   4. In parallel, compress BGZF_MAX_BLOCK_SIZE chunks, wrapping in BGZF headers.
        //   5. Write to output stream.
       

    }
    else
    {
        // For BAM file input:

        // 1. Base layer: Slurp up raw chunks of the file without
        // doing any processing.

        // 2. Traverse the binary, finding start positions of the BGZF
        // blocks by reading the EXTRA fields.

        // 3. Do parallel decompression into a buffer (using ULEN, the
        // uncompressed length field to determine output buffer size

        // 4. Serially, traverse the uncompressed output, finding line
        // starts based on the BAM format.

        // 5. In parallel, convert to SAM structure, copying to memory.
        // 6. In serial, write out the converted structure.
        

    }

    close_if_present(in_fh);
    close_if_present(out_fh);

    return 0;
}

