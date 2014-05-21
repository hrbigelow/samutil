// test to see whether scanf is suitable for reading binary files
#include <stdint.h>

// taken from http://sourceforge.net/p/predef/wiki/Endianness/
// also see http://en.wikipedia.org/wiki/Endianness#Endianness_in_files_and_byte_swap

enum {
    ENDIAN_UNKNOWN,
    ENDIAN_BIG,
    ENDIAN_LITTLE,
    ENDIAN_BIG_WORD,   /* Middle-endian, Honeywell 316 style */
    ENDIAN_LITTLE_WORD /* Middle-endian, PDP-11 style */
};

int endianness(void)
{
    uint32_t value;
    uint8_t *buffer = (uint8_t *)&value;

    buffer[0] = 0x00;
    buffer[1] = 0x01;
    buffer[2] = 0x02;
    buffer[3] = 0x03;

    switch (value)
    {
    case UINT32_C(0x00010203): return ENDIAN_BIG;
    case UINT32_C(0x03020100): return ENDIAN_LITTLE;
    case UINT32_C(0x02030001): return ENDIAN_BIG_WORD;
    case UINT32_C(0x01000302): return ENDIAN_LITTLE_WORD;
    default:                   return ENDIAN_UNKNOWN;
    }
}




int main(int argc, char ** argv)
{
    FILE * naked_bam_fh = fopen(argv[1]);
    if (naked_bam_fh == NULL)
    {
        fprintf(stderr, "Error: Couldn't open file %s\n", argv[1]);
        exit(23);
    }
    size_t mem = static_cast<size_t>(atoi(argv[2]));


    // assume, for the purposes of argument, that this is an unzipped
    // bam file.  This is NOT the same as the format that samtools -u
    // outputs.  That format, though it is called 'uncompressed bam'
    // does in fact contain gzip headers, but the data blocks are not
    // compressed.

    // The format that you get if you take a traditional bam file and
    // do gunzip -c on it, I'm calling 'naked bam' format.

    char * header_text;
    char * read_name;
    uint32_t * cigar;
    uint8_t * seq;
    char * qual;
    char * tag;
    char val_type;
    char ** ref_names;
    int32_t * ref_lengths;
    int32_t l_text, n_ref, l_name, l_ref, block_size, ref_id, pos, l_seq, next_ref_id, next_pos, tlen;
    uint32_t bin_mq_nl, flag_nc;

    // slurp up the first part of the file
    char * chunk_buffer_in = new char[mem + 1];
    fread(chunk_buffer_in, 1, mem, naked_bam_fh);

    char * read_pointer = chunk_buffer_in;

    sscanf(read_pointer, "BAM\1%i", &l_text);
    read_pointer += 8; // magic plus 4-byte l_text field

    header_text = read_pointer; // warning: header_text is not expected to be NULL-terminated

    sscanf(read_pointer, "%i", &n_ref);
    read_pointer += 4;

    // assert(n_ref < 100000000); // this should be a reasonable number of reference sequences

    // read the alignment reference information
    ref_names = new char* [n_ref];
    ref_lengths = new int32_t[n_ref];
    for (size_t r = 0; r != n_ref; ++r)
    {
        sscanf(read_pointer, "%i", &l_name);
        ref_names[r] = read_pointer;
        read_pointer += l_name + 1;
        sscanf(read_pointer, "%i", &ref_lengths[r]);
        read_pointer += 4;
    }

    // find the line starts.
    // hack
    while(read_pointer < chunk_buffer_in + mem)
    {
        line_starts.push_back(read_pointer);
        sscanf(read_pointer, "%i", &block_size);
        read_pointer += 4 + block_size;
    }
    
    //now read the records
    // make sure to comply with little-endian enforcement
    for (size_t l = 0; l != line_starts.size(); ++l)
    {
        read_pointer = line_starts[l];
        sscanf(read_pointer, "%i%i%i%u%u%i%i%i%i",
               &block_size, &ref_id, &pos, &bin_mq_ml, &flag_nc, &l_seq, &next_ref_id, &next_pos, &tlen);
        n_cigar_op = flag_nc ^ 65535<<16;
        flag = flag_nc>>16;

        read_pointer += 36; // 9 int32 fields (4 bytes each)
        cigar = read_pointer;

        read_pointer += n_cigar_op * 4;

        for (size_t c = 0; c != n_cigar_op; ++c)
        {
            // parse the cigar?
        }

        seq = read_pointer;
        read_pointer += (l_seq + 1) / 2;
        qual = read_pointer;
    }

    delete ref_names;
    delete ref_lengths;
    delete chunk_buffer_in;

    fclose(naked_bam_fh);
    return 0;
}
