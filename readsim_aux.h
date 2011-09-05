#include <string>
#include <set>
#include <utility>
#include <gsl/gsl_rng.h>

#include "cisortho/memory.h"
#include "cigar_ops.h"
#include "cisortho/dna.h"
#include "sam_helper.h"
#include "gtf.h"
#include "file_utils.h"


class Locus
{
 public:
    memory::cstore<std::string> species;
    memory::cstore<std::string> contig_name;
    size_t position;
    char original_base;
    char mutated_base;
 Locus(char const* _species, char const* _contig, 
       size_t _pos, char _orig, char _mut) :
    position(_pos), original_base(_orig), mutated_base(_mut)
    {
        species(std::string(_species));
        contig_name(std::string(_contig));
    }
    friend bool operator<(Locus const& a, Locus const& b);
};



struct LessLocus
{
    bool operator()(Locus const& a, Locus const& b) const
    {
        return a.species() < b.species() ||
            (a.species() == b.species() &&
             (a.contig_name() < b.contig_name() ||
              (a.contig_name() == b.contig_name() &&
               a.position < b.position)));
    }
};

typedef std::set<Locus, LessLocus> LOCUS_SET;



std::set<SequenceProjection> 
load_sequence_projection(char const* coord_transform_file);


LOCUS_SET parse_somatic_mutations(char const* somatic_mutation_file);


void get_transformed_sequence(cis::dna_t const* dna,
                              Cigar::CIGAR_VEC const& target2read,
                              LOCUS_SET::const_iterator mutation_locus_start,
                              LOCUS_SET::const_iterator mutation_locus_end,
                              char * read_buffer,
                              bool zero_terminate);

size_t const NUM_QUAL_SCORES = 50;

typedef double QUALSCORE_DIST[NUM_QUAL_SCORES][4][4];

void fill_quality_distribution(QUALSCORE_DIST * qualscore_dist);

size_t simulate_errors(char const* founder_read,
                       char const* qual_string,
                       QUALSCORE_DIST qualscore_dist,
                       char zero_quality_code,
                       gsl_rng * rand_gen,
                       char * called_read);


enum ReadType {
    FIRST_READ,
    SECOND_READ,
    SINGLE_READ
};


void simulate_read(ReadType read_type);


void set_paired_read_flags(bool sense_strand,
                           int * flag_first_in_pair,
                           int * flag_second_in_pair);



//print first and second reads to fastq file handles.
//assume the two reads conform to SAM format spec and that the
//seq field corresponds to the sense strand of the reference.
//left and right reads are assumed paired in sequencing and mapped in
//the proper pair, 
void print_paired_fastq_entries(FILE * read1_fh, FILE * read2_fh,
                                SamLine const* left_read,
                                SamLine const* right_read);


// class for sampling reads from a 

class ReadSampler
{
    LOCUS_SET somatic_mutations;
    LOCUS_SET::const_iterator somatic_mutations_start;
    LOCUS_SET::const_iterator somatic_mutations_end;
    size_t read_counter;
    size_t read_length;
    size_t qual_string_index;
    size_t num_qual_strings;
    char zero_quality_code;
    char * qual_buffer;
    char * seq_raw[2];
    char * seq_sim[2];
    char * qual[2];
    char * seq_print[2];
    char * qual_print[2];
    char * perfect_qual;
    bool do_reverse_second_quals;
    bool do_blind_read_names;
    QUALSCORE_DIST qualscore_distribution;
    gsl_rng * rand_gen;

    static char const* fastq_fmt;
    static char const* qseq_fmt;

 public:

    BufferedFile qual1_buf_file;
    BufferedFile qual2_buf_file;
    char const* qual_buf_fmt;
    size_t entry_num_lines; // 1 for qseq, 4 for fastq

    ReadSampler(LOCUS_SET const& somatic_mutations,
                char const* q1_file, char const* q2_file,
                size_t read_length, char zero_quality_code,
                size_t max_mem,
                bool do_reverse_second_quals,
                bool do_blind_read_names);

    ~ReadSampler();

    void Initialize();

    std::pair<SamLine *, SamLine *>
        sample_pair(SequenceProjection const& transcript_proj,
                    cis::dna_t const* target_dna,
                    size_t transcript_start_pos,
                    size_t transcript_end_pos,
                    size_t min_median_qual_score,
                    bool simulate_errors,
                    bool do_sample_sense_strand);

    void set_somatic_mutation_range(LOCUS_SET::const_iterator start,
                                    LOCUS_SET::const_iterator end);
    
    void load_paired_quals(char min_median_qual_code);
};
