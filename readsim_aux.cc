#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cassert>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <sys/timeb.h>

#include "readsim_aux.h"
#include "dep/nucleotide_stats.h"
#include "dep/tools.h"
#include "dep/stats_tools.h"






/*
SequenceProjection InvertProjection(SequenceProjection const& sp)
{
    char inv_cigar_string[10240];
    Cigar::CIGAR_VEC inv_cigar_vec = Cigar::Invert(sp.cigar, 0);
    Cigar::ToString(inv_cigar_vec.begin(),
                    inv_cigar_vec.end(), 
                    inv_cigar_string);

    char const* strand_string = sp.same_strand ? "+" : "-";

    return SequenceProjection(sp.species.c_str(),
                              sp.target_dna.c_str(),
                              sp.source_dna.c_str(),
                              strand_string,
                              0, inv_cigar_string);
}
*/


bool operator<(Locus const& a, Locus const& b)
{
    return a.species() < b.species() ||
        (a.species() == b.species() &&
         (a.contig_name() < b.contig_name() ||
          (a.contig_name() == b.contig_name() &&
           a.position < b.position)));
}


//loads a gtf file, generating sequence projections for each transcript
//Assumes gtf file obeys UCSC Specification:
/*  	
        From: http://genome.ucsc.edu/FAQ/FAQformat.html#format4

        "
        GTF (Gene Transfer Format) is a refinement to GFF that
        tightens the specification. The first eight GTF fields are the
        same as GFF. The group field has been expanded into a list of
        attributes. Each attribute consists of a type/value
        pair. Attributes must end in a semi-colon, and be separated
        from any following attribute by exactly one space.
        
        The attribute list must begin with the two mandatory attributes:
        
        gene_id value - A globally unique identifier for the genomic
        source of the sequence.

        transcript_id value - A globally unique identifier for the
        predicted transcript.
        "
        
        And from: http://mblab.wustl.edu/GTF2.html

        Textual attributes should be surrounded by doublequotes.

        In addition, the following content restrictions are assumed:

        1.  transcript_id is unique within a contig, strand combination

        2.  within a (transcript_id, contig, strand), exons are
        non-overlapping

  */


LOCUS_SET parse_somatic_mutations(char const* somatic_mutation_file)
{
    LOCUS_SET somatic_mutations;

    FILE * somatic_mutation_fh = open_if_present(somatic_mutation_file, "r");
    if (somatic_mutation_fh != NULL)
    {
        while (! feof(somatic_mutation_fh))
        {
            char species[1000];
            char contig_name[1000];
            size_t position;
            char orig_base;
            char mut_base;
            fscanf(somatic_mutation_fh, "%s\t%s\t%zu\t%c\t%c\n", 
                   species, contig_name, &position, &orig_base, &mut_base);

            somatic_mutations.insert(somatic_mutations.end(), 
                                     Locus(species, contig_name, 
                                           position, orig_base, mut_base));
        }
    }
    close_if_present(somatic_mutation_fh);

    return somatic_mutations;
}


char const* ReadSampler::fastq_fmt = "@%*s\n%*s\n%*s\n%s\n";
char const* ReadSampler::qseq_fmt = "%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%s\n";



void ReadSampler::load_paired_quals(char min_median_qual_code)
{
    char qual1[1000];
    char qual2[1000];
    size_t num_quals_above[2];
    size_t last_rnum_at_reload = SIZE_MAX;

    if (! (this->qual1_buf_file.valid && this->qual2_buf_file.valid))
    {
        fprintf(stderr, "ReadSampler::load_paired_quals: cannot load quals.  "
                "No qual file sources provided\n");
        exit(1);
    }

    char * qual1_line;
    char * qual2_line;
    char * buffer_ptr = this->qual_buffer;

    bool advanced_chunk1, advanced_chunk2;
    bool reloaded_file1, reloaded_file2;

    for (size_t rnum = 0; rnum != this->num_qual_strings; ++rnum)
    {

        qual1_line = this->qual1_buf_file.next_n_lines(this->entry_num_lines, &advanced_chunk1, &reloaded_file1);
        qual2_line = this->qual2_buf_file.next_n_lines(this->entry_num_lines, &advanced_chunk2, &reloaded_file2);

        if (reloaded_file1 != reloaded_file2)
        {
            fprintf(stderr, "ReadSampler::load_paired_quals: "
                    "Qual source files have non-matching lines\n");
            exit(1);
        }

        size_t nscanned1 = sscanf(qual1_line, this->qual_buf_fmt, qual1);
        size_t nscanned2 = sscanf(qual2_line, this->qual_buf_fmt, qual2);

        if (nscanned1 != 1 || nscanned2 != 1)
        {
            if (nscanned1 != 1)
            {
                fprintf(stderr, "Error parsing line in %s:\n%s\n", this->qual1_buf_file.file, qual1_line);
            }            
            if (nscanned2 != 1)
            {
                fprintf(stderr, "Error parsing line in %s:\n%s\n", this->qual2_buf_file.file, qual2_line);
            }        
            exit(1);
        }

        if (reloaded_file1)
        {
            //only use file 1 as the metric.  
            //Assume file 2 is the same size.  It will be flagged if it is not.

            if (rnum == last_rnum_at_reload)
            {
                fprintf(stderr, "No valid qual strings in provided fastq files.\n");
                exit(1);
            }
            last_rnum_at_reload = rnum;
        }

        
        //optionally, test this qualstring pair.
        if (min_median_qual_code > this->zero_quality_code)
        {
            num_quals_above[0] = 0;
            num_quals_above[1] = 0;
            for (size_t q = 0; q != this->read_length; ++q)
            {
                num_quals_above[0] += qual1[q] >= min_median_qual_code ? 1 : 0;
                num_quals_above[1] += qual2[q] >= min_median_qual_code ? 1 : 0;
            }
            if (num_quals_above[0] < this->read_length / 2 ||
                num_quals_above[1] < this->read_length / 2)
            {
                --rnum;
                continue;
            }
        }

        if (strlen(qual1) < this->read_length || strlen(qual2) < this->read_length)
        {
            fprintf(stderr, "Error: load_paired_quals: requested read length of %Zu \n"
                    "exceeds quality string length of provided fastq files %Zu and %Zu\n",
                    this->read_length, strlen(qual1), strlen(qual2));
            exit(1);
        }

        //check 
        strncpy(buffer_ptr, qual1, this->read_length);
        buffer_ptr[this->read_length] = '\0';
        buffer_ptr += this->read_length + 1;
        if (this->do_reverse_second_quals)
        {
            std::reverse(qual2, qual2 + this->read_length);
        }
        strncpy(buffer_ptr, qual2, this->read_length);
        buffer_ptr[this->read_length] = '\0';
        buffer_ptr += this->read_length + 1;
    }

}


void add_somatic_mutations(cis::dna_t const* dna,
                           LOCUS_SET::const_iterator & mutation_locus_start,
                           LOCUS_SET::const_iterator & mutation_locus_end,
                           size_t genome_pos, 
                           size_t genome_block_length,
                           size_t read_pos,
                           char * read_buffer)
{
                           
    assert(mutation_locus_start != mutation_locus_end);
    assert(dna->name == (*mutation_locus_start).contig_name());

    Locus low_query(dna->species().c_str(), dna->name.c_str(), 
                    genome_pos, 'A', 'A');

    Locus high_query(dna->species().c_str(), dna->name.c_str(), 
                     genome_pos + genome_block_length, 'A', 'A');

    LOCUS_SET::const_iterator low_iter =
        std::lower_bound(mutation_locus_start, mutation_locus_end, low_query);

    LOCUS_SET::const_iterator high_iter =
        std::upper_bound(low_iter, mutation_locus_end, high_query);

    for (LOCUS_SET::const_iterator iter = low_iter;
         iter != high_iter; ++iter)
    {
        read_buffer[read_pos + (*iter).position - genome_pos] = 
            (*iter).mutated_base;
    }
}


void get_transformed_sequence(cis::dna_t const* dna,
                              Cigar::CIGAR_VEC const& genome2read,
                              LOCUS_SET::const_iterator mutation_locus_start,
                              LOCUS_SET::const_iterator mutation_locus_end,
                              char * read_buffer,
                              bool zero_terminate)
{

    size_t read_pos = 0;
    size_t genome_pos = 0;
    //int64_t genome_pos = Cigar::LeftOffset(genome2read, true);

    for (Cigar::CIGAR_VEC::const_iterator it = genome2read.begin();
         it != genome2read.end(); ++it)
    {
        if ((*it).op.code == Cigar::M)
        {
            dna->sequence(genome_pos, genome_pos + (*it).length, 
                          read_buffer + read_pos);
        }

        if (mutation_locus_start != mutation_locus_end)
        {
            add_somatic_mutations(dna,
                                  mutation_locus_start, 
                                  mutation_locus_end, 
                                  genome_pos,
                                  (*it).length,
                                  read_pos,
                                  read_buffer);
        }

        genome_pos += Cigar::UnitLength((*it), true);
        read_pos += Cigar::UnitLength((*it), false);
    }
    if (zero_terminate)
    {
        read_buffer[read_pos] = '\0';
    }
}


void fill_quality_distribution(QUALSCORE_DIST * qualscore_dist)
{
    float error_prob;
    double * data_prob;

    for (size_t q = 0; q != NUM_QUAL_SCORES; ++q)
    {
        error_prob = QualityToErrorProb(q);
        for (size_t fb = 0; fb != 4; ++fb) // founder base
        {
            data_prob = (*qualscore_dist)[q][fb];
            std::fill(data_prob, data_prob + 4, (error_prob / 3.0));

            data_prob[fb] = 1.0 - error_prob;
        }
    }
}


//generates a 'called_read' with simulated errors,
//returning the number of errors.
size_t simulate_errors(char const* founder_read,
                       char const* qual_string,
                       QUALSCORE_DIST qualscore_dist,
                       char zero_quality_code,
                       gsl_rng * rand_gen,
                       char * called_read)
{
    char const* fbase;
    char const* qual;
    char * cbase;

    size_t fbase_index;
    size_t cbase_index;

    size_t num_errors = 0;

    int qualscore;

    for (fbase = founder_read, qual = qual_string, cbase = called_read;
         *fbase != '\0'; ++fbase, ++qual, ++cbase)
    {
        fbase_index = 
            Nucleotide::base_to_index[static_cast<size_t>(*fbase)];
    
        qualscore = static_cast<int>(*qual - zero_quality_code);

        assert(qualscore < static_cast<int>(NUM_QUAL_SCORES));

        cbase_index =
            SampleDiscreteDistribution(rand_gen, 
                                       qualscore_dist[qualscore][fbase_index], 1.0);

        *cbase = Nucleotide::bases_upper[cbase_index];
        num_errors += cbase_index == fbase_index ? 0 : 1;
    }
    return num_errors;
}



void set_paired_read_flags(bool sense_strand,
                           int * left_flag,
                           int * right_flag)
{
    *left_flag = 
        SamFlags::MULTI_FRAGMENT_TEMPLATE |
        SamFlags::ALL_FRAGMENTS_MAPPED |
        SamFlags::NEXT_FRAGMENT_ON_NEG_STRAND;

    *right_flag =
        SamFlags::MULTI_FRAGMENT_TEMPLATE |
        SamFlags::ALL_FRAGMENTS_MAPPED |
        SamFlags::THIS_FRAGMENT_ON_NEG_STRAND;

    if (sense_strand)
    {
        //transcript is POS stranded
        *left_flag |= SamFlags::FIRST_FRAGMENT_IN_TEMPLATE;
        *right_flag |= SamFlags::LAST_FRAGMENT_IN_TEMPLATE;
    }
    else
    {
        //transcript is NEG stranded
        *right_flag |= SamFlags::FIRST_FRAGMENT_IN_TEMPLATE;
        *left_flag |= SamFlags::LAST_FRAGMENT_IN_TEMPLATE;

    }
}


/*
void print_paired_fastq_entries(FILE * first_fh, FILE * second_fh,
                                SamLine const* left_read,
                                SamLine const* right_read)
{
    assert(left_read->multi_fragment_template() &&
           left_read->all_fragments_mapped() &&
           right_read->multi_fragment_template() &&
           right_read->all_fragments_mapped());

    assert(left_read->first_fragment_in_template() != 
           right_read->first_fragment_in_template());

    size_t read_length = left_read->raw_read_length();

    char * seq_buffer = new char[read_length + 1];
    char * qual_buffer = new char[read_length + 1];

    seq_buffer[read_length] = '\0';
    qual_buffer[read_length] = '\0';

    SamLine const* first_read = 
        left_read->first_fragment_in_template() ? left_read : right_read;

    SamLine const* second_read =
        left_read->second_read_in_pair() ? left_read : right_read;

    assert(first_read != second_read);

    FILE * second_target_fh = second_fh == NULL ? first_fh : second_fh;

    char const* seq_to_print;
    char const* qual_to_print;

    if (! first_read->this_fragment_on_pos_strand())
    {
        cis::ReverseComplement(first_read->seq, read_length, seq_buffer);
        strcpy(qual_buffer, first_read->qual);
        std::reverse(qual_buffer, qual_buffer + read_length);
        seq_to_print = seq_buffer;
        qual_to_print = qual_buffer;
    }
    else
    {
        seq_to_print = first_read->seq;
        qual_to_print = first_read->qual;
    }

    fprintf(first_fh, "@%s\n%s\n+\n%s\n", first_read->qname,
            seq_to_print, qual_to_print);

    if (! second_read->this_fragment_on_pos_strand())
    {
        cis::ReverseComplement(second_read->seq, read_length, seq_buffer);
        strcpy(qual_buffer, second_read->qual);
        std::reverse(qual_buffer, qual_buffer + read_length);
        seq_to_print = seq_buffer;
        qual_to_print = qual_buffer;
    }
    else
    {
        seq_to_print = second_read->seq;
        qual_to_print = second_read->qual;
    }

    fprintf(second_target_fh, "@%s\n%s\n+\n%s\n", second_read->qname,
            seq_to_print, qual_to_print);
    
    delete seq_buffer;
    delete qual_buffer;

}
*/

ReadSampler::ReadSampler(LOCUS_SET const& _sm,
                         char const* q1_file, char const* q2_file,
                         size_t _rl, char _zqc,
                         size_t max_mem,
                         bool _do_reverse_second_quals,
                         bool _do_blind_read_names) :
    somatic_mutations(_sm),
    read_counter(0),
    read_length(_rl),
    qual_string_index(0),
    num_qual_strings(10000),
    zero_quality_code(_zqc),
    do_reverse_second_quals(_do_reverse_second_quals),
    do_blind_read_names(_do_blind_read_names),
    qual1_buf_file(q1_file, max_mem / 2, 1000),
    qual2_buf_file(q2_file, max_mem / 2, 1000)
{ 
    this->qual_buffer = new char[(this->read_length + 1) * 2 * this->num_qual_strings];
    for (size_t p = 0; p != 2; ++p)
    {
        this->seq_raw[p] = new char[this->read_length + 1];
        this->seq_sim[p] = new char[this->read_length + 1];
    }
    this->perfect_qual = new char[this->read_length + 1];
    std::fill(this->perfect_qual, 
              this->perfect_qual + this->read_length, 
              static_cast<char>(this->zero_quality_code + 40));
    this->perfect_qual[this->read_length] = '\0';

    somatic_mutations_start = somatic_mutations.end();
    somatic_mutations_end = somatic_mutations.end();
}

void ReadSampler::Initialize()
{
    fill_quality_distribution(& this->qualscore_distribution);    
    this->rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(this->rand_gen, millitime.millitm);

    //initialize buffered qual files
    bool file1_success = this->qual1_buf_file.initialize();
    bool file2_success = this->qual2_buf_file.initialize();

    if (! (file1_success && file2_success))
    {
        return;
    }

    //determine file format
    bool advanced_chunk, reloaded_file;
    char * test_line1 = this->qual1_buf_file.next_n_lines(4, &advanced_chunk, &reloaded_file);
    this->qual2_buf_file.next_n_lines(4, &advanced_chunk, &reloaded_file);

    char qual[1000];
    size_t fastq_nfields_read = sscanf(test_line1, fastq_fmt, qual);
    size_t qseq_nfields_read = sscanf(test_line1, qseq_fmt, qual);

    if (qseq_nfields_read == 1)
    {
        this->qual_buf_fmt = ReadSampler::qseq_fmt;
        this->entry_num_lines = 1;
    }
    else if (fastq_nfields_read == 1)
    {
        this->qual_buf_fmt = ReadSampler::fastq_fmt;
        this->entry_num_lines = 4;
    }
    else
    {
        fprintf(stderr, "Error: Unknown qual source format.  Must be fastq or qseq\n");
        exit(1);
    }
}

ReadSampler::~ReadSampler()
{
    delete this->qual_buffer;
    delete this->seq_raw[0];
    delete this->seq_raw[1];
    delete this->seq_sim[0];
    delete this->seq_sim[1];
    delete this->perfect_qual;
    gsl_rng_free(this->rand_gen);
}

void 
ReadSampler::set_somatic_mutation_range(LOCUS_SET::const_iterator start,
                                        LOCUS_SET::const_iterator end)
{
    this->somatic_mutations_start = start;
    this->somatic_mutations_end = end;
}

std::pair<SamLine *, SamLine *>
ReadSampler::sample_pair(SequenceProjection const& transcript_proj,
                         cis::dna_t const* target_dna,
                         size_t transcript_start_pos,
                         size_t transcript_end_pos,
                         size_t min_median_qual_score,
                         bool do_simulate_errors,
                         bool do_sample_sense_strand)
{
    size_t mapping_quality = 255;
    size_t outer_mate_dist = transcript_end_pos - transcript_start_pos;

    const bool zero_terminate = true;
    char * left_read_seq_print;
    char * right_read_seq_print;

    char * left_read_qual_print;
    char * right_read_qual_print;

    size_t num_left_errors;
    size_t num_right_errors;

    char min_median_qual_code = this->zero_quality_code + min_median_qual_score;

    Cigar::CIGAR_VEC transcript2left_read;
    transcript2left_read.push_back(Cigar::Unit(Cigar::Ops[Cigar::D], transcript_start_pos));
    transcript2left_read.push_back(Cigar::Unit(Cigar::Ops[Cigar::M], read_length));

    Cigar::CIGAR_VEC genome2left_read =
        Cigar::TransitiveMerge(transcript_proj.cigar, 
                               transcript_proj.cigar_index,
                               transcript2left_read, false, false);

    get_transformed_sequence(target_dna, genome2left_read,
                             this->somatic_mutations_start, 
                             this->somatic_mutations_end,
                             this->seq_raw[0], zero_terminate);

    Cigar::CIGAR_VEC transcript2right_read;
    transcript2right_read.push_back
        (Cigar::Unit(Cigar::D, transcript_end_pos - read_length));

    transcript2right_read.push_back(Cigar::Unit(Cigar::Ops[Cigar::M], read_length));

    Cigar::CIGAR_VEC genome2right_read = 
        Cigar::TransitiveMerge(transcript_proj.cigar, 
                               transcript_proj.cigar_index,
                               transcript2right_read, false, false);

    get_transformed_sequence(target_dna, genome2right_read,
                             this->somatic_mutations_start, 
                             this->somatic_mutations_end,
                             this->seq_raw[1], zero_terminate);

    //generate mutations if we have valid fastq file inputs.
    if (do_simulate_errors)
    {
        if (this->qual_string_index == num_qual_strings)
        {
            this->qual_string_index = 0;
        }

        if (this->qual_string_index == 0)
        {
            this->load_paired_quals(min_median_qual_code);
            ++this->read_counter;
        }
        size_t qual_buffer_offset = this->qual_string_index * 2 * (read_length + 1);

        this->qual[0] = this->qual_buffer + qual_buffer_offset;
        this->qual[1] = this->qual_buffer + qual_buffer_offset + (read_length + 1);
                        
        num_left_errors = 
            simulate_errors(this->seq_raw[0], this->qual[0],
                            this->qualscore_distribution,
                            zero_quality_code,
                            this->rand_gen, this->seq_sim[0]);
        
        num_right_errors = 
            simulate_errors(this->seq_raw[1], this->qual[1],
                            this->qualscore_distribution, 
                            zero_quality_code,
                            this->rand_gen, this->seq_sim[1]);

        left_read_qual_print = this->qual[0];
        right_read_qual_print = this->qual[1];

        left_read_seq_print = this->seq_sim[0];
        right_read_seq_print = this->seq_sim[1];

        ++this->qual_string_index;
        ++this->read_counter;
    }
    else
    {
        left_read_seq_print = this->seq_raw[0];
        right_read_seq_print = this->seq_raw[1];
        left_read_qual_print = perfect_qual;
        right_read_qual_print = perfect_qual;
        ++this->read_counter;
    }

    Cigar::CIGAR_VEC genome2left_read_trimmed = Cigar::Trim(genome2left_read, false);

    char cigar_string1[1000];
    Cigar::ToString(genome2left_read_trimmed.begin(), 
                    genome2left_read_trimmed.end(), 
                    cigar_string1);

    Cigar::CIGAR_VEC genome2right_read_trimmed = 
        Cigar::Trim(genome2right_read, false);

    char cigar_string2[1000];
    Cigar::ToString(genome2right_read_trimmed.begin(), 
                    genome2right_read_trimmed.end(), 
                    cigar_string2);
                    
    int64_t position1 = Cigar::LeftOffset(genome2left_read, true);
    int64_t position2 = Cigar::LeftOffset(genome2right_read, true);

    //buffer the sam lines, keeping them in sorted order.
    int64_t isize = outer_mate_dist;

    bool sampling_strand = do_sample_sense_strand ? 
        transcript_proj.same_strand :
        (! transcript_proj.same_strand);

    int left_read_flag, right_read_flag;

    set_paired_read_flags(sampling_strand, &left_read_flag, &right_read_flag);

    char qname[1000];

    if (this->do_blind_read_names)
    {
        sprintf(qname, "%zu", this->read_counter);
    }
    else
    {
        char left_strand = 
            (left_read_flag & SamFlags::THIS_FRAGMENT_ON_NEG_STRAND) == 0 ? '+' : '-';

        char right_strand = 
            (right_read_flag & SamFlags::THIS_FRAGMENT_ON_NEG_STRAND) == 0 ? '+' : '-';

        char const* left_read_string = sampling_strand ? "read1" : "read2";
        char const* right_read_string = sampling_strand ? "read2" : "read1";

            sprintf(qname, "%zu:%s:%s:%c:%zu:%s:%zu:%s:%s:%c:%zu:%s:%zu:fragment_size:%zu",
                    this->read_counter, 
                    left_read_string,
                    target_dna->name.c_str(),
                    left_strand, position1 + 1, cigar_string1, num_left_errors,
                    right_read_string,
                    target_dna->name.c_str(),
                    right_strand, position2 + 1, cigar_string2, num_right_errors,
                    isize);
    }

    char const* tag_string1 = "";
    SamLine * left_read = 
        new SamLine(DATA_LINE, qname, left_read_flag, target_dna->name.c_str(),
                    position1, mapping_quality, cigar_string1,
                    "=", position2, isize,
                    left_read_seq_print, left_read_qual_print, tag_string1);

    char const* tag_string2 = "";
    SamLine * right_read =
        new SamLine(DATA_LINE, qname, right_read_flag, target_dna->name.c_str(),
                    position2, mapping_quality, cigar_string2,
                    "=", position1, -1 * isize,
                    right_read_seq_print, right_read_qual_print, tag_string2);

    assert(left_read->qname != NULL);
    return std::make_pair(left_read, right_read);
}
