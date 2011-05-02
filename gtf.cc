#include "gtf.h"

#include <set>
#include <cstring>

#include "dep/tools.h"

bool GTFEntry::get_next_record(FILE * gtf_fh)
{
    int num_parsed = fscanf(gtf_fh, "%[^\n]\n", this->line);

    if (num_parsed != 1)
    {
        if (feof(gtf_fh))
        {
            return false;
        }
        else
        {
            fprintf(stderr, "GTFEntry::get_next_record: error parsing line\n");
            exit(1);
        }
    }

    if (strncmp(line, "##", 2) == 0)
    {
        //skip comment lines
        this->clear_record();
        this->is_data_line = false;
        return true;
    }

    char score_string[10];

    int parsed_main_fields = 
        sscanf(this->line, "%s\t%s\t%s\t%zu\t%zu\t%s\t%c\t%c\t%[^\n]\n",
               this->seqname, this->source, this->feature, 
               &this->start, &this->end, score_string, 
               &this->strand, &this->frame, this->attribute_string);

    if (parsed_main_fields != 9)
    {
        fprintf(stderr, "Error: GTF file doesn't have 9 tab-separated fields\n"
                "Please see specification at:\n"
                "http://genome.ucsc.edu/FAQ/FAQformat.html#format4"
                "Error at this line:\n%s\n", this->line);
        exit(1);
    }
        
    int parsed_attributes = 
        sscanf(this->attribute_string, "gene_id \"%[^;\"]\"; transcript_id \"%[^;\"]\";",
               this->gene_id, this->transcript_id);

    if (score_string[0] == '.')
    {
        this->score = -1.0;
    }
    else
    {
        sscanf(score_string, "%f", &this->score);
    }
        
    if (parsed_attributes != 2)
    {
        fprintf(stderr, "Error: GTF line has bad attribute string:\n%s\n"
                "Please see specification at %s",
                this->line,
                "http://genome.ucsc.edu/FAQ/FAQformat.html#format4");
        exit(1);
    }
    this->is_data_line = true;
    return true;
}


void GTFEntry::print(FILE * gtf_fh)
{
    fprintf(gtf_fh, "%s\t%s\t%s\t%Zu\t%Zu\t%f\t%c\t%c\t%s\n",
            this->seqname, this->source, this->feature, 
            this->start, this->end, this->score, 
            this->strand, this->frame, this->attribute_string);
}


SequenceProjection::SequenceProjection(char const* _species,
                                       char const* _source_dna,
                                       char const* _target_dna,
                                       char const* _strand_string,
                                       int64_t offset,
                                       char const* cigar_string)
{
    species = std::string(_species);
    source_dna = std::string(_source_dna);
    target_dna = std::string(_target_dna);
    same_strand = _strand_string[0] == '+';

    cigar = Cigar::FromString(cigar_string, offset);
    cigar_index = Cigar::ComputeOffsets(cigar);
}


//order by coordinates of 
bool SequenceProjection::operator<(SequenceProjection const& b) const
{
    
    int64_t a_left_offset = Cigar::LeftOffset(this->cigar, true);
    int64_t b_left_offset = Cigar::LeftOffset(b.cigar, true);

    int64_t a_right_neg_offset = - Cigar::RightOffset(this->cigar, true);
    int64_t b_right_neg_offset = - Cigar::RightOffset(b.cigar, true);

    return this->source_dna < b.source_dna ||
        (this->source_dna == b.source_dna &&
         (a_left_offset < b_left_offset ||
          (a_left_offset == b_left_offset &&
           (a_right_neg_offset < b_right_neg_offset ||
            (a_right_neg_offset == b_right_neg_offset &&
             (this->target_dna < b.target_dna ||
              (this->target_dna == b.target_dna &&
               (this->same_strand < b.same_strand))))))));
}


//use this if source = tx, target = genome.
//orders elements by contig on genome, then by coordinates relative to
//genome, then by transcript name.  This ordering optimizes buffer seeking
//in the genome sequence file.
bool LessSequenceProjectionTarget::operator()(SequenceProjection const& a,
                                              SequenceProjection const& b) const
{
    int64_t a_left_offset = Cigar::LeftOffset(a.cigar, false);
    int64_t b_left_offset = Cigar::LeftOffset(b.cigar, false);
    
    int64_t a_right_offset = Cigar::Length(a.cigar, false);
    int64_t b_right_offset = Cigar::Length(b.cigar, false);

    return a.target_dna < b.target_dna ||
        (a.target_dna == b.target_dna &&
         (a_left_offset < b_left_offset ||
          (a_left_offset == b_left_offset &&
           (a_right_offset < b_right_offset ||
            (a_right_offset == b_right_offset &&
             (a.source_dna < b.source_dna ||
              (a.source_dna == b.source_dna &&
               a.same_strand < b.same_strand)))))));
}



unique_transcript::unique_transcript(char const* _cn, char const* _tid, char _st) :
    contig_name(_cn), transcript_id(_tid), strand(_st) { }

unique_transcript::unique_transcript() : contig_name(""), transcript_id(""), strand('.') { }



bool unique_transcript::operator<(unique_transcript const& ut) const
{
    return this->contig_name < ut.contig_name ||
        (this->contig_name == ut.contig_name &&
         (this->transcript_id < ut.transcript_id ||
          (this->transcript_id == ut.transcript_id &&
           (this->strand < ut.strand))));
}



//ensures the set of exons in a given transcript are non-overlapping
struct non_overlapping_interval
{
    int64_t start_boundary;
    int64_t end_boundary;
    non_overlapping_interval(int64_t _start, int64_t _end) : 
        start_boundary(_start), end_boundary(_end) { }

    non_overlapping_interval() : 
        start_boundary(0), end_boundary(0) { }

    bool operator<(non_overlapping_interval const& interval) const
    {
        return (this->end_boundary <= interval.start_boundary);
    }
};


void GTFEntry::clear_record()
{
    this->line[0] = '\0';
    this->seqname[0] = '\0';
    this->source[0] = '\0';
    this->feature[0] = '\0';
    this->start = 0;
    this->end = 0;
    this->score = 0;
    this->strand = '.';
    this->frame = 0;
    this->attribute_string[0] = '\0';
    this->gene_id[0] = '\0';
    this->transcript_id[0] = '\0';
}


std::map<std::string, std::string>
gtf_to_transcript_gene_map(char const* gtf_file)
{
    FILE * gtf_fh = open_or_die(gtf_file, "r", "Input GTF file");
    GTFEntry ge;
    std::map<std::string, std::string> transcript_gene_map;
    std::string t, g;
    while (ge.get_next_record(gtf_fh))
    {
        t = std::string(ge.transcript_id);
        g = std::string(ge.gene_id);
        transcript_gene_map[t] = g;
    }
    return transcript_gene_map;
}

std::set<SequenceProjection>
gtf_to_sequence_projection(char const* gtf_file, char const* species)
{
    FILE * gtf_fh = open_if_present(gtf_file, "r");
    if (gtf_fh == NULL)
    {
        fprintf(stderr, "Error: couldn't open gtf file '%s' "
                "or none provided\n", gtf_file);
        exit(1);
    }

    typedef std::set<non_overlapping_interval> EXONS;
    typedef std::map<unique_transcript, EXONS> TRANSCRIPTS;
    typedef std::pair<EXONS::iterator, bool> EXONS_INSERTER;

    TRANSCRIPTS transcripts;

    EXONS_INSERTER exon_ins;
    
    GTFEntry ge;

    while (ge.get_next_record(gtf_fh))
    {

        if (strcmp(ge.feature, "exon") != 0)
        {
            continue;
        }

        unique_transcript ut(ge.seqname, ge.transcript_id, ge.strand);
        non_overlapping_interval ni(ge.start - 1, ge.end);

        EXONS & exons = transcripts[ut];
       
        exon_ins = exons.insert(ni);

        if (! exon_ins.second)
        {
            fprintf(stderr, "Error: Found duplicate or overlapping exon\n"
                    "transcript %s on %s (%c strand)\n"
                    "duplicate exon: %Zu to %Zu\n",
                    ge.transcript_id, ge.seqname, ge.strand, ge.start, ge.end);
            exit(1);
        }

    }
    fclose(gtf_fh);

    typedef std::set<SequenceProjection> SP_SET;
    SP_SET sequence_projection;
    std::pair<SP_SET::iterator, bool> sp_inserter;

    TRANSCRIPTS::iterator tx_iter;
    for (tx_iter = transcripts.begin(); tx_iter != transcripts.end(); ++tx_iter)
    {
        unique_transcript const& transcript = (*tx_iter).first;
        EXONS const& exons = (*tx_iter).second;

        assert(! exons.empty());

        //construct a cigar_string from a well-ordered set of exons
        Cigar::CIGAR_VEC cigar;
        EXONS::iterator eit = exons.begin();
        EXONS::iterator last_eit = exons.begin();

        size_t start_on_genome = (*eit).start_boundary;
        cigar.push_back(Cigar::Unit(Cigar::M, 
                                    (*eit).end_boundary - 
                                    (*eit).start_boundary));
        ++eit;
        for (; eit != exons.end(); ++eit, ++last_eit)
        {
            cigar.push_back(Cigar::Unit(Cigar::D,
                                        (*eit).start_boundary 
                                        - (*last_eit).end_boundary));
            cigar.push_back(Cigar::Unit(Cigar::M,
                                        (*eit).end_boundary - 
                                        (*eit).start_boundary));
        }
        
        char * cigar_string = new char[cigar.size() * 10]; //(assume each op can have a 9-digit number
        Cigar::ToString(cigar.begin(), cigar.end(), cigar_string);

        char strand_string[2] = ".";
        strand_string[0] = transcript.strand;

        SequenceProjection 
            sp(species, 
               transcript.contig_name.c_str(), 
               transcript.transcript_id.c_str(),
               strand_string, start_on_genome, cigar_string);

        sp_inserter = sequence_projection.insert(sp);

        if (! sp_inserter.second)
        {
            fprintf(stderr, "Found duplicate transcript:\n"
                    "%s:%c (%s)\n", 
                    sp.source_dna.c_str(),
                    (sp.same_strand ? '+' : '-'),
                    sp.target_dna.c_str());
            exit(1);
        }
        delete cigar_string;
    }

    return sequence_projection;
    
}
