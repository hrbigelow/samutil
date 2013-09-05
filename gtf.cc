#include "gtf.h"

#include "seq_projection.h"

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

    if (strncmp(line, "#", 1) == 0)
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
        

    char * att1 = strstr(this->attribute_string, "gene_id");
    char * att2 = strstr(this->attribute_string, "transcript_id");
    int parsed_attributes = 0;

    parsed_attributes += att1 == NULL ? 0 : sscanf(att1, "gene_id \"%[^;\"]\";", this->gene_id);
    parsed_attributes += att2 == NULL ? 0 : sscanf(att2, "transcript_id \"%[^;\"]\";", this->transcript_id);
    
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


//use this if source = tx, target = genome.
//orders elements by contig on genome, then by coordinates relative to
//genome, then by transcript name.  This ordering optimizes buffer seeking
//in the genome sequence file.
bool LessSequenceProjectionTarget::operator()(SequenceProjection const& a,
                                              SequenceProjection const& b) const
{
    int64_t a_left_offset = a.transformation.empty() ? 0 : a.transformation[0].jump_length;
    int64_t b_left_offset = b.transformation.empty() ? 0 : b.transformation[0].jump_length;
    
    return a.target_dna < b.target_dna ||
        (a.target_dna == b.target_dna &&
         (a_left_offset < b_left_offset ||
          (a_left_offset == b_left_offset &&
           (a.source_dna < b.source_dna ||
            (a.source_dna == b.source_dna &&
             a.same_strand < b.same_strand)))));
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
        EXONS::const_iterator it;

        assert(! exons.empty());

        SequenceProjection 
            sp(species, 
               transcript.transcript_id.c_str(),
               transcript.contig_name.c_str(), 
               transcript.strand, 
               std::vector<block_offsets>());

        size_t last_end_boundary = 0;

        for (it = exons.begin(); it != exons.end(); ++it)
        {
            sp.transformation.push_back
                (block_offsets((*it).start_boundary - last_end_boundary,
                               (*it).end_boundary - (*it).start_boundary));

            last_end_boundary = (*it).end_boundary;
            sp.total_block_length += (*it).end_boundary - (*it).start_boundary;
        }
        
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
    }
    return sequence_projection;
}
