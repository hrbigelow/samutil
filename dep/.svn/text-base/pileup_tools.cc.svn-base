#include "pileup_tools.h"
#include "tools.h"
#include "nucleotide_stats.h"

#include <cassert>

char const PileupSummary::nucleotides[] = "ACGTN.acgtn,";

//maps every character to itself except the letters 'ACGTN.acgtn,'
//are mapped to 'N' (78)
char PileupSummary::code_to_redux[] = {
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 78, 45, 78, 47,
 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
 64, 78, 66, 78, 68, 69, 70, 78, 72, 73, 74, 75, 76, 77, 78, 79,
 80, 81, 82, 83, 78, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
 96, 78, 98, 78,100,101,102, 78,104,105,106,107,108,109, 78,111,
112,113,114,115, 78,117,118,119,120,121,122,123,124,125,126,127,
128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255
};


//load and parse a single line.  Assumes the line has all required fields.
//consumes the current line including newline character.
bool PileupSummary::load_line(FILE * file)
{

    char lookahead;
    lookahead = getc(file);
    ungetc(lookahead, file);
    if (lookahead == EOF)
    {
        return false;
    }
    

    size_t total_read_depth;

    int num_desired_fields = this->_with_soapsnp ? 8 : 4;

    int scanned_fields;
    if (this->_with_soapsnp)
    {
        scanned_fields =
            fscanf(file, "%s\t%i\t%c",
                   this->_reference, &this->_position, 
                   &this->_reference_base);

        if (this->_reference_base == '*')
        {
            //this is an indel site
            scanned_fields += 
                fscanf(file, "\t%*s\t%i\t%i\t%i\t%zu\t",
                       &this->_consensus_quality, &this->_snp_quality,
                       &this->_max_mapping_quality, &total_read_depth);
            scanned_fields++;
            this->_is_indel_line = true;
        }
        else
        {
            //this is a normal consensus pileup site
            scanned_fields += 
                fscanf(file, "\t%c\t%i\t%i\t%i\t%zu\t",
                       &this->_consensus_base,
                       &this->_consensus_quality, &this->_snp_quality,
                       &this->_max_mapping_quality, &total_read_depth);

            this->_is_indel_line = false;
        }            
    }

    else 
    {
        scanned_fields =
            fscanf(file, "%s\t%i\t%c\t%zu\t", this->_reference, 
                   &this->_position, &this->_reference_base, 
                   &total_read_depth);

        this->_is_indel_line = false;
    }


    //to comply with zero or one-based coordinates
    if (this->_ones_offset == 1)
    {
        --this->_position;
    }
    
    if (scanned_fields != num_desired_fields)
    {
        fprintf(stderr, "PileupSummary::load_line: Warning: badly formatted line\n");
        return false;
    }

    if (_quality_codes != NULL &&
        _bases != NULL &&
        _bases_upper != NULL)
    {
        delete _quality_codes;
        delete _bases;
        delete _bases_upper;
    }

    if (this->_do_load_alignment)
    {
        return this->load_read_data(file, total_read_depth);
    }
    else
    {
        //ignore rest of line
        fscanf(file, "%*[^\n]\n");
        return true;
    }
}

bool PileupSummary::load_read_data(FILE * file, size_t total_read_depth)
{
    this->_quality_codes = new char[total_read_depth + 1];
    this->_bases = new char[total_read_depth + 1];
    this->_bases_upper = new char[total_read_depth + 1];
    
    const int max_indel_size = 1000;
    char indel_sequence[max_indel_size];

    char pileup_ccode;
    int indel_size;
    int indel_bin;
    int pileup_value;
    int current_read = 0;
    //int deletion_read = 0;

    while(fscanf(file, "%c", &pileup_ccode) != 0)
    {
        
        //reduce the pileup code
        int pileup_code = pileup_ccode;
        char pileup_redux = PileupSummary::code_to_redux[pileup_code];

        if (pileup_redux == '\n' || pileup_redux == EOF)
        {
            //end of the line or file
            break;
        }
            
        indel_size = 0;
        pileup_value = Nucleotide::base_to_index[pileup_code];
            
        if (pileup_redux == 'N' || pileup_redux == '+' || pileup_redux == '-')
        {
            //actual sequence

            if (pileup_redux == 'N')
            { 
                // ACGTN.acgtn, are mapped to 'N'
                char real_base;
                if (pileup_ccode == '.')
                {
                    real_base = toupper(this->_reference_base);
                }
                else if (pileup_ccode == ',')
                {
                    real_base = tolower(this->_reference_base);
                }
                else
                {
                    real_base = pileup_ccode;
                }

                pileup_value = Nucleotide::base_to_index[static_cast<int>(real_base)];

                this->base_counts[pileup_value]++;
                this->sum_of_counts++;
                indel_bin = 0;

                this->_bases_upper[current_read] = toupper(real_base);
                this->_bases[current_read] = real_base;

                ++current_read;
            }
            
            else 
            {
                //indel
                fscanf(file, "%i", &indel_size);
                fgets(indel_sequence, indel_size + 1, file);
                for (int i = 0; i != indel_size + 1; ++i)
                {
                    indel_sequence[i] = toupper(indel_sequence[i]);
                }
                indel_bin = std::min(this->_indel_histo_size, indel_size) * 
                    (pileup_redux == '+' ? 1 : -1);

                std::map<std::string, int> & indel = this->indel_seqs[indel_bin];
                std::string indel_str(indel_sequence);
                if (indel.find(indel_str) == indel.end())
                {
                    indel.insert(std::make_pair(indel_str, 0));
                }
                indel[indel_str]++;
            }

            //indel size histogram considers matches to be 'zero-length indels'
            this->indel_counts[indel_bin]++;
                
        }

        else if (pileup_redux == '^')
        {
            //beginning of a read
            //ignore the next one character (a character-encoded mapping quality)
            fscanf(file, "%*c");
        }
        else if (pileup_redux == '$')
        {
            //the end of a read.  ignored
        }
        else if (pileup_redux == '*')
        {
            //what is this?  '*' is not documented, but according to
            //Heng Li, represents a deletion.  fprintf(stderr,
            //"asterisk found!\n");
            this->_bases[current_read++] = '*';
        }
        else if (pileup_redux == '\t' ||
                 pileup_redux == ' ')
        {
            //end of read_string, get the qual string
            fscanf(file, " "); //eat white space
            fgets(this->_quality_codes, total_read_depth + 1, file);

        }
        else 
        {
            fprintf(stderr, "Don't know this samtools pileup character code: %c\n", 
                    pileup_redux);
            return false;
        }
            
    }

    //now iterate through the bases and qualities, packing them
    size_t num_seen_gaps = 0;
    for (size_t read_pos = 0; read_pos != total_read_depth; ++read_pos)
    {
        if (this->_bases[read_pos] == '*')
        {
            ++num_seen_gaps;
            continue;
        }
        size_t write_pos = read_pos - num_seen_gaps;
        this->_bases[write_pos] = this->_bases[read_pos];
        this->_bases_upper[write_pos] = this->_bases_upper[read_pos];
        this->_quality_codes[write_pos] = this->_quality_codes[read_pos];
    }

    this->_read_depth = total_read_depth - num_seen_gaps;
    this->_bases[this->_read_depth] = '\0';
    this->_bases_upper[this->_read_depth] = '\0';
    this->_quality_codes[this->_read_depth] = '\0';

    return true;
    
}
