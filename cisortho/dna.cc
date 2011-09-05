#include <algorithm>
#include <sstream>
#include <cctype>
#include <cassert>
#include <cstdio>
#include <cstring>

#include "dna.h"

using std::string;

namespace cis 
{

    string dna_t::iupac =      "XACMGRSVTWYHKDBN";
    string dna_t::iupac_comp = "XTGKCYSBAWRDMHVN";

    char base_to_complement[] =
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XTXGXXXCXXXXXXNXXXXXAXXXXXXXXXXX"
        "XtXgXXXcXXXXXXnXXXXXaXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

    char base_to_uppercase_complement[] = 
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XTXGXXXCXXXXXXNXXXXXAXXXXXXXXXXX"
        "XTXGXXXCXXXXXXNXXXXXAXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
	
    //operators for dna_t
    bool operator==(dna_t const& d1, dna_t const& d2) {
        return 
            d1.species() == d2.species() &&
            d1.name == d2.name &&
            d1.length() == d2.length();
    }


    bool operator!=(dna_t const& d1, dna_t const& d2) { return ! (d1 == d2); }


    bool operator<(dna_t const& d1, dna_t const& d2) { 
        return 
            d1.species() < d2.species() || 
            (d1.species() == d2.species() && 
             (d1.name < d2.name));
    }



    dna_t::dna_t(std::string const& _path, 
                 string const& _species, 
                 string const& _name, 
                 int64_t const _seek_start_pos, 
                 int64_t const _seek_end_pos) : 
        seek_start_pos(_seek_start_pos), 
        seek_end_pos(_seek_end_pos) {

        species(_species);
        source(named_stream(_path));

        if (seek_start_pos > seek_end_pos) {
            cerr<<"start position ("<<seek_start_pos
                <<") should be less than end position ("
                <<seek_end_pos<<")!"<<endl;
            exit(86);
        }
        name = _name;
    }


    dna_t::dna_t(istream & i, 
                 std::string const& dna_directory){
        string sp;
        std::string filename;
        i>>sp
         >>name
         >>seek_start_pos
         >>seek_end_pos
         >>filename;
        i.seekg(1,std::ios::cur); //go past the newline
		
        species(sp);
        source(dna_t::named_stream
               (std::string(dna_directory 
                            + std::string("/")
                            + filename)));

        //source().safe_open();
		
    }
		
    dna_t::dna_t(std::string const& _species,
                 std::string const& _contig_name,
                 int64_t const _seek_start_pos, 
                 int64_t const _seek_end_pos) :
        seek_start_pos(_seek_start_pos),
        seek_end_pos(_seek_end_pos)
    {
        this->species(_species);
        this->name = _contig_name;
    }
    
                 


    bool less_dna_ptr::operator()(cis::dna_t const* d1, 
                                  cis::dna_t const* d2) const { 
        return *d1 < *d2; 
    }

    bool less_dna::operator()(dna_t const& d1, 
                              dna_t const& d2) const { 
        return d1 < d2; 
    }

    ostream& operator<<(ostream &o, dna_t const& d){
        o<<d.species()<<'\t'
         <<d.name<<'\t'
         <<d.seek_start_pos<<'\t'
         <<d.seek_end_pos<<'\t'
         <<d.source().p;
        return o;
    }




    istream& operator>>(istream &i, dna_t & d){

        string sp;
        std::string path;
        i>>sp
         >>d.name
         >>d.seek_start_pos
         >>d.seek_end_pos
         >>path;

        d.species(sp);
        d.source(dna_t::named_stream(path));
        d.source().safe_open();
        return i;
    }


		

    cis::dnastring RevCompDNA(cis::dnastring const& seq){

        cis::dnastring rev(seq);
        std::reverse(rev.begin(), rev.end());

        int len=seq.length(),i;

        for (i=0; i < len; i++) rev[i] = CompNUC(rev[i]);

        return rev;
    }


    //preserves case while reverse complementing
    string RevCompStr(string const& seq){
        string rev(seq);
        std::reverse(rev.begin(), rev.end());
        int len=seq.length(),i;
        for (i=0; i < len; i++) rev[i] = CompChar(rev[i]);
        return rev;
    }





    //This reverses the bit order of first four bits
    NUC CompNUC(NUC n) { 
        static NUC comp[] = { X,T,G,K,C,Y,S,B,A,W,R,D,M,H,V,N };
        return comp[n];
    }



    char Nuc2char(NUC d){
        if (d >= 16){
            std::fprintf(stderr, "Nuc2char: non-nucleotide stored.");
            assert(false);
        }
        static char chars[] = "XACMGRSVTWYHKDBN*";
        return chars[static_cast<int>(d)];
    }


    char Nuc2lcchar(NUC d){
        if (d >= 16){
            std::fprintf(stderr, "Nuc2char: non-nucleotide stored.");
            assert(false);
        }
        static char chars[] = "xacmgrsvtwyhkdbn*";
        return chars[static_cast<int>(d)];
    }

    NUC char2NUC(unsigned char c){

        static NUC nucs[] = {
            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z,
            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z,
            z,A,B,C,D,z,z,G, H,z,z,K,z,M,N,z, z,z,R,S,T,z,V,W, X,Y,z,z,z,z,z,z,
            z,A,B,C,D,z,z,G, H,z,z,K,z,M,N,z, z,z,R,S,T,z,V,W, X,Y,z,z,z,z,z,z,

            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z,
            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z,
            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z,
            z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z, z,z,z,z,z,z,z,z
        };

        static std::set<char> seen_chars;
        if (nucs[static_cast<int>(c)] == z){
            if (seen_chars.find(c) == seen_chars.end()){
                seen_chars.insert(c);
                std::fprintf(stderr, "char2NUC: nuc %c converted to X.  Further warnings disabled...\n", c);
            }
            return X;
        }
        else return nucs[static_cast<int>(c)];
    }


    //complements the nucleotide or consensus, preserving case.
    char CompChar(char c){
        int lower_adjust = 
            static_cast<unsigned int>(c - 'a') <= 
            static_cast<unsigned int>('z' - 'a') ?
            'a' - 'A' : 0;
        char comp = Nuc2char(CompNUC(char2NUC(c))); //will be upper case
        return comp + lower_adjust;
    }


    void TranslateSeq(NUC *const iseq, char const*const seq, int const size){
        for (int i=0; i < size; i++) iseq[i] = char2NUC(seq[i]);
    }


    string FromDNAString(cis::dnastring const& d){

        string s;
        s.resize(d.size());
        for (unsigned int i=0; i < s.size(); i++) s[i] = Nuc2char(d[i]);
        return s;
    }

    string FromDNAString(cis::NUC const*const sequence, int size){

        string s;
        s.resize(size);
        for (size_t i=0; i < s.size(); i++) s[i] = Nuc2char(sequence[i]);
        return s;
    }
    cis::dnastring ToDNAString(string const s){

        cis::dnastring d;
        d.resize(s.size());
        for (unsigned int i=0; i < d.size(); i++) d[i] = char2NUC(s[i]);
        return d;
    }


    cis::dnastring ToDNAString(char const*const s, int sz){
        cis::dnastring d;
        d.resize(sz);
        for (int i=0; i < sz; i++) d[i] = char2NUC(s[i]);
        return d;
    }

    //output the actual sequence associated with the region
    string dna_t::sequence(uint64_t start, uint64_t end) const {
	
        std::ostringstream os;

        int64_t size = end - start;
        if (size <= 0){
            cerr<<"end position of "<<end<<" must be strictly greater than start position "<<start<<endl;
            exit(63);
        }
        else if (size > 10000000) {
            cerr<<"A maximum size of 10,000,000 is allowed for outputting raw dna sequence"<<endl;
            exit(57);
        }

        else {
            char *buf = new char[size];
            source().seekg(seek_start_pos + start, litestream::POS_BEGIN);
            source().read(buf, size);
            os.write(buf, size);
            delete buf;
        }
        return os.str();
    }


    //output the actual sequence associated with the region
    void dna_t::sequence(uint64_t start, uint64_t end, char * sequence) const {
	
        if (end <= start){
            cerr<<"end position of "<<end<<" must be strictly greater than start position "<<start<<endl;
            exit(63);
        }
        uint64_t size = end - start;
        if (size > 10000000) {
            cerr<<"A maximum size of 10,000,000 is allowed for outputting raw dna sequence"<<endl;
            exit(57);
        }
        if (seek_start_pos + start + size > static_cast<uint64_t>(seek_end_pos))
        {
            fprintf(stderr, "Error: dna_t::sequence: Requested region %Zu to %Zu on "
                    "%s exceeds length %Zu\n", 
                    start, end, this->name.c_str(), this->length());
            exit(1);
        }
        
        else {
            source().seekg(seek_start_pos + start, litestream::POS_BEGIN);
            source().read(sequence, size);
        }
    }




    string PrintDNAStrand(dna_strand d){
        string r;
        switch(d){
        case POS: r = "POS"; break;
        case NEG: r = "NEG"; break;
        case POS_NEG: r = "BOTH"; break;
        case NON_STRANDED: r = "NON_STRANDED"; break;
        }
        return r;
    }


    dna_strand ComplementStrand(dna_strand source){

        switch (source){
        case POS: return NEG; break;
        case NEG: return POS; break;
        case POS_NEG:
        case NON_STRANDED: return source; break;
        default: return source; break;
        }
    }


    CONS Nuc2Cons(NUC n){ return n == cis::z ? cis::mNULL : (CONS)(1<<n); }
    NUC Cons2Nuc(CONS c){ 
        if (c == cis::mNULL){ return cis::z; } 
        else {
            int i;
            for (i=0; i < 16; i++) if (c>>i == 1) { break; }
            return (NUC)i;
        }
    }

    void ReverseComplementUppercase(char const* seq, size_t size, char * revcomp_upper)
    {
        if (seq != revcomp_upper)
        {
            strncpy(revcomp_upper, seq, size);
        }
        for (size_t r = 0; r != size; ++r)
        {
            revcomp_upper[r] = base_to_uppercase_complement[static_cast<int>(seq[r])];
        }
        std::reverse(revcomp_upper, revcomp_upper + size);
    }

    void ReverseComplement(char const* seq, size_t size, char * revcomp)
    {
        if (seq != revcomp)
        {
            strncpy(revcomp, seq, size);
        }
        for (size_t r = 0; r != size; ++r)
        {
            revcomp[r] = base_to_complement[static_cast<int>(seq[r])];
        }
        std::reverse(revcomp, revcomp + size);
    }

} //namespace cis
