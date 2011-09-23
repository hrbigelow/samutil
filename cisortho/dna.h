#ifndef _DNA
#define _DNA

#include <string>
#include <iostream>
#include <utility>
#include <fstream>
#include <cstdlib>
#include <stdint.h>

#include "memory.h"
#include "dna_types.h"
#include "litestream.h"

using std::string;
using std::ostream;
using std::istream;
using std::cerr;
using std::endl;
using std::pair;




namespace cis {
    class dna_t {
    public:
        static string iupac;
        static string iupac_comp;

        static char base_to_complement[256];
        static char base_to_uppercase_complement[256];


        class named_stream : public litestream {
        public:
            std::string p;
            std::ifstream base;

        named_stream(std::string const& p) : litestream(base), p(p) { }

            //danger--does not copy!  you must open the stream *after* it's stored...
        named_stream(named_stream const& n) : litestream(base), p(n.p) { }
            void safe_open(){
                if (is_open()) return;
                open(p.c_str());
                if (! this->good()){
                    cerr<<"Couldn't open dna stream "<<p<<endl;
                    exit(50);
                }
            }
            ~named_stream(){ if (is_open()) close(); }
      
            friend bool operator<(named_stream const &s1, named_stream const& s2){ return s1.p < s2.p; }
        };
    
        memory::cstore<string> species;
        memory::store<named_stream> source;
	

        string name;
        int64_t length() const { return seek_end_pos - seek_start_pos; };

        int64_t seek_start_pos; //position in the source file of the start position
        int64_t seek_end_pos;

        //std::string const source_path; //source file path containing the dna

        dna_t(std::string const& path, 
              string const& species, 
              string const& _name, 
              int64_t const _seek_start_pos = 0, 
              int64_t const _seek_end_pos = 0);

        dna_t(istream &, std::string const&);

        dna_t(std::string const& _species,
              std::string const& _contig_name,
              int64_t const _seek_start_pos = 0, 
              int64_t const _seek_end_pos = 0);
            
        //dna_t(string const& s) : species(s), name(""), seek_start_pos(0), seek_end_pos(0) {}
        //dna_t(dna_t const&);
        dna_t(dna_t const& d) : name(d.name), 
                                seek_start_pos(d.seek_start_pos), 
                                seek_end_pos(d.seek_end_pos)
        { 
            species(d.species());
        }

        //dna_t &operator=(dna_t const&) { assert(false); string s; return dna_t(s); }
        //dna_t(string="",int=0, );
        friend ostream& operator<<(ostream&,dna_t const&);
        //friend istream& operator>>(istream&,dna_t&);
        friend bool operator==(dna_t const&, dna_t const&);
        friend bool operator!=(dna_t const&, dna_t const&);
        friend bool operator<(dna_t const&, dna_t const&);

        string sequence(uint64_t, uint64_t) const;
        void sequence(uint64_t, uint64_t, char * sequence) const;
    };


    CONS Nuc2Cons(NUC n);
    NUC Cons2Nuc(CONS c);

 
    string PrintDNAStrand(dna_strand);
 

    dna_strand ComplementStrand(dna_strand source);


    NUC CompNUC(NUC);
    char CompChar(char);
    dnastring RevCompDNA (dnastring const& seq);
    string RevCompStr(string const& s);
    char Nuc2char(NUC d);
    NUC char2NUC(unsigned char c);
    //CONS TranslateCons(char c);
    void TranslateSeq(NUC *const iseq, char const*const seq, int const size);
    string FromDNAString(dnastring const& d);
    string FromDNAString(NUC const*const sequence, int size);
    //dnastring ToDNAString(string const& s);
    dnastring ToDNAString(string const s);
    dnastring ToDNAString(char const*const, int);
  
    void ReverseComplementUppercase(char const* seq, size_t size,
                                    char * revcomp_upper);

    void ReverseComplement(char const* seq, size_t size, char * revcomp);
  
} // end namespace cis
  
#endif //_DNA
  
