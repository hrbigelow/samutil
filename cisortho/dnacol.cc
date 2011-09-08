#include <sstream>
#include <cassert>
#include <climits>
#include <cstdio>

#include "dnacol.h"
#include "dna.h"

typedef cis::dna_collection DC;

DC::dna_collection(DNAS const& d) : DNAS(d) {}


void DC::open_dnas(){
    for (cis::DNAS::iterator dna_it = (*this).begin(); 
         dna_it != (*this).end(); ++dna_it)
        (*dna_it)->source().safe_open();
}
    
void DC::calc_offsets(){

    int64_t gp = 0;
    for (cis::DNAS::iterator dna_it = (*this).begin(); 
         dna_it != (*this).end(); ++dna_it){
        offsets[(*dna_it)] = gp;
        offsets_inv[gp] = (*dna_it);
        gp += (*dna_it)->length();
    }
    total_length_ = gp;
}


void DC::make_ready(char const* dna_index_file, std::string dna_directory)
{
    std::ifstream dna_index_stream(dna_index_file);
    if (! dna_index_stream)
    {
        fprintf(stderr, "dna_collection::make_ready: "
                "dna index file %s couldn't be opened for reading\n", 
                dna_index_file);
        exit(1);
    }

    dna_index_stream.seekg(0, std::ios::end);
    int endpos = dna_index_stream.tellg();
    dna_index_stream.seekg(0, std::ios::beg);

    while ((int)dna_index_stream.tellg() != endpos)
    {
        this->insert(new cis::dna_t(dna_index_stream, dna_directory));
    }
    dna_index_stream.close();

    this->calc_offsets();
    this->open_dnas();
    
}

//return the index of a dna / locus pair
int64_t DC::ToIndex(cis::dna_t const* dna, int64_t locus) const {
    OIT it = offsets.find(dna);
    assert(it != offsets.end());
    return (*it).second + locus;
}


int64_t DC::ToIndex(DNAS::const_iterator dna_iterator, int64_t locus) const {
    if (dna_iterator == this->end()) return total_length_;
    else return offsets.find(*dna_iterator)->second + locus;
}


int64_t DC::TotalLength(DNAS::const_iterator start_dna,
                        DNAS::const_iterator end_dna) const {
    return ToIndex(end_dna, 0) - ToIndex(start_dna, 0);
}


//return the dna / locus pair from an index
std::pair<cis::dna_t const*, int64_t> DC::FromIndex(int64_t index) const {
	
    OIIT it = offsets_inv.upper_bound(index); //first element whose key > k.  we want the one before it
    it--;
    return std::make_pair((*it).second, index - (*it).first);
}		



//scan the fasta file and measure the length of each sequence
void DC::add(string const species, 
             string const& fasta_dir,
             string const& fasta_file){

    std::ostringstream errmsg;

    std::string fasta_path(fasta_dir + std::string("/") + fasta_file);
    FILE *ff = fopen(fasta_path.c_str(), "r");
  
    //std::ifstream ff(fasta_path.c_str());
    if (ff == NULL) {
        errmsg<<"Couldn't open Fasta sequence file (option -s): "<<fasta_path<<endl;
        cerr<<errmsg.str();
        exit(58);
    }

    const unsigned int HEADER = 1000;
    char chunk[HEADER];
    char name[HEADER];

    //int64_t maxline = LLONG_MAX;
    int64_t start, end;
    DNAS pieces;

    while (1) {
    
        if (feof(ff)) break;
    
        //ff.getline(chunk, HEADER);
        //if (ff.eof()) break;
    
        if (1 != fscanf(ff, ">%s", name)){
            //if (1 != sscanf(chunk, ">%s", name)){
            cerr<<"Error reading header line:\n"<<chunk
                <<"\nin sequence file "<<fasta_path
                <<".  Must be, for example:\n>chrI\n"<<endl;
            exit(36);
        }

        while (fgetc(ff) != '\n'){
            if (feof(ff)) {
                fprintf(stderr, "Missing a sequence for id %s\n", name);
                exit(36);
            }
        } // eat up the rest until 
        //fscanf(ff, " %*[^\n]\n"); //skip ahead to the next line

        assert(name[0] != 0);
        start = ftello(ff);
        //start = ff.tellg();
        fscanf(ff, "%*s\n");
        //ff.ignore(maxline, '\n'); //assumes fasta format with no newlines in the sequence
        end = ftello(ff) - 1;
        //int s = ff.tellg();
        //end = ff.tellg();// - 1; //because we don't want the newline counted
        //end--; //do not count the newline
    
        this->insert(new dna_t(fasta_path, species, name, start, end));
    }
    //ff.close();
    fclose(ff);

}


int64_t DC::num_bases() const {

    int64_t sz = 0;
    for (cis::DNAS::const_iterator dna_it = (*this).begin(); 
         dna_it != (*this).end(); ++dna_it) 
        sz += (*dna_it)->length();
    return sz;
}


namespace cis {

    cis::dna_t const* GetDNAByName(DNAS const& dnas, string const& species, string const& sdna){

        dna_t const dna("", species, sdna);
        DNAS::const_iterator dna_i = dnas.find(&dna);
        cis::dna_t const* dnap = NULL;

        if (dna_i == dnas.end()){
            static int warn_count = 0;
            if (warn_count < 10) {
                cerr<<"Warning: annotation is on "
                    <<dna.species()<<" "<<dna.name
                    <<" but this piece doesn't exist in the sequence input.  Ignoring."<<endl;
            }
            warn_count++;
        }
        else dnap = *dna_i;
        return dnap;
    }


} // end namespace cis
