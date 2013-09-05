#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <getopt.h>
#include <cstring>
#include <parallel/algorithm>

#include "dep/tools.h"
#include "md5.h"

int const base_to_index[] =
    {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };


char const* bases = "ACGT";

int nmer_spectrum_usage(size_t ndef, size_t idef, size_t Idef, char const* sdef, size_t tdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "nmer_spectrum [OPTIONS] input.cfa output.spectrum\n\n"
            "Options:\n\n"
            "-n      INT    nmer length to tally [%zu]\n"
            "-i      INT    insert size min [%zu]\n"
            "-I      INT    insert size max [%zu]\n"
            "-s      STRING identification string [%s]\n"
            "-t      INT    number threads to use [%zu]\n"
            "-c      FLAG   if present, only test collisions [absent]\n",
            ndef, idef, Idef, sdef, tdef);
    return 1;
}

struct digest
{
    unsigned char val[16];
    digest(unsigned char * v, size_t bytes)
    {
        std::fill(this->val, this->val + 16, 0);
        memcpy(this->val, v, bytes);
    }
    digest(digest const& d)
    {
        memcpy(this->val, d.val, 16);
    }
    bool operator<(digest const& d) const
    {
        return memcmp(this->val, d.val, 16) < 0;
    }
};

typedef std::map<digest, int> DMAP;
typedef std::map<digest, DMAP> DDMAP;

struct less_ddmap_iter
{
    bool operator()(DDMAP::iterator const& d1,
                    DDMAP::iterator const& d2)
    {
        return (*d1).first < (*d2).first;
    }
};

typedef std::pair<DDMAP::iterator, std::vector<std::pair<size_t, size_t> > > RMAP_VAL;
typedef std::map<DDMAP::iterator, std::vector<std::pair<size_t, size_t> >, less_ddmap_iter> RMAP;
typedef std::basic_string<unsigned char> USTRING;
struct md5_collision
{
    size_t bin_length;
    std::map<digest, unsigned char* > * collisions;
    std::set<USTRING> * collided_bstrings;
    size_t * repeat_count;
    md5_collision(size_t bl, std::map<digest, unsigned char* > * c,
                  std::set<USTRING> * cb,
                  size_t * rc) : 
        bin_length(bl), collisions(c), collided_bstrings(cb), 
        repeat_count(rc) { }
    unsigned char * operator()(unsigned char * bin_nmer,
                               unsigned char * digest);
};


unsigned char * md5_collision::operator()(unsigned char * bin_nmer,
                                          unsigned char * this_digest)
{
    std::map<digest, unsigned char* >::iterator it = 
        (*this->collisions).insert
        (std::make_pair(digest(this_digest, 16UL), 
                        static_cast<unsigned char *>(NULL))).first;

    if ((*it).second == static_cast<unsigned char *>(NULL))
    {
        // insertion succeeded.  we haven't seen this before
        (*it).second = new unsigned char[this->bin_length];
        memcpy((*it).second, bin_nmer, this->bin_length);
    }
    else
    {
        (*this->repeat_count)++;
        if (memcmp((*it).second, bin_nmer, this->bin_length) != 0)
        {
            (*this->collided_bstrings).insert(USTRING((*it).second, this->bin_length));
            (*this->collided_bstrings).insert(USTRING(bin_nmer, this->bin_length));
            // fprintf(stderr, "collision\n");
            // fwrite((*it).first.val, 1, 16, stderr);
            // fprintf(stderr, "\t");
            // fwrite(in_string, 1, this->bin_length, stderr);
            // fprintf(stderr, "\t");
            // fwrite((*it).second, 1, this->bin_length, stderr);
            // fprintf(stderr, "\n");
        }
        else
        {
            // fprintf(stderr, "pass\n");
        }
    }
    return this_digest;
}


struct calc_md5
{
    size_t bin_length; // the number of bytes of binary length

    calc_md5(size_t length) : bin_length(length) { }

    unsigned char * operator()(unsigned char *in_string, 
                               unsigned char * digest);
};


unsigned char * calc_md5::operator()(unsigned char * in_string, 
                                     unsigned char * digest)
{
    MD5_CTX ctx;
    MD5Init(&ctx);
    MD5Update(&ctx, in_string, this->bin_length);
    MD5Final(&ctx);
    memcpy(digest, ctx.digest, 16);
    return digest;
}


struct append_maps
{
    unsigned char ** digests;
    size_t digest_bytes;
    append_maps(unsigned char ** d, size_t db) : digests(d), digest_bytes(db) { }
    void operator()(RMAP_VAL d);
};


// data type is a std::pair<DDMAP::iterator, std::vector<std::pair<size_t, size_t> > >
void append_maps::operator()(RMAP_VAL d)
{
    DDMAP::iterator & dit = d.first;
    std::vector<std::pair<size_t, size_t> > & ranges = d.second;
    std::vector<std::pair<size_t, size_t> >::iterator rit;
    for (rit = ranges.begin(); rit != ranges.end(); ++rit)
    {
        size_t E_min = (*rit).first;
        size_t E_max = (*rit).second;
        for (size_t E = E_min; E != E_max; ++E)
        {
            (*dit).second[digest(this->digests[E], this->digest_bytes)]++;
        }
    }
}

void seq2bseq(char const* seq, size_t n, unsigned char **bseq, size_t bn)
{
    std::fill((*bseq), (*bseq) + bn, 0);
    char const* this_seq = seq;
    for (size_t chr = 0; chr != n; ++chr, ++this_seq)
    {
        (*bseq)[chr>>2] |= base_to_index[static_cast<unsigned int>(*this_seq)]<<((chr % 4) * 2);
    }
}


void bseq2seq(unsigned char const* bseq, char ** seq, size_t n)
{
    std::fill((*seq), (*seq) + n, 0);
    char * this_seq = (*seq);
    for (size_t chr = 0; chr != n; ++chr, ++this_seq)
    {
        size_t ind = bseq[chr / 4]>>((chr % 4) * 2) & 3;
        *this_seq = bases[ind];
    }
}


int main(int argc, char ** argv)
{

    size_t ndef = 50;
    size_t idef = 250;
    size_t Idef = 350;
    char sdef[256] = "none";
    size_t tdef = 1;

    size_t nmer_length = ndef;
    size_t insert_size_min = idef;
    size_t insert_size_max = Idef;
    char * id_string = sdef;
    size_t num_threads = tdef;
    bool just_collisions = false;

    char c;

    while ((c = getopt(argc, argv, "n:i:I:s:t:c")) >= 0)
    {
        switch(c)
        {
        case 'n': nmer_length = static_cast<size_t>(atof(optarg)); break;
        case 'i': insert_size_min = static_cast<size_t>(atof(optarg)); break;
        case 'I': insert_size_max = static_cast<size_t>(atof(optarg)); break;
        case 's': strcpy(id_string, optarg); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'c': just_collisions = true; break;
        default: return nmer_spectrum_usage(ndef, idef, Idef, sdef, tdef); break;
        }
    }

    char const* fasta_file;
    char const* spectrum_file;

    if (argc != optind + 2)
    {
        return nmer_spectrum_usage(ndef, idef, Idef, sdef, tdef);
    }

    fasta_file = argv[optind];
    spectrum_file = argv[optind + 1];

    FILE * fasta_fh = open_or_die(fasta_file, "r", "Fasta file");
    FILE * spectrum_fh = open_or_die(spectrum_file, "w", "Output spectrum file");

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    __gnu_parallel::_Settings psettings;
    psettings.for_each_minimal_n = 20;
    __gnu_parallel::_Settings::set(psettings);
    

    size_t nmer_binlength = (nmer_length / 4) + (nmer_length % 4 == 0 ? 0 : 1) ;
    fprintf(stderr, "nmer_binlength = %zu\n", nmer_binlength);

    char header[1024];

    // allocate
    size_t const max_transcript_length = 120000;

    char * seq = new char[max_transcript_length];
    unsigned char * bseq[4];
    unsigned char term_mask;

    switch(nmer_length % 4)
    {
    // case 0: term_mask = 255; break; //  11111111
    // case 1: term_mask = 192; break; //  11111100
    // case 2: term_mask = 240; break; //  11110000
    // case 3: term_mask = 252; break; //  11000000

    case 0: term_mask = 255; break; //  11111111
    case 1: term_mask = 3; break;   //  00000011
    case 2: term_mask = 15; break;  //  00001111
    case 3: term_mask = 63; break; //   00111111
    }

    DDMAP fragments;
    DDMAP::iterator fit1;
    DMAP::const_iterator fit2;
    RMAP range_map;
    std::vector<RMAP_VAL> range_vec;

    unsigned char * bbuf = new unsigned char[max_transcript_length * nmer_binlength * 2];
    unsigned char * dbuf = new unsigned char[max_transcript_length * 16];

    unsigned char * trunc_bseq[max_transcript_length];
    unsigned char * digests[max_transcript_length];
    unsigned char * bptr = bbuf;
    unsigned char * dptr = dbuf;

    for (size_t b = 0; b != max_transcript_length; ++b, bptr += nmer_binlength * 2)
    {
        trunc_bseq[b] = bptr;
    }
    for (size_t d = 0; d != max_transcript_length; ++d, dptr += 16)
    {
        digests[d] = dptr;
    }

    // map to monitor possible md5 collisions
    size_t digest_bytes = std::min(16UL, nmer_binlength);
    size_t num_distinct_points = 0;
    size_t num_total_points = 0;

    std::map<digest, unsigned char* > collisions;
    std::set<USTRING> collided_bstrings;
    size_t collision_count = 0;
    
    md5_collision collision_op(nmer_binlength, & collisions,
                               & collided_bstrings, 
                               & collision_count);

    while (fscanf(fasta_fh, "%[^\n]\n%s\n", header, seq) == 2)
    {
        size_t L = strlen(seq);
        size_t N_nmer = nmer_length < L ? L - nmer_length : 0;
        size_t S_max = insert_size_min < L ? L - insert_size_min : 0;

        if (N_nmer == 0 || S_max == 0)
        {
            continue;
        }

        // pre 1. convert entire seq to b-string in 4 offsets
        for (size_t reg = 0; reg != 4; ++reg)
        {
            bseq[reg] = new unsigned char[L / 4 + 1];
            seq2bseq(seq + reg, L - reg, & bseq[reg], L / 4 + 1);
        }

        // pre 2. generate list of truncated b-strings.  Must zero out remaining bits
        for (size_t S = 0; S != N_nmer; ++S)
        {
            memcpy(trunc_bseq[S], bseq[S % 4] + (S / 4), nmer_binlength);
            trunc_bseq[S][nmer_binlength - 1] &= (term_mask);
        }
        
        // transform all binary sequences to md5 checksums inplace
        if (nmer_binlength > 16)
        {
            __gnu_parallel::transform(trunc_bseq, trunc_bseq + N_nmer,
                                      digests, digests, calc_md5(nmer_binlength));

            if (just_collisions)
            {
                std::transform(trunc_bseq, trunc_bseq + N_nmer,
                               digests, digests, collision_op);
            }
        }
        else
        {
            for (size_t b = 0; b != N_nmer; ++b)
            {
                memcpy(digests[b], trunc_bseq[b], nmer_binlength);
            }
        }
        
        range_map.clear();

        fprintf(stderr, "%s\n", header);

        if (just_collisions)
        {
            continue;
        }

        for (size_t S = 0; S != S_max; ++S)
        {
            fit1 = fragments.insert(std::make_pair(digest(digests[S], digest_bytes), 
                                                    std::map<digest, int>())).first;

            size_t E_min = S + insert_size_min - nmer_length;
            size_t E_max = std::min(S + insert_size_max - nmer_length, N_nmer);

            range_map[fit1].push_back(std::make_pair(E_min, E_max));
            
        }

        range_vec.clear();
        range_vec.reserve(range_map.size());
        for (RMAP::iterator rit = range_map.begin(); rit != range_map.end(); ++rit)
        {
            range_vec.push_back(*rit);
        }
        // take each iterator, and for each range, increment entry in the map
        
        __gnu_parallel::for_each(range_vec.begin(), range_vec.end(), 
                                 append_maps(digests, digest_bytes));

        for (size_t reg = 0; reg != 4; ++reg)
        {
            delete bseq[reg];
        }
    }

    fclose(fasta_fh);

    //consolidate the spectrum
    size_t spectrum_hist[101];
    size_t trunc;
    std::fill(spectrum_hist, spectrum_hist + 101, 0);

    for (fit1 = fragments.begin(); fit1 != fragments.end(); ++fit1)
    {
        for (fit2 = (*fit1).second.begin(); fit2 != (*fit1).second.end(); ++fit2)
        {
            trunc = std::min((*fit2).second, 100);
            spectrum_hist[trunc]++;
            num_distinct_points++;
            num_total_points += (*fit2).second;
        }
    }

    for (size_t bin = 0; bin != 101; ++bin)
    {
        fprintf(spectrum_fh, "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n", 
                id_string, nmer_length, insert_size_min, insert_size_max, 
                bin, spectrum_hist[bin], num_distinct_points,
                num_total_points);
    }

    if (just_collisions)
    {
        fprintf(stderr, "repeat_count: %zu\ncollision_count: %zu\n",
                *collision_op.repeat_count, 
                (*collision_op.collided_bstrings).size());
    }

    fclose(spectrum_fh);
    delete seq;
    delete bbuf;
    delete dbuf;

    return 0;

}
