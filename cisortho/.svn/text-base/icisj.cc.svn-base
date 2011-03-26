#include <sstream>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <utility>
#include <sstream>
#include <getopt.h>

#include "string_tools.h"
#include "search.h"
#include "pattern.h"
#include "pwm.h"
#include "dna.h"
#include "hit.h"
#include "index_trie_scan.h"
#include "cluster.h"
#include "parse_example.h"
#include "input_grammar.hpp"
#include "parse_funcs.hpp"

#include "boost/spirit/tree/parse_tree.hpp"
#include "boost/spirit/tree/tree_to_xml.hpp"
#include "boost/spirit/utility/confix.hpp"


float MarkovLogProb(NUC * pfx, int sz)
{

    static float logfreqs[] = {
        MIN_SCORE                     ,
        log(0.33                     ),
        log(       0.17              ),
        log(0.33 + 0.17              ),
        log(              0.17       ),
        log(0.33 +        0.17       ),
        log(       0.17 + 0.17       ),
        log(0.33 + 0.17 + 0.17       ),
        log(                     0.33),
        log(0.33 +               0.33),
        log(       0.17 +        0.33),
        log(0.33 + 0.17 +        0.33),
        log(              0.17 + 0.33),
        log(0.33 +        0.17 + 0.33),
        log(       0.17 + 0.17 + 0.33),
        log(0.33 + 0.17 + 0.17 + 0.33)
    };

    float frac = 0.0;
    for (int i=0; i < sz; ++i) 
    {
        frac += logfreqs[static_cast<int>(pfx[i])];
    }

    return frac;
}


//here we need the log-sum-exp trick...
//given t1 = log(f1), t2 = log(f2), updates t1 := log(f1 + f2);
void MarkovLogSum(void *accu, NUC * pfx, int sz)
{
    float & t1 = *(float *)accu; //t is log(a), log(b) = MarkovProb(), we want log(a + b)
    float t2 = MarkovLogProb(pfx, sz);
    float m = std::max(t1, t2);
    float tp1 = t1 - m, tp2 = t2 - m;
    float s = exp(tp1) + exp(tp2);
    t1 = log(s) + m;
}


void Increment(void *accumulator, NUC * prefix, int size)
{
    int & counter = * static_cast<int *>(accumulator);
    counter++;
}


struct score_below : public std::unary_function<LOCUS, bool> {
public:
    int thresh;
    score_below(int t) : thresh(t) {}
    bool operator()(LOCUS const& l)
    {
        return l.first->score(l.second.sequence) < thresh;
    }
};



char *input = NULL;
char *hits_file = NULL;
char *cluster_file = NULL;

long lTrieBuffer = 100000;
long lGPosBuffer = 100000;

//should be longer than longest possible chromosome,
//but shorter than roughly the maximum memory to hold
//1/100 that number of hits
const int g_max_dna_search_length = 1000000000;


char const* usage_message =
 "Usage: icisj [options]\n"
 "Options (short or long may be used)\n"
 "\n"
 "-i  --input                input file or stream\n"
 "-t  --trie-buffer-suze     size of buffer for buffered reading of trie file\n"
 "-g  --gpos-buffer-suze     size of buffer for buffered reading of gpos file\n"
 "-h  --hits-file-output     hits output file\n"
 "-c  --cluster-file-output  cluster output file\n";


static struct option long_options[] = {
    {"input", required_argument, 0, 'i'},
    {"trie-buffer-size", required_argument, 0, 't'},
    {"gpos-buffer-size", required_argument, 0, 'g'},
    {"hits-file-output", required_argument, 0, 'h'},
    {"cluster-file-output", required_argument, 0, 'c'},
    {0, 0, 0, 0}
};


int main(int argc, char ** argv)
{


    int c;
    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "i:t:g:h:c:",
                        long_options, &option_index);
    
        if (c == -1) 
        {
            printf("%s\n", usage_message);
            return 0;
            break;

            switch (c)
            {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                printf (" with arg %s", optarg);
                printf ("\n");
                break;
            
            case 'i': 
                input = new char[strlen(optarg) + 1];
                strcpy(input, optarg);
                break;
            case 't': 
                lTrieBuffer = atol(optarg);
                break;
            case 'g': 
                lGPosBuffer = atol(optarg);
                break;
            case 'h':
                hits_file = new char[strlen(optarg) + 1];
                strcpy(hits_file, optarg);
                break;
            case 'c':
                cluster_file = new char[strlen(optarg) + 1];
                strcpy(cluster_file, optarg);
                break;
    
            case '?':
                /* getopt_long already printed an error message. */
                break;
     
            default:
                abort ();
            }
        }
    }

    if (input == NULL)
    {
        fprintf(stdout, "You must provide the name of a properly formatted input file\n");
        fprintf(stdout, "For example:\n\n\n");
        fprintf(stdout, "%s", INPUT_EXAMPLE);
        return 0;
    }

    litestream::postype GPosBuffer = lGPosBuffer;
    litestream::postype TrieBuffer = lTrieBuffer;

    std::ifstream input_stream(input);
    if (! input_stream)
    {
        cerr<<"Couldn't open input "<<input<<endl;
        exit(86);
    }

    input_stream.unsetf(std::ios::skipws); //  Turn off white space skipping on the stream
  
    std::vector<char> input_contents;
    std::copy(std::istream_iterator<char>(input_stream),
              std::istream_iterator<char>(),
              std::back_inserter(input_contents));

    input_stream.close();

    //  ensure that we have enough trailing newlines to eliminate
    //  the need to check for end of file in the grammar.
    input_contents.push_back('\n');
    input_contents.push_back('\n');
  
    std::vector<char>::const_iterator begin_of_input = input_contents.begin();
    std::vector<char>::const_iterator end_of_input = input_contents.end();

    parsed Parsed;
    cis_grammar<parsed> Parsing(Parsed);

    tree_parse_info<std::vector<char>::const_iterator> parse_tree = 
    boost::spirit::ast_parse(begin_of_input, end_of_input, Parsing, 
                             (space_p | comment_p("#")));


    FILE * cluster_stream = fopen(cluster_file, "w");
    if (cluster_stream == NULL)
    {
        fprintf(stderr, "Couldn't open output clusters file %s", cluster_file);
        exit(31);
    }

    FILE * hits_stream = fopen(hits_file, "w");
    if (hits_stream == NULL)
    {
        fprintf(stderr, "Couldn't open output hits file %s", hits_file);
        exit(31);
    }


    //tree_to_xml(std::cout, parse_tree.trees, "");

    if (! parse_tree.full)
    {
        fprintf(stderr, "Couldn't parse input.");
        exit(30);
    }

    set<string> included_species;
    std::vector<std::string>::iterator species_it;
    for (species_it = Parsed.species.begin();
         species_it != Parsed.species.end(); ++species_it)
    {
        included_species.insert(*species_it);
    }

    //Reverse complement all patterns
    std::vector<cis::pattern *> patterns_rc;
    for (size_t i=0; i != Parsed.patterns.size(); ++i)
    {
        if (Parsed.patterns[i]->kind == "matrix" &&
            (Parsed.patterns[i]->label.strand == cis::POS ||
             Parsed.patterns[i]->label.strand == cis::NEG))
        {
            patterns_rc.push_back(cis::MakeRCPattern(Parsed.patterns[i]));
        }
    }

    Parsed.prepend_paths();


    std::vector<string> pattern_names;
    std::vector<cis::pattern *>::iterator pat_it;
    for (pat_it = Parsed.patterns.begin();
         pat_it != Parsed.patterns.end();
         ++pat_it)
    {
        pattern_names.push_back((*pat_it)->label.id);
    }


    int max_patlength = 0;
    for (pat_it = Parsed.patterns.begin();
         pat_it != Parsed.patterns.end();
         ++pat_it)
    {
        max_patlength = std::max(max_patlength, (*pat_it)->max_length);
    }

    //parse the species file to get a list of species names that correspond to the
    std::string dna_index_path(Parsed.dnaindexfile);

    std::ifstream species_stream(dna_index_path.c_str());
    if (! species_stream)
    {
        cerr<<"Couldn't open species list file "<<Parsed.dnaindexfile<<endl;
        exit(50);
    }

    cis::dna_collection dnac;

    species_stream.seekg(0, std::ios::end);
    int endpos = species_stream.tellg();
    species_stream.seekg(0, std::ios::beg);

    while ((int)species_stream.tellg() != endpos) 
    dnac.insert(new cis::dna_t(species_stream, Parsed.data_directory));
    species_stream.close();

    dnac.calc_offsets();
    dnac.open_dnas();

    std::cout<<"Read "<<dnac.size()<<" pieces of DNA, "
             <<dnac.num_bases()<<" total bases."<<endl;

    std::ifstream trie_base;
    litestream trie_stream(trie_base, TrieBuffer);

    trie_stream.open(Parsed.treefile.c_str());
    if (! trie_stream.good()) { 
        fprintf(stderr, "Couldn't open index trie file %s\n", Parsed.treefile.c_str());
        exit(50);
    }
  
    std::ifstream gpos_base;
    litestream gpos_stream(gpos_base, GPosBuffer);

    gpos_stream.open(Parsed.posfile.c_str());

    if (! gpos_stream.good()) { 
        fprintf(stderr, "Couldn't open index gpos file %s\n", Parsed.posfile.c_str());
        exit(50);
    }
  
  
    index_trie_scan<nuc_trie<hit_info> > 
    scantrie(dnac, 
             std::string(Parsed.treefile), 
             std::string(Parsed.posfile), 
             max_patlength);
  
  
    std::cout<<"Parsed index trie info, with "<<scantrie.dnac.size()
             <<" pieces of DNA, "<<scantrie.dnac.num_bases()<<" total bases, "
             <<TotalSize(scantrie, trie_stream)<<" total suffixes."<<endl;
  
    int GenomeSize = TotalSize(scantrie, trie_stream);
    float logGenomeSize = log((float)GenomeSize);
  
    map<cis::pattern *, LOCI> loci;
  
    //now add the reverse complemented patterns.
    Parsed.patterns.insert(Parsed.patterns.end(),
                           patterns_rc.begin(),
                           patterns_rc.end());
  
    //display the matrices of the parsed patterns
    for (size_t p = 0; p != Parsed.patterns.size(); ++p)
    {
        if (Parsed.patterns[p]->kind == "matrix")
        {
            printf("Matrix for %s %s\n", Parsed.patterns[p]->label.id.c_str(),
                   cis::PrintDNAStrand(Parsed.patterns[p]->label.strand).c_str());
            PrintMatrix(* static_cast<IMATRIX const*>(Parsed.patterns[p]->ptr2));
            PrintPWMMatrix(* static_cast<pwm_trie::MAT const*>(Parsed.patterns[p]->ptr));
        }
    }



    //set score thresholds
    for (pat_it = Parsed.patterns.begin();
         pat_it != Parsed.patterns.end();
         ++pat_it)
    {
        cis::pattern * pat = (*pat_it);

        if (pat->kind == "matrix")
        {

            std::pair<int, int> threshold = std::make_pair(0,0);
            int outer_max_matrix_score = 0;
            int max_matrix_score = 0;
            int min_matrix_score = 0;
            int num_estimated_hits = 0;
            float logfreq = -1000000.0;
            int pattern_count = 0;
            int num_hits_to_find = pat->total_hits;
            pwm_trie * ptrie = static_cast<pwm_trie *>(pat->t);

            while (num_hits_to_find > 0 && min_matrix_score >= pat->min_score)
            {

                if (min_matrix_score != pat->min_score)
                {

                    ptrie->set_scores(min_matrix_score, min_matrix_score);
                    ptrie->set_curscore(0);
                    pattern_count = 0;
                    ptrie->apply_to_leaves(&Increment, &pattern_count);
                    if (! pattern_count) {
                        --min_matrix_score;
                        continue;
                    }
          
                    ptrie->set_scores(min_matrix_score, max_matrix_score);
                    ptrie->apply_to_leaves(&MarkovLogSum, &logfreq);
                    num_estimated_hits = 
                    static_cast<int>(std::floor(exp(logfreq + logGenomeSize)));
          
                    printf("For %s %s, estimated %i hits of %i desired at scores [%i, %i], log val = (%f)\n", 
                           pat->label.id.c_str(), cis::PrintDNAStrand(pat->label.strand).c_str(),
                           num_estimated_hits, num_hits_to_find, 
                           min_matrix_score, max_matrix_score,
                           logfreq + logGenomeSize);
          
                    //if the number of estimated hits that would be found by this
                    //score range is less than the number needed, lower the score
                    //band (min,max) by one and continue.
                    if (num_estimated_hits < num_hits_to_find)
                    {
                        --min_matrix_score;
                        max_matrix_score = min_matrix_score;
                        continue;
                    }
                }

                //if we get here it means we have met the desired number of hits
                printf("Searching in score range [%i, %i]\n",
                       min_matrix_score, outer_max_matrix_score);

                ptrie->set_scores(min_matrix_score, outer_max_matrix_score);
                ptrie->set_curscore(0);

                //         trie_stream.seekg(0, litestream::POS_BEGIN);
                //         node_info node = readNode(trie_stream);

                LOCI loci_tmp = 
                Search(scantrie, pat->t, trie_stream, gpos_stream, pat->max_length);
        
                loci[pat].insert(loci[pat].end(), loci_tmp.begin(), loci_tmp.end());
                FindThresholdScore(loci_tmp, num_hits_to_find, &threshold);
                num_hits_to_find -= threshold.second;
                num_estimated_hits = 0;
                logfreq = -1000000.0;
                --min_matrix_score;
                outer_max_matrix_score = min_matrix_score;
            }

            score_below pred(threshold.first);
            LOCI & ploc = loci[pat];
            std::cout<<CountHits(ploc)<<" for "
                     <<pat->label.id<<"("<<pat->cluster_group<<") before remove_if"<<endl;
            ploc.erase(remove_if(ploc.begin(), ploc.end(), pred), ploc.end());
			
        }

        else loci[pat] = Search(scantrie, pat->t,
                                trie_stream, gpos_stream, pat->max_length);

        std::cout<<CountHits(loci[pat])<<" for "
                 <<pat->label.id<<"("<<pat->cluster_group<<")"<<endl;

    }

    //dna_t const* dna;
    cis::REG_MAP hits;

    int cluster_num = 0;

    std::vector<std::string> cluster_names;
    std::vector<cis::cluster *>::iterator cluster_iter;
    for (cluster_iter = Parsed.clusters.begin();
         cluster_iter != Parsed.clusters.end(); ++cluster_iter)
    {
        cluster_names.push_back((*cluster_iter)->name);
    }
    
    std::string ctlist = join_token(cluster_names, ",");

    std::string patlist = join_token(pattern_names, ",");
    std::string splist = join_token(Parsed.species, ",");

    //since it may be memory prohibitive to make all hits at once,
    //make them for a fraction of dna pieces and then filter them by
    //clustering...  traverse the pieces of dna in groups, advancing
    cis::dna_collection::const_iterator start_dna = dnac.begin();
    cis::dna_collection::const_iterator end_dna = dnac.begin();
    cis::dna_collection::const_iterator dna_iter;


    while (end_dna != dnac.end())
    {
        while (end_dna != dnac.end() &&
               dnac.TotalLength(start_dna, end_dna) < g_max_dna_search_length)
        {
            end_dna++; 
        }


        // 		if (included_species.find(dna->species()) == included_species.end()) continue;

        for (pat_it = Parsed.patterns.begin();
             pat_it != Parsed.patterns.end();
             ++pat_it)
        {
            cis::pattern * pat = (*pat_it);
            cis::REG_MAP hits_tmp = makeHits(loci[pat], pat, gpos_stream, dnac, 
                                             start_dna, end_dna,
                                             pat->max_length);

            for (cis::REG_MAP::const_iterator rit = hits_tmp.begin();
                 rit != hits_tmp.end(); ++rit)
            {
                hits[(*rit).first].insert(rit->second.begin(), rit->second.end());
            }

        }

        for (cis::REG_MAP::iterator rit = hits.begin(); rit != hits.end(); ++rit)
        {

            cis::dna_t const* this_dna = rit->first;
            cis::REGIONS_MULTI & this_dna_hits = rit->second;

            
            for (cluster_iter = Parsed.clusters.begin();
                 cluster_iter != Parsed.clusters.end(); ++cluster_iter)
            {
                cis::cluster * clust = (*cluster_iter);

                clust->set_hits_container(&this_dna_hits);
                clust->clear_tags(this_dna_hits.begin(), this_dna_hits.end());
                clust->assign_clusters();

                if (Parsed.merge_overlap) cis::MergeOverlappingClusters(this_dna_hits);
                cluster_num = cis::ReduceClusterNumbers(this_dna_hits, cluster_num);
                PrintHitClusters(this_dna_hits, cluster_stream);
            }

            //filter clusters if there are any clustering criteria 
            if (! Parsed.clusters.empty())
            for (cis::RIT rit = this_dna_hits.begin(); rit != this_dna_hits.end(); rit++)
            {
                cis::hit const* h = static_cast<cis::hit const*>(*rit);
                if (h->is_cluster == 0) { delete h; this_dna_hits.erase(rit); }
            }

            //now append the output to the output file 
            cis::PrintHits(this_dna_hits, *this_dna, hits_stream);

            for (cis::RIT rit = this_dna_hits.begin(); rit != this_dna_hits.end(); rit++)
            {
                cis::hit const* h = static_cast<cis::hit const*>(*rit);
                delete h;
            }
            this_dna_hits.clear();
        }
        hits.clear();
    }

    fclose(cluster_stream);
    fclose(hits_stream);

    for (size_t p = 0; p != Parsed.patterns.size(); ++p) 
    {
        delete Parsed.patterns[p];
    }
		
    return 0;
}
