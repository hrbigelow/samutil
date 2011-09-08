#include <iostream>
#include <numeric>
#include <string>
#include <cstdio>
#include <cassert>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "region.h"
#include "hit.h"
#include "misc.h"
#include "enum.h"
#include "dnacol.h"


using namespace std;
using namespace cis;
using namespace cis_enum;

namespace cis {

    int StrandToEnum(char s)
    {

        int es;

        switch (s)
        {
        case '+': case '1': case 'p': case 'P': es = POS; break;
        case '-': case '0': case 'n': case 'N': es = NEG; break;
        default:
            cerr<<"Warning: don't understand strand type "
                <<s<<".  Assigning 'NON_STRANDED'"<<endl;
            es = NON_STRANDED; break;
        }
        return es;

    }

    //assumes two regions are on the same dna.
    //returns the displacement necessary for the to be
    //adjacent, either negative or positive.
    int64_t RegionDistance(REG_CR r1, REG_CR r2)
    {

        assert(r1.dna == r2.dna);
        int d1 = r1.start - r2.end;
        int d2 = r2.start - r1.end;
        return (abs(d1) < abs(d2)) ? d1 : d2;

    }


    int64_t RegionOverlap(REG_CR r1, REG_CR r2)
    {

        assert(r1.dna == r2.dna);

        int s = max(r1.start, r2.start);
        int e = min(r1.end, r2.end);
        int l = e - s;

        return max(0, l);
    }


    int region::region_id = 0;

    const int region::NONE = 0;
    const int region::LEFT_MOST = 1<<0;
    const int region::RIGHT_MOST = 1<<1;

    region::region(dna_t const& d, string n, int64_t s, 
                   int64_t e, dna_strand ds, int i, int g,
                   int gr, void const* _payload) :
        dna(d), group_id(g), start(s), end(e), strand(ds), 
        name(n), group_relation_(gr), payload(_payload)
    {
        id = (i == -1) ? region_id++ : i;
        if (start > dna.length() || end > dna.length())
        {
            cerr<<"Warning: region "<<dna.species()<<": "<<dna.name<<" "<<dna.length()
                <<" long, region from "<<start<<" to "<<end<<endl;
        }
        if (start > end)
        {
            cerr<<"Error: Trying to create a backwards region."<<endl;
            exit(75);
        }

    }
	
	
    istream& operator>>(istream&i, region& r)
    {
        char strand;
        i>>r.start>>r.end>>strand;
        switch (strand)
        {
        case '+': case '1': case 'p': case 'P': r.strand = POS; break;
        case '-': case '0': case 'n': case 'N': r.strand = NEG; break;
        default:
            cerr<<"Warning: don't understand strand type "
                <<strand<<".  Assigning 'NON_STRANDED'"<<endl;
            r.strand = NON_STRANDED; break;
        }
	
        if (r.start >= r.end)
        {
            cerr<<"Error: region has start="<<r.start<<" >= end="<<r.end<<endl;
            exit(1);
        }
        return i;
    }

    //structural comparison
    bool operator==(REG_CR r1, REG_CR r2)
    {
        return 
        // 		(r1.dna == r2.dna) &&
        (r1.start == r2.start) &&
        (r1.end == r2.end) &&
        (r1.strand == r2.strand);
    }



    bool less_region_ptr::operator()(REG_P r1, REG_P r2) const 
    { 
        return *r1 < *r2; 
    }

    bool region::contains(REG_CR r) const 
    {
        return 
        // 		dna == r.dna &&
        start <= r.start &&
        end >= r.end;
    }


    bool region::overlaps(REG_CR r) const
    {
        assert(dna == r.dna);
        return end > r.start && r.end > start;
    }

    bool operator<(REG_CR r1, REG_CR r2)
    {
        return 
        // 		r1.dna < r2.dna ||
        // 		(r1.dna == r2.dna &&
        (r1.start < r2.start ||
         (r1.start == r2.start &&
          (r1.end < r2.end ||
           (r1.end == r2.end &&
            (r1.strand < r2.strand)))));
    }

    //   string lower_case_regions(string const& seq, REG_PC r, 
    //                             cis::region_tree * regs)
    //{

    // 	REGIONS_MULTI overlaps = cis::IntervalOverlap(*regs, *r);

    // 	int start, end, size;
    // 	string seq_mask(seq);

    //  REGIONS_MULTI::iterator overlaps_iter;
    //  for (overlaps_iter = overlaps.begin();
    //  overlaps_iter != overlaps.end(); ++overlaps_iter)
    //{
    //  REG_PC o = (*overlaps_iter);
    //       start = max((int)o->start - (int)r->start, 0);
    //       end = min((int)r->end, (int)o->end - (int)r->start);
    //       size = end - start;
    //       string sub = seq_mask.substr(start, size);
    //       boost::to_lower(sub);
    //       seq_mask.replace(start, size, sub);
    // 	}
    // 	cout<<seq<<endl;
    // 	cout<<seq_mask<<endl;
    // 	return seq_mask;
    //   }


    //   offset_type relative_pos(REG_CR sub, REG_CR main)
    //{

    // 	assert(sub.dna == main.dna);
    // 	offset_type r;

    // 	if (main.strand == POS)
    //{
    //       if (sub.end <= main.start) r = UPSTREAM;
    //       else if (sub.start >= main.end) r = DOWNSTREAM;
    //       else r = IN;
    // 	}
    // 	else {
    //       if (sub.end <= main.start) r = DOWNSTREAM;
    //       else if (sub.start >= main.end) r = UPSTREAM;
    //       else r = IN;
    // 	}

    // 	return r;
    //   }


    //a kludge for the lack of random access iterators...
    RIT_M_R min_iter(REGIONS_MULTI * c, RIT_M_R a, RIT_M_R b)
    {

        if (a == c->end()) return b;
        else if (b == c->end()) return a;
        else {
            assert((*a)->dna == (*b)->dna);
            return (c->key_comp()(*a, *b)) ? a : b;
        }
    }



    REGIONS_MULTI MergeRegions(REG_MAP r)
    {
        REGIONS_MULTI reg;
        REG_MAP::iterator reg_iter;
        for (reg_iter = r.begin(); reg_iter != r.end(); ++reg_iter)
        {
            REG_MAP::value_type & d = (*reg_iter);
            reg.insert(d.second.begin(), d.second.end());
        }

        return reg;
    }


    REG_MAP SplitRegions(REGIONS_MULTI const& regions)
    {
        REG_MAP contig_to_regions;
        REGIONS_MULTI::const_iterator reg_iter;
        for (reg_iter = regions.begin(); reg_iter != regions.end(); ++reg_iter)
        {
            dna_t const* contig = &(*reg_iter)->dna;
            contig_to_regions[contig].insert(*reg_iter);
        }
        return contig_to_regions;
    }

    //here we don't depend on the source information to tell us what numbers the exons have...
    REG_P MakeRegion(cis::dna_t const* dna, 
                     int64_t start,
                     int64_t end,
                     string const& strand,
                     string const& name,
                     int start_adjust,
                     int five_prime_flank,
                     int three_prime_flank,
                     int id,
                     int group_id)
    {
	
        dna_strand bstrand = POS;
        ostringstream errmsg;
	
        bstrand = (dna_strand)cis_enum::from_string(strand);
	
        REG_P reg = NULL;
	
        if (start == 0 && start_adjust == -1)
        {
            cerr<<"this file appears to be 0-based, but was specified as ones-based"<<endl;
            exit(852);
        }

        static int num_out_of_range = 0;

        if (start + start_adjust > dna->length() ||
            end > dna->length())
        {
            if (num_out_of_range < 5)
            {
                fprintf(stderr,
                        "Warning: Taking subset of region %"PRId64" -> %"PRId64
                        "on dna %s %s (%"PRId64" bases long)\n",
                        start, end, dna->species().c_str(),
                        dna->name.c_str(), dna->length()
                        );
                num_out_of_range++;
                if (num_out_of_range == 5)
                fprintf(stderr, "Further warnings omitted\n");
            }
        }      

        int64_t start_final = bstrand == POS ?
        max(int64_t(0), (int64_t)(start + start_adjust - five_prime_flank)) :
        max(int64_t(0), (int64_t)(start + start_adjust - three_prime_flank));
        int64_t end_final = bstrand == POS ?
        min((int64_t)dna->length(), (int64_t)(end + three_prime_flank)) :
        min((int64_t)dna->length(), (int64_t)(end + five_prime_flank));
    
        reg = new region(*dna, name, start_final, end_final, bstrand, id, group_id,
                         region::NONE, NULL);
    
        return reg;
    
    }
  


    //load regions from an rdbfile.
    REG_MAP RDBToRegions(DNAS& dnas, char const* rdbfile,
                         bool with_group_field)
    {

        int region_id, group_id = -1;
        int64_t bound1, bound2;
        char species[256], dnaname[256], strand[16], name[256];
        REG_MAP regions;

        FILE * rdbstream = fopen(rdbfile, "r");
        if (rdbstream == NULL)
        {
            fprintf(stderr, "Couldn't open RDB regions file %s\n", rdbfile);
            exit(50);
        }

        int num_fields = with_group_field ? 7 : 6;
        int num_assign = 0;

        while (! feof(rdbstream))
        {

            if (with_group_field) 
            num_assign = 
            std::fscanf(rdbstream, "%i %i %s %s %"PRId64" %"PRId64" %s\n",
                        &region_id, &group_id, species, dnaname, 
                        &bound1, &bound2, strand);
      
            else 
            num_assign = 
            std::fscanf(rdbstream, "%i %s %s %"PRId64" %"PRId64" %s\n",
                        &region_id, species, dnaname, 
                        &bound1, &bound2, strand);
      
            if (num_assign == num_fields)
            { //success reading this line
                std::string sspecies(species);
                std::string sdnaname(dnaname);
                cis::dna_t const* dna = GetDNAByName(dnas, sspecies, sdnaname);

                if (dna == NULL)
                {
                }
                else
                {
                    int64_t start = min(bound1,bound2);
                    int64_t end = max(bound1,bound2);
                    REG_P reg = MakeRegion(dna, start, end, std::string(strand), 
                                           std::string(name),
                                           0, 0, 0, region_id, group_id);
          
                    if (regions.find(dna) == regions.end()) regions[dna] = REGIONS_MULTI();
                    if (reg != NULL) regions[dna].insert(reg);
                }
            } 
            else
            {

                if (! feof(rdbstream))
                {
                    char line[1000];
                    fgets(line, 1000, rdbstream);
                    fprintf(stderr,
                            "Found %i fields (should be %i) in RDB regions file. (Next line is:\n%s\n",
                            num_assign, num_fields, line);
                    exit(35);
                }
            }
      
      
        }
    
        fclose(rdbstream);
    
        return regions;
    
    }      
  
  

    //find the left and rightmost intervals among each with the same group_id
    //and set their group_relation_ mask
    void SetGroupRelation(RIT const& begin, RIT const& end)
    {
        set<int> group_ids;
        int group_id;
        for (RIT rit = begin; rit != end; ++rit)
        {
            group_id = (*rit)->group_id;
            if (group_ids.find(group_id) == group_ids.end())
            {
                group_ids.insert(group_id);
                const_cast<region *>(*rit)->group_relation_ |= region::LEFT_MOST;
            }
        }
        group_ids.clear();
        for (RIT_REV rit = RIT_REV(end); rit != RIT_REV(begin); ++rit)
        {
            group_id = (*rit)->group_id;
            if (group_ids.find(group_id) == group_ids.end())
            {
                group_ids.insert(group_id);
                const_cast<region *>(*rit)->group_relation_ |= region::RIGHT_MOST;
            }
        }      
    }

    //returns a tab-separated representation of the region and its DNA sequence,
    //reverse complemented if region is negative stranded
    void PrintDNARegion(REG_PC dna_region, FILE * outfile)
    {
        std::fprintf(outfile, "%i\t%s\t%s\t%"PRId64"\t%"PRId64"\t%s\n",
                     dna_region->id,
                     dna_region->dna.species().c_str(),
                     dna_region->dna.name.c_str(),
                     dna_region->start, 
                     dna_region->end,
                     PrintDNAStrand(dna_region->strand).c_str());
    }




    //simply prints the region's id and sequence (which is reverse complemented if
    //the region is on the negative strand)
    void PrintDNARegionSimple(REG_PC dna_region, int max_sequence_length, 
                              FILE * outfile)
    {
        if (dna_region->length() > max_sequence_length) return;
        string dna_sequence = 
        dna_region->dna.sequence(dna_region->start, dna_region->end);
        if (dna_region->strand == cis::NEG) dna_sequence = RevCompStr(dna_sequence);
        std::fprintf(outfile, "%i\t%s\n", dna_region->id, dna_sequence.c_str());
    }


  

  

} // namespace cis
