//functions for calculating INTRON, EXON and other relationships from a region's genomic context

#include "region_association.h"

#include <sstream>
#include <vector>
#include <cassert>

#include "region.h"
#include "nested.h"
#include "misc.h"
#include "enum.h"

/* 
   
   overall:  
   1. build an NCList of all regions on a piece of DNA
   2. for each query region:
   a.  pad it with max_assoc_distance
   b.  query NCList, getting a set of overlapping regions
   c.  build a region association, recording all required output fields.
   DNA annotation regions are assumed to be non-overlapping.
   
   output:
   association_id query_region_id distance same_strand 
   number_intervening target_region_id

   find and print all query to target region associations based on
   the association distance and maximum number of regions between.
   use <first_association_id> as the start value of the counter for
   the association records. 
   returns the last association id used + 1
   prints:
   association_id, query_region_id, region_distance, 
   same_strand, num_regions_between, target_region_id

*/
int cis::PrintRegionAssociations(cis::REG_MAP const& query_regions_on_dna,
                                 cis::REG_MAP const& target_regions_on_dna,
                                 int max_assoc_distance, /* negative value has the effect of requiring a minimum overlap */
                                 int min_assoc_distance, /* negative value has the effect of setting a maximum overlap */
                                 int max_regions_between,
                                 int first_association_id,
                                 FILE * outstream)
{
  
    cis::REGIONS_MULTI overlapping;
    int region_distance, region_overlap;
    int num_regions_between;
    bool same_strand;
    int association_id = first_association_id;
  
    //left and right distance ordering are for keeping track of the
    //number of intervening regions to the left and right
    std::set<int> left_distance_ordering;
    std::set<int> right_distance_ordering;
  
    for (cis::REG_MAP::const_iterator rit = query_regions_on_dna.begin(); 
         rit != query_regions_on_dna.end(); ++rit)
    {
        dna_t const* dna = (*rit).first;
        assert(dna->species() != "");
        cis::REGIONS_MULTI const& query_regions = (*rit).second;
        cis::REGIONS_MULTI const& target_regions = 
            (*target_regions_on_dna.find(dna)).second;
    
        region_tree * region_tree = BuildTree(target_regions);

        //print_tree(region_tree, 0);

        for (cis::REGIONS_MULTI::const_iterator qr_iter = query_regions.begin(); 
             qr_iter != query_regions.end(); ++qr_iter)
        {
            region const* query_region = *qr_iter;
        
            region const* & qr = query_region;
            region * query_region_padded = 
                MakeRegion(dna, qr->start, qr->end,
                           PrintDNAStrand(qr->strand), qr->name, 0, 
                           std::max(0, max_assoc_distance),
                           std::max(0, max_assoc_distance), 
                           qr->id, 0);

            overlapping = IntervalOverlap(*region_tree, *query_region_padded);


            //get the ordering of region distances of regions to the left and right
            left_distance_ordering.clear();
            right_distance_ordering.clear();

            for (RIT_M overlap_iter = overlapping.begin();
                 overlap_iter != overlapping.end(); ++overlap_iter)
            {
                REG_PC touching_region = *overlap_iter;
                region_distance = RegionDistance(*query_region, *touching_region);
                            
                if (query_region->end > touching_region->end)
                {
                    right_distance_ordering.insert(region_distance);
                }

                if (query_region->start < touching_region->start)
                {
                    left_distance_ordering.insert(region_distance);
                }
            }

            //do it all again, this time outputting the data
            for (cis::REGIONS_MULTI::const_iterator overlap_iter = overlapping.begin();
                 overlap_iter != overlapping.end(); ++overlap_iter)
            {
                region const* touching_region = *overlap_iter;
                region_distance = RegionDistance(*query_region, *touching_region);
                region_overlap = RegionOverlap(*query_region, *touching_region);
                if (region_distance < min_assoc_distance ||
                    region_distance > max_assoc_distance) 
                {
                    continue;
                }

                same_strand = query_region->strand == touching_region->strand;

                if (query_region->end > touching_region->end)
                {
                    //query bounds target on the right
                    num_regions_between = 
                        std::distance(right_distance_ordering.begin(),
                                      right_distance_ordering.find(region_distance));
                }

                else if (query_region->start < touching_region->start)
                {
                    //query bounds target on the left
                    num_regions_between = 
                        std::distance(left_distance_ordering.begin(),
                                      left_distance_ordering.find(region_distance));
                }

                else 
                {
                    //target contains query
                    num_regions_between = 0;
                }

                if (num_regions_between <= max_regions_between)
                {
                    fprintf(outstream, "%i\t%i\t%i\t%i\t%i\t%s\t%i\n",
                            association_id++, query_region->id, touching_region->id,
                            region_distance, region_overlap, 
                            PrintStrandBool(same_strand).c_str(), 
                            num_regions_between);
                }
            }
            delete query_region_padded;
        }
        delete region_tree;
    }
    return association_id;
}
