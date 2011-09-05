#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "meta_gene.h"
#include "gtf.h"


#include <cassert>

/*
  Class for organizing and transforming exon coordinates
*/

MetaGene::MetaGene() : gene_extent_tree(NULL), meta_exon_tree(NULL), 
                       bounds_initialized(false) { }


MetaGene::~MetaGene()
{
    if (this->gene_extent_tree != NULL)
    {
        free_payload<GENE_ITERS>(this->gene_extent_tree);
        delete this->gene_extent_tree;
    }
    if (this->meta_exon_tree != NULL)
    {
        delete this->meta_exon_tree;
    }
}



void MetaGene::Initialize(FILE * gtf_fh, size_t pseudo_intron_length)
{
    GTFEntry gtf_entry;
    std::map<std::string, FeatureJumps> exon_bounds;
    std::map<std::string, FeatureJumps>::iterator bounds_iter;

    this->gene_extent_tree = new IntervalTree();
    //IntervalTree * meta_exon_tree = new IntervalTree();

    size_t cumul_height = 0;
    size_t start = 0;
    size_t end = 0;

    //parse GTF file, loading the exon_bounds structure
    //initialize a temporary exon_bounds structure
    //initialize this->transcripts
    while (gtf_entry.get_next_record(gtf_fh))
    {
        if (strcmp(gtf_entry.feature, "exon") != 0)
        {
            continue;
        }

        this->transcripts[gtf_entry.gene_id].insert(gtf_entry.transcript_id);

        bounds_iter = exon_bounds.find(gtf_entry.seqname);
        if (bounds_iter == exon_bounds.end())
        {
            bounds_iter = exon_bounds.insert
                (exon_bounds.begin(),
                 std::make_pair(gtf_entry.seqname, FeatureJumps()));
        }

        FeatureJumps & jumps = (*bounds_iter).second;

        jumps.insert(gtf_entry.start_bound(), true, gtf_entry.gene_id);
        jumps.insert(gtf_entry.end_bound(), false, gtf_entry.gene_id);

    }

    //set of feature names for the current region
    std::set<unique_gene_description> current_genes;
    std::set<unique_gene_description>::const_iterator gene_iter;

    //initialize this->meta_exons
    //initialize this->gene_bounds
    for (bounds_iter = exon_bounds.begin(); 
         bounds_iter != exon_bounds.end(); 
         ++bounds_iter)
    {
        std::string contig = (*bounds_iter).first;
        FeatureJumps & jumps = (*bounds_iter).second;

        if (this->meta_exon.find(contig) == this->meta_exon.end())
        {
            this->meta_exon[contig] = Cigar::CIGAR_VEC();
            end = 0; // this signifies we haven't seen a previous end
                     // on this contig.
        }

        Cigar::CIGAR_VEC & cigar_vec = this->meta_exon[contig];

        //process this contig
        for (jumps.init(); ! jumps.is_last(); jumps.next())
        {
            size_t boundary = jumps.boundary();
            int jump_height = jumps.height();

            assert(cumul_height + jump_height >= 0);

            //record current genes for this contiguous feature
            std::set<std::string>::const_iterator feature;
            std::set<std::string> const& features = jumps.features();
            for (feature = features.begin(); feature != features.end(); ++feature)
            {
                current_genes.insert(unique_gene_description(contig.c_str(), 
                                                             (*feature).c_str(), '+'));
            }

            /* at the start of a region (when height == 0 before the jump)
               1. record the start boundary
               2. optionally add a meta_intron
               3. add a pseudo_spacer
            */
            if (cumul_height == 0)
            {
                //at start of region
                assert(jump_height > 0);
                start = boundary;
                if (end == 0)
                {
                    //we're at the start of the very FIRST region in the contig
                    cigar_vec.push_back(Cigar::Unit(Cigar::Ops[Cigar::D], start));
                }
                else
                {
                    //this isn't the first region
                    assert(end <= start);
                    size_t meta_intron_length = start - end;
                    //bool do_shrink = pseudo_intron_length < meta_intron_length;

                    cigar_vec.push_back(Cigar::Unit(Cigar::Ops[Cigar::I], pseudo_intron_length));
                    cigar_vec.push_back(Cigar::Unit(Cigar::Ops[Cigar::D], meta_intron_length));

                    // if (do_shrink)
                    // {
                    //     cigar_vec.push_back
                    //         (Cigar::Unit(Cigar::D, 
                    //                      meta_intron_length - pseudo_intron_length));
                    // }
                    // else
                    // {
                    //     cigar_vec.push_back
                    //         (Cigar::Unit(Cigar::I, 
                    //                      pseudo_intron_length - meta_intron_length));
                    // }
                }
            }

            /*
              at the end of a region:
              1. record the end boundary
              2. from start and end boundary, create a 'M' state and add to cigar_vec
              3. for all feature_names, update gene_bounds for both start and end.
            */

            if (cumul_height + jump_height == 0)
            {
                end = boundary;
                cigar_vec.push_back(Cigar::Unit(Cigar::Ops[Cigar::M], end - start));
                //this->meta_exon_tree->insert(contig.c_str(), start, end, NULL);

                size_t cigar_index = cigar_vec.size() - 1;
                
                //printf("cigar_index: %Zu\n", cigar_index);

                for (gene_iter = current_genes.begin();
                     gene_iter != current_genes.end();
                     ++gene_iter)
                {
                    unique_gene_description const& gene = (*gene_iter);
                    
                    GENE_INDEX::iterator bounds_iter = this->gene_bounds.find(gene);
                    if (bounds_iter == this->gene_bounds.end())
                    {
                        bounds_iter = this->gene_bounds.insert
                            (this->gene_bounds.begin(),
                             std::make_pair(gene, BoundsInfo(UINT64_MAX, 0, UINT64_MAX)));
                    }
                    
                    BoundsInfo & bounds = (*bounds_iter).second;

                    bounds.start_boundary = std::min(bounds.start_boundary, start);
                    bounds.end_boundary = std::max(bounds.end_boundary, end);
                    bounds.start_cigar_index = std::min(bounds.start_cigar_index, cigar_index);
                }
                current_genes.clear();
            }
            cumul_height += jump_height;
        }
        assert(cumul_height == 0);
    }

    //initialize this->meta_exon_offsets
    std::map<std::string, Cigar::CIGAR_VEC>::const_iterator exon_iter;
    for (exon_iter = this->meta_exon.begin(); 
         exon_iter != this->meta_exon.end();
         ++exon_iter)
    {
        std::string const& contig = (*exon_iter).first;
        Cigar::CIGAR_VEC const& cigar_vec = (*exon_iter).second;

        this->meta_exon_offsets[contig] = Cigar::ComputeOffsets(cigar_vec);
    }

    std::pair<BY_START::iterator, bool> insert_result;

    //for looking up gene information from intervals
    GENE_ITERS * gene_iters;
    GENE_INDEX::iterator gene_bounds_iter;
     
    for (gene_bounds_iter = this->gene_bounds.begin();
         gene_bounds_iter != this->gene_bounds.end();
         ++gene_bounds_iter)
    {
        gene_iters = new GENE_ITERS();
        insert_result =
            this->gene_extent_tree->insert((*gene_bounds_iter).first.contig.c_str(),
                                           (*gene_bounds_iter).second.start_boundary,
                                           (*gene_bounds_iter).second.end_boundary,
                                           gene_iters);
        if (! insert_result.second)
        {
            delete gene_iters;
            gene_iters = static_cast<GENE_ITERS *>((*insert_result.first).second->payload);
        }
        (*gene_iters).push_back(gene_bounds_iter);
        
    }
    assert(check_tree(this->gene_extent_tree));
}



size_t MetaGene::transcript_number(GTFEntry const& gtf_entry) const
{
    CHAR_SET::iterator transcript_iter;
    GENE_TRANSCRIPT_MAP::const_iterator transcripts_iter;

    transcripts_iter = this->transcripts.find(gtf_entry.gene_id);
    assert(transcripts_iter != this->transcripts.end());
   
    transcript_iter = (*transcripts_iter).second.find(gtf_entry.transcript_id);
    assert(transcript_iter != (*transcripts_iter).second.end());
    
    return std::distance(transcript_iter, (*transcripts_iter).second.end());
}



bool unique_gene_description::operator<(unique_gene_description const& u) const
{
    return this->contig < u.contig ||
        (this->contig == u.contig &&
         (this->gene < u.gene ||
          (this->gene == u.gene &&
           (this->strand < u.strand))));
}
    

/*

//project a query alignment onto a collapsed contig, where all
//stretches of non-exonic region are shrunk to a given size
Cigar::CIGAR_VEC
MetaGene::CollapseToContig(char const* query_contig, 
                           Cigar::CIGAR_VEC const& query_alignment)
{
    //1. find all overlapping gene intervals on current contig.
    std::vector<IntervalTree const*> overlaps =
        this->meta_exon_tree->overlap(query_contig, 
                                      query_start_bound, 
                                      query_end_bound);

    std::vector<IntervalTree const*>::iterator overlap_iter;

    if (overlaps.empty())
    {
        //this locus falls outside of any annotated exon
        return false
    }
    IntervalTree const* overlap = *overlap_iter;
    GENE_ITERS & genes = *static_cast<GENE_ITERS *>(overlap->payload);

    assert(! genes.empty());

    return ProjectSingleRegion(genes.begin(), 
                               query_start_bound,
                               query_end_bound,
                               use_absolute_coords,
                               proj_start_bound,
                               proj_end_bound);
}



//project a query region onto each containing collapsed meta_gene
std::map<unique_gene_description, std::pair<size_t, size_t> >
MetaGene::CollapseToGenes(char const* query_contig, 
                          size_t query_start_bound,
                          size_t query_end_bound)
{
    
    std::map<unique_gene_description, std::pair<size_t, size_t> > results;
    
    //1. find all overlapping gene intervals on current contig.
    std::vector<IntervalTree const*> overlaps =
        this->gene_extent_tree->overlap(query_contig, 
                                        query_start_bound, 
                                        query_end_bound);
    
    if (overlaps.empty())
    {
        //this locus falls outside of any annotated exon
        return results;
    }
    
    //projects the locus into the collapsed meta-exon structure of
    //each containing gene
    bool success;
    
    std::vector<IntervalTree const*>::iterator overlap_iter;
    for (overlap_iter = overlaps.begin(); 
         overlap_iter != overlaps.end(); ++overlap_iter)
    {
        IntervalTree const* overlap = *overlap_iter;
        GENE_ITERS & genes = *static_cast<GENE_ITERS *>(overlap->payload);
        
        for (size_t i = 0; i != genes.size(); ++i)
        {
            success = ProjectRegion(genes[i], 
                                    query_start_bound,
                                    query_end_bound,
                                    false,
                                    proj_start_bound,
                                    proj_end_bound);
            if (success)
            {
                results[(*genes[i]).first] = 
                    std::make_pair(proj_start_bound, proj_end_bound);
            }
        }
    }
    return results;
}


bool MetaExon::ProjectSingleRegion(GENE_INDEX::iterator gene_iter,
                                   size_t query_start_bound,
                                   size_t query_end_bound,
                                   bool absolute_coords,
                                   size_t * proj_start_bound,
                                   size_t * proj_end_bound)
{
    unique_gene_description const& gene_desc = (*gene_iter).first;
    BoundsInfo const& gene_bounds = (*gene_iter).second;

    assert(this->meta_exon.find(gene_desc.contig) != meta_exon.end());
    Cigar::CIGAR_VEC const& contig_cigar = this->meta_exon[gene_desc.contig];
    std::vector<size_t> const& contig_offsets = this->contig_offsets[gene_desc.contig];

    size_t start_index = gene_bounds.start_cigar_index;
    Cigar::CIGAR_ITER cigar_start_iter = contig_cigar.begin() + start_index;
    assert((*cigar_start_iter).op == Cigar::M);

    size_t proj_offset = use_absolute_coords ? contig_offsets[start_index] : 0;

    *proj_start_bound = proj_offset +
        Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                            query_start_bound,
                            true, &start_is_missing);

    *proj_end_bound = proj_offset +
        Cigar::ProjectCoord(cigar_start_iter, contig_cigar.end(),
                            query_end_bound
                            true, &end_is_missing);

    
    return ! (start_is_missing || end_is_missing);
}

*/
