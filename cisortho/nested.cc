#include "nested.h"

#include <cstdlib>
#include <cassert>
#include <cstdio>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>


namespace cis {

    //compares two BY_START iterators safely such that if either points
    //to the end, handles this correctly
    bool less_iterator(BY_START::iterator & a, BY_START::iterator & b, 
                       BY_START::iterator & end)
    {
    
        if (a == b) 
        {
            return false;
        }
        else 
        {
            if (a == end)
            {
                return false;
            }
            else 
            {
                if (b == end) 
                {
                    return true;
                }
                else 
                {
                    return (*a).first < (*b).first;
                }
            }
        }

    }


    //   bool less_interval_start::operator()(region_tree const* a, 
    //                                        region_tree const* b) const {
    //     return a->interval->start < b->interval->start;
    //   }

    //   bool less_interval_end::operator()(region_tree const* a, 
    //                                      region_tree const* b) const {
    //     return a->interval->end < b->interval->end;
    //   }


    region_tree::region_tree(region const* r) : interval(r) 
    {
        children = new BY_START();
        assert(children != NULL);
    }
    

    // region_tree::region_tree(region_tree const& r){
    // 	cout<<"Error shouldn't be here"<<endl;
    // }

    region_tree::~region_tree() 
    {
        
        assert(children != NULL);
        
        BY_START::iterator s_iter;
        for (s_iter = (*children).begin();
             s_iter != (*children).end(); ++ s_iter)
        {
            delete (*s_iter).second;
        }

        delete children;

        children = NULL;

        //cout<<"In DTor"<<endl;
    }


    //iterates through the tree by depth-first search and in start order
    //check well-ordered-ness and that all intervals are on the same contig
    void check_tree(region_tree const* t){

        if (t->children == NULL) {
            printf("children is NULL\n");
        }

        if (t->children->empty()) return;

        BY_START::iterator 
        cur = t->children->begin(), 
        pre = t->children->begin();

        check_tree((*cur).second);
        cur++;
        while (cur != t->children->end()){

            region_tree const * child = (*cur).second;
            check_tree(child);
            region const& c = *child->interval;
            region const& p = *((*pre).second)->interval;

            if (! (p.start < c.start && p.end < c.end)){
                printf("intervals not well ordered: ");
                printf("%" PRId64 "", c.start);
                printf("%"PRId64"=>%"PRId64" and %"PRId64"=>%"PRId64"\n", 
                       c.start, c.end, p.start, p.end);
            }

            cur++;
            pre++;
        }
    }


    void print_tree(region_tree const* node, int indent_level)
    {
        if (node->interval != NULL)
        {
            region const& iv = *node->interval;
            for (int i = 0; i != indent_level; ++i)
            {
                putchar('|');
                putchar(' ');
            }
            printf("%i[%s:%"PRId64"->%"PRId64"]\n", iv.id, iv.dna.name.c_str(), iv.start, iv.end);
        }
        assert(node->children != NULL);
        for (BY_START::iterator child = (*node->children).begin();
             child != (*node->children).end(); ++child)
        {
            print_tree((*child).second, indent_level+1);
        }
    }

    /*
      
      Creates a new node with <r> as the interval, and inserts it in the
      proper place in the tree.

      Rule RIGHTMOST_PARENT: The tree is organized so that the unique
      parent of a node shall be the right-most interval (if traversed
      depth-first, sibling order) that contains the node's interval.
      
      Finds the lower-bounding subnode by searching the node's children,
      which are sorted on the interval start position.  This is the
      'first node that is not less than the query node'.  Call this 'R',
      for 'right'.  Find its sibling to the left, call that 'L'.

      Find the last node L such that L.start < Q.start
      Take R as the successor to L (hence, Q.start < R.start)
      
      Then, one of three things is true:
      @ means 'contains'
      ^ means 'is a sibling of'

      1.   L @ Q (Q.end <= L.end)
      2.   L ^ Q (L.end < Q.end)
         2a.  Q ^ R (Q.end < R.end)
         2
      
      Actions:
      1.   insert_rec Q into L
      2a.  map::insert Q at R
      2b.  scan forward until the first sibling of Q is found.
           map::insert the range into Q.
           delete from current node
           map::insert Q into current node at R
      
      Comments:
      By construction, Q cannot be a child of R since Q.start < R.start
      Even though there may be other nodes before L that also contain Q,
      we don't consider them as possible parents, for efficiency reasons.

    */

    void region_tree::insert_rec(region const* query)
    {


        BY_START::iterator const& end = children->end();
        //find the first node 'right' such that query.start < right->start
        BY_START::iterator right = children->lower_bound(query->start);
        BY_START::iterator left = right;
        if (left != children->begin())
        {
            --left;
        }
        else
        {
            left = end;
        }

        if (left != end && query->end <= (*left).second->interval->end)
        {
            //1.  L @ Q
            (*left).second->insert_rec(query);
        }
        
        else 
        { 
            //L ^ Q

            if (right != end && query->start == (*right).second->interval->start)
            {
                //containment relationship, either Q @ R or R @ Q
                if ((*right).second->interval->end <= query->end)
                {
                    // Q @ R
                    BY_START::iterator end_range = right;
                    for ( ; end_range != end; ++end_range)
                    {
                        if (query->end < (*end_range).second->interval->end)
                        {
                            break;
                        }
                    }
                    //create a new node and store appropriately
                    region_tree * query_node = new region_tree(query);
                    query_node->children->insert(right, end_range);
                    children->erase(right, end_range);
                    children->insert(std::make_pair(query->start, query_node));
                }
                else
                {
                    // R @ Q
                    (*right).second->insert_rec(query);
                }
            }
            else
            {
                // Q ^ R
                children->insert
                    (right, std::make_pair(query->start, 
                                           new region_tree(query)));
            }
        }
    }


    void region_tree::insert(region const* r){
        assert(interval == NULL);
        insert_rec(r);
    }


    //What information do you need to determine whether two intervals pass
    //a secondary filter?

    //recursively load <found> with intervals that overlap the query interval <q>
    //in the intervals in root and its descendants
    void IntervalOverlap_rec(region_tree const* root, REG_CR query_region, 
                             REGIONS_MULTI & found)
    {
        std::pair<BY_START::const_iterator, 
        BY_START::const_iterator> bounds = 
        RangeQuery(root->children, query_region);
        
        for (BY_START::const_iterator n = bounds.first; n != bounds.second; ++n) 
        {
            found.insert((*n).second->interval);
            IntervalOverlap_rec((*n).second, query_region, found);
        }
    }
    

    REGIONS_MULTI IntervalOverlap(region_tree const& root, REG_CR query_region){
        assert(root.interval == NULL);
        REGIONS_MULTI found;
        IntervalOverlap_rec(&root, query_region, found);
        return found;
    }



    //find all nodes whose interval overlaps query, by reverse linear search
    std::pair<BY_START::const_iterator, BY_START::const_iterator> 
    RangeQuery(BY_START const* nodes, region const& query){

        assert(nodes != NULL);

        //'right' is the first node that does NOT overlap query
        BY_START::const_iterator right = nodes->lower_bound(query.end);
        BY_START::const_iterator left = right;
        
        while (left != nodes->begin())
        {
            --left;
            if ((*left).second->interval->end < query.start)
            {
                ++left;
                break;
            }
        }

        return std::make_pair(left, right);
    }


    //build a tree from a set of regions...
    region_tree * BuildTree(REGIONS_MULTI const& regs){
        if (regs.empty()) 
        {
            return NULL;
        }

        REGIONS_MULTI::const_iterator r = regs.begin();
        dna_t const& contig = (*r)->dna;

        region_tree * root = new region_tree(NULL);

        while (r != regs.end()) 
        {
            if (contig != (*r)->dna)
            {
                fprintf(stderr, "Error: BuildTree: All regions must be on the same contig\n"
                        "First region is on %s, current region on %s\n",
                        contig.name.c_str(), (*r)->dna.name.c_str());
                exit(1);
            }
            root->insert(*r++);
        }

        check_tree(root);
        return root;
    }


    //build a map of trees, each on separate contigs
    cis::TREE_MAP BuildTrees(REGIONS_MULTI const& regions)
    {
        cis::REG_MAP contig_to_regions = cis::SplitRegions(regions);
        cis::REG_MAP::iterator iter;
        cis::TREE_MAP trees;
        for (iter = contig_to_regions.begin(); iter != contig_to_regions.end(); ++iter)
        {
            dna_t const* contig = (*iter).first;
            REGIONS_MULTI const& regions_on_contig = (*iter).second;
            trees[contig] = BuildTree(regions_on_contig);
        }
        return trees;
    }

} // namespace cis
