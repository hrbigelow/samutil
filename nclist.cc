#include "nclist.h"

#include <cstdlib>
#include <cassert>

IntervalTree::IntervalTree(char const* _contig, size_t _start, size_t _end, void * _payload) : 
    has_interval(true), start(_start), end(_end), payload(_payload)
{
    contig = new char[strlen(_contig) + 1];
    strcpy(contig, _contig);

    children = new BY_START();
    assert(children != NULL);
}
    

IntervalTree::IntervalTree() : has_interval(false), contig(NULL), 
                               start(0), end(0), payload(NULL)
{

    children = new BY_START();
    assert(children != NULL);
}


// IntervalTree::IntervalTree(IntervalTree const& r){
// 	cout<<"Error shouldn't be here"<<endl;
// }

IntervalTree::~IntervalTree() 
{
        
    assert(children != NULL);

    if (contig != NULL)
    {
        delete contig;
    }

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


IntervalTree::IntervalTree(IntervalTree const&)
{
    assert(false);
}
IntervalTree & IntervalTree::operator=(IntervalTree const&)
{
    assert(false);
}
IntervalTree::IntervalTree(IntervalTree &)
{
    assert(false);
}
IntervalTree & IntervalTree::operator=(IntervalTree &)
{
    assert(false);
}


bool check_warn_contains(IntervalTree const* parent,
                         IntervalTree const* child)
{
    bool contains = true;
    if (parent->has_interval)
    {
        contains = node_contains(parent, child);
        if (! contains)
        {
            fprintf(stderr, "Error: parent node interval %s:%Zu-%Zu "
                    "doesn't contain its child interval %s:%Zu-%Zu.\n",
                    parent->contig, parent->start, parent->end,
                    child->contig, child->start, child->end);
        }
    }
    return contains;
}




//iterates through the tree by depth-first search and in start order
//check well-ordered-ness and that all intervals are on the same contig
bool check_tree(IntervalTree const* root)
{

    bool valid_tree = true;

    if (root->children == NULL) {
        printf("children is NULL\n");
    }

    if (root->children->empty()) return valid_tree;

    BY_START::iterator 
        cur = root->children->begin(), 
        pre = root->children->begin();

    //check that the root, if it actually has an interval, the
    //interval contains all child intervals
    valid_tree = valid_tree && check_warn_contains(root, (*cur).second);

    valid_tree = valid_tree && check_tree((*cur).second);


    //check that all child intervals are well ordered and on same contig.
    cur++;
    while (cur != root->children->end()){

        IntervalTree const * child = (*cur).second;
        valid_tree = valid_tree && check_warn_contains(root, child);

        valid_tree = valid_tree && check_tree(child);
        IntervalTree const* ptree = (*pre).second;

        if (! well_ordered(ptree, child))
        {
            fprintf(stderr, "intervals not well ordered: ");
            fprintf(stderr, "%s: %Zu=>%Zu and %s: %Zu=>%Zu\n", 
                   (*child).contig, (*child).start, (*child).end, 
                   (*ptree).contig, (*ptree).start, (*ptree).end);
            valid_tree = false;
        }

        cur++;
        pre++;
    }
    return valid_tree;
}


void print_tree(IntervalTree const* node, int indent_level)
{
    if (node->has_interval)
    {
        for (int i = 0; i != indent_level; ++i)
        {
            putchar('|');
            putchar(' ');
        }
        printf("%s: %Zu->%Zu\n", (*node).contig, (*node).start, (*node).end);
    }
    assert(node->children != NULL);
    for (BY_START::iterator child = (*node->children).begin();
         child != (*node->children).end(); ++child)
    {
        print_tree((*child).second, indent_level+1);
    }
}

/*
      
  Creates a new node with (contig, start, end) as the interval, and inserts it
  in the proper place in the tree.

  Rule RIGHTMOST_PARENT: The tree is organized so that the unique
  parent of a node shall be the right-most interval (if traversed
  depth-first, sibling order) that *contains* the node's interval.
      
  Finds the lower-bounding subnode by searching the node's children,
  which are sorted on the interval start position.  This is the 'first
  node that is not less than the query node'.  Call this 'R', for
  'right'.  Find its sibling to the left, call that 'L'.  Since 'R' is
  the first node not-less-than the query, L is guaranteed to be less
  than the query, or not exist.

  Then, one of four things is true:
  @ means 'contains'
  ^ means 'is a sibling of'

Derivation:

1.   (Q.end <= L.end)         : L @ Q
2.   (L.end < Q.end)          : L ^ Q
  2a. (Q.start < R.start)
    2a1. (Q.end < R.end)      : Q ^ R
    2a3. (R.end <= Q.end)     : Q @ R 
  2b. (Q.start == R.start)
    2b1. (R.end <= Q.end)     : Q @ R
    2b2. (Q.end < R.end)      : R @ Q


Refinement (this is how the code is organized)

0a       empty                : insert Q
0b.      only R exists        : Regard L ^ Q as true.
0c.      only L exists        : Regard Q ^ R as true.

1.       (Q.end <= L.end)     : L @ Q
2.       (L.end < Q.end)      : L ^ Q
  2a.    (R.end <= Q.end)     : Q @ R
  2b.    (Q.end < R.end)
    2b1. (Q.start < R.start)  : Q ^ R
    2b2. (Q.start == R.start) : R @ Q


Now, having added in the 'contig' element in these intervals, this is
how the code should be:

R will be the first (left-most) node that is NOT less than Q.  So, this is true:

Q.contig < R.contig || (Q.contig == R.contig && Q.start <= R.start)

(R will be 'end' if no such node is found)
L will be the node just to the left of R, or 'end' if it exists.

If R is not found, children.last() is less than Q by definition.








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

//insert the requested interval into the tree.
//return an iterator pointing to the new node, and/or a bool denoting
//whether successfully inserted.
std::pair<BY_START::iterator, bool>
IntervalTree::insert_rec(char const* query_contig, size_t query_start, 
                         size_t query_end, void * payload)
{
    BY_START::iterator const& end = children->end();
    //find the first node, 'R' such that Q.start <= R.start
    //find (left, right) such that: L.start < Q.start <= R.start

    BY_START::iterator right = 
        children->lower_bound(ContigCoord(query_contig, query_start));

    BY_START::iterator left = right;
    if (left != children->begin() && ! children->empty())
    {
        --left;
    }
    else
    {
        //'right' is at the beginning. there is no preceding interval.  
        //set left = end to signal this.
        assert(right == children->begin());
        left = end;
    }

    IntervalTree * left_node = (left == children->end()) ? NULL : (*left).second;
    IntervalTree * right_node = (right == children->end()) ? NULL : (*right).second;
    IntervalTree * query_node = new IntervalTree(query_contig, query_start, query_end, payload);

    if (node_contains(left_node, query_node)) // L @ Q
    {
        delete query_node;
        return (*left).second->insert_rec(query_contig, query_start, query_end, payload);
    }
        
    else // L does not contain Q
    { 
        if (node_contains(query_node, right_node)) // Q @ R
        {
            BY_START::iterator end_range = right;
            for ( ; end_range != end; ++end_range)
            {
                if (! node_contains(query_node, (*end_range).second))
                {
                    break;
                }
            }
            //splice Q between this node and all of the nodes it contains, starting at R
            query_node->children->insert(right, end_range);
            children->erase(right, end_range);
            return children->insert(std::make_pair(ContigCoord(query_contig, query_start), 
                                                   query_node));
        }
        else // Q does not contain R
        {
            if (node_contains(right_node, query_node)) // R @ Q
            {
                delete query_node;
                return (*right).second->insert_rec(query_contig, query_start, query_end, payload);
            }
            else // Q ^ R
            {
                return children->insert
                    (std::make_pair(ContigCoord(query_contig, query_start), query_node));
            }
        }
    }
}


std::pair<BY_START::iterator, bool>
IntervalTree::insert(char const* contig, size_t start, size_t end, 
                     void * payload){
    assert(! this->has_interval);
    return insert_rec(contig, start, end, payload);
}


//What information do you need to determine whether two intervals pass
//a secondary filter?

//recursively load <found> with intervals that overlap the query interval <q>
//in the intervals in root and its descendants
void IntervalOverlap_rec(IntervalTree const* root, 
                         char const* contig, 
                         size_t start, size_t end,
                         std::vector<IntervalTree const*> * found)
{
    std::pair<BY_START::const_iterator, 
              BY_START::const_iterator> bounds = 
        RangeQuery(root->children, contig, start, end);
        
    for (BY_START::const_iterator n = bounds.first; n != bounds.second; ++n) 
    {
        found->push_back((*n).second);
        IntervalOverlap_rec((*n).second, contig, start, end, found);
    }
}
    

std::vector<IntervalTree const*>
IntervalTree::overlap(char const* contig, size_t start, size_t end){
    assert(! this->has_interval);
    std::vector<IntervalTree const*> found;
    IntervalOverlap_rec(this, contig, start, end, &found);
    return found;
}



//find all nodes whose interval overlaps query, by reverse linear search
std::pair<BY_START::const_iterator, BY_START::const_iterator> 
RangeQuery(BY_START const* nodes, 
           char const* contig, size_t start, size_t end){

    assert(nodes != NULL);

    //'right' is the first node that does NOT overlap query
    BY_START::const_iterator right = nodes->lower_bound(ContigCoord(contig, end));
    BY_START::const_iterator left = right;
        
    while (left != nodes->begin())
    {
        --left;
        if (strcmp((*left).second->contig, contig) != 0
            || (*left).second->end <= start)
        {
            //aha!  we've just found our first 'left' node that does not overlap query.
            //take its immediate neighbor to the right (++left)
            ++left;
            break;
        }
    }

    return std::make_pair(left, right);
}



//tell whether a contains b.  By convention, if either are NULL,
//returns false
bool node_contains(IntervalTree const* a, IntervalTree const* b)
{
    return a != NULL 
        && b != NULL
        && ((strcmp(a->contig, b->contig) == 0)
            && a->start <= b->start
            && b->end <= a->end);
}


//tell whether a neighbors b to its left.  By convention, if either
//are NULL, returns false.
bool well_ordered(IntervalTree const* a, IntervalTree const* b)
{
    if (a == NULL || b == NULL)
    {
        return false;
    }
    else
    {
        int cmp = strcmp(a->contig, b->contig);
    return (cmp < 0
            || (cmp == 0
                && a->start < b->start
                && a->end < b->end));
    }
}
