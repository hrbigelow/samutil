#ifndef _PERMUTATION_COLLECTION_H
#define _PERMUTATION_COLLECTION_H

#include <utility>
#include <vector>

namespace cis {


  class PermutationCollection {
    //weights for combining the four component scores into one

  public:
    typedef std::pair<float, float> POINT;
    typedef std::vector<POINT> POINTS;

    static POINTS offset_curve;
    static POINTS permutation_curve;
    static POINTS orientation_curve;
    static POINTS hits_similarity_curve;


    static void Initialize(POINTS _offset_curve,
                           POINTS _permutation_curve,
                           POINTS _orientation_curve,
                           POINTS _hits_similarity_curve) {
   
      offset_curve = _offset_curve;
      permutation_curve = _permutation_curve;
      orientation_curve = _orientation_curve;
      hits_similarity_curve = _hits_similarity_curve;
    
    }
    
    
  };


//see module.txt For a function f() expressed as a series of points,
//calculate f(x), where x is a point to be interpolated.  returns
//default_value if out of bounds assumes 'points' (x,y) are sorted
//increasing in x
  float InterpolatePoint(PermutationCollection::POINTS const& points, 
                         float query_point, float default_value);
    


} // end namespace cis

#endif // _PERMUTATION_COLLECTION_H
