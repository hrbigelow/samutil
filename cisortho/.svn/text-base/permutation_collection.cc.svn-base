#include "permutation_collection.h"

namespace cis {

  float InterpolatePoint(PermutationCollection::POINTS const& points, 
                         float query_point, float default_value){
    
    if (points.size() == 0) return default_value;
    
    PermutationCollection::POINT last_point(points[0]);
    for (size_t i = 1; i != points.size(); ++i){
      if (query_point >= last_point.first &&
          query_point < points[i].first) {
        float x_fraction = 
          query_point - last_point.first / 
          points[i].first - last_point.first;
      
        float y_fraction =
          x_fraction * (points[i].second - last_point.second);

        return last_point.second + y_fraction;
      }
    }

    return default_value;
    
  }

} // end namespace cis
  
