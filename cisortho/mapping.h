#include <cmath>

#include "permutation_collection.h"
#include "conservation.h"

class MappingData {
  
 public:
  int target_index; //temporary placeholder during extend_map
  bool target_used;
  float similarity;
  int offset;
  bool same_orientation;
  bool in_order;  

 MappingData() : 
  target_index(0), target_used(false), similarity(0), offset(0),
    same_orientation(true), in_order(true) {}
};


class ComponentScoreData {

 public:
  float permutation_score;
  float offset_score;
  float orientation_score;
  float similarity_score;
  ComponentScoreData(float ps, float os, float ors, float ss) :
    permutation_score(ps), offset_score(os), orientation_score(ors),
    similarity_score(ss) { }

};




template <typename Collection>
float ApplyCPDCurve(ComponentScoreData const& components){
  
  float offset_odds = 
    cis::InterpolatePoint(Collection::offset_curve,
                          components.offset_score, -HUGE_VAL);
  
  float permutation_odds = 
    cis::InterpolatePoint(Collection::permutation_curve,
                          components.permutation_score, -HUGE_VAL);
  
  float similarity_odds = 
    cis::InterpolatePoint(Collection::hits_similarity_curve,
                          components.similarity_score, -HUGE_VAL);
  
  float orientation_odds = 
    cis::InterpolatePoint(Collection::orientation_curve,
                          components.orientation_score, -HUGE_VAL);
  
  return
    offset_odds + permutation_odds + orientation_odds + similarity_odds;
  
}





template <typename Collection>
void MaximalScores(ComponentScoreData const& component_scores,
                   MappingData * mapping, 
                   int query_size,
                   MappingData * best_mapping,
                   float * best_summary_score){
  
  float summary_score = ApplyCPDCurve<Collection>(component_scores);

  //record this as the best mapping if it is
  if (summary_score > *best_summary_score){
    for (int qi = 0; qi != query_size; ++qi)
      best_mapping[qi] = mapping[qi];

    *best_summary_score = summary_score;
  }

}


template <typename Collection>
void MaximalScores(ComponentScoreData const& component_scores,
                   MappingData * mapping, 
                   int query_size,
                   MappingData * /*best_mapping*/,
                   ComponentScoreData * best_components){

  if (best_components->permutation_score < component_scores.permutation_score)
    best_components->permutation_score = component_scores.permutation_score;
  
  if (best_components->offset_score < component_scores.offset_score)
    best_components->offset_score = component_scores.offset_score;
  
  if (best_components->orientation_score < component_scores.orientation_score)
    best_components->orientation_score = component_scores.orientation_score;
  
  if (best_components->similarity_score < component_scores.similarity_score)
    best_components->similarity_score = component_scores.similarity_score;
  
}



//should this function be switched over to one in which mapping is iteratively
//updated as each pairing is set?
//in that case, the pairwise fields 'offset', 'similarity', and 'same_orientation'
//would be updated every time the target index was updated
//the field 'in_order' is more complex.  it should stand for whether the component
//is in order with the previous element.
//this means that the 

template <typename Collection, typename BestScore>
  void extend_mapping(Collection const& query, Collection const& target,
                      int qi, MappingData * mapping, 
                      MappingData * best_mapping,
                      BestScore * best_scores){
  
  if (qi >= query.size()) {
    //condition:  'mapping' contains a valid mapping
     // of query and target sites.  the mapping will have -1's in places,
     // these are unmapped query sites.
     
    ComponentScoreData current_scores = cis::ScoreSummary(mapping, query.size());

     //This will optionally compute a log-odds
    MaximalScores<Collection>(current_scores, mapping,
                              query.size(), best_mapping, best_scores);
    
   } else {
     
     for (int ti = 0; ti != target.size(); ++ti) { 

      if ((! mapping[ti].target_used) &&
          Collection::ValidElementPair(query.sites[qi], 
                                       target.sites[ti])) {

        mapping[ti].target_used = true;
        mapping[qi].target_index = ti;
        mapping[qi].similarity =
          OrthologyScore(query.sites[qi], target.sites[ti],
                         &mapping[qi].same_orientation);
        
        mapping[qi].offset =
          target.sites[ti]->start - query.sites[qi]->start;

        mapping[qi].in_order = qi == 0 ? false :
          mapping[qi-1].target_index < mapping[qi].target_index;

        extend_mapping(query, target, qi+1, mapping, best_mapping, best_scores);

        mapping[ti].target_used = false;

      }
    }
  }
}
