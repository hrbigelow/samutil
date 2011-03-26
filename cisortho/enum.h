#ifndef _ENUM_H
#define _ENUM_H

#include <string>


namespace cis_enum {

  enum rel_position { HEAD, TAIL, EXON, PARTIAL_EXON, INTRON };
  const int kNumRelPositions = 5;
  
  /* enum syn_type { HEADS = 11, TAILS, TANDEM, INTRONIC, EXONIC, MIXED_INTRONIC,  */
  /* 				INTERGENIC_INTRONIC, */
  /* 				UNKNOWN, NA }; */
  
  /* enum offset_type { UPSTREAM = 20, IN, DOWNSTREAM }; */
  
  int from_string(std::string const& s);
  std::string print(rel_position);
/*   std::string print(syn_type); */
/*   std::string print(offset_type); */
  

} // namespace cis_enum

#endif // _ENUM_H
