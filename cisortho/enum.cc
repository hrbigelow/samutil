#include "enum.h"
#include "dna_types.h"

namespace ce = cis_enum;

namespace cis_enum {

  int from_string(std::string const& s){
    
    int r = 0;
    if (s == "POS" || s == "+") r = cis::POS;
    else if (s == "NEG" || s == "-") r = cis::NEG;
    else if (s == "POS_NEG") r = cis::POS_NEG;
    else if (s == "NON_STRANDED" || s == ".") r = cis::NON_STRANDED;
    else if (s == "HEAD") r = ce::HEAD;
    else if (s == "TAIL") r = ce::TAIL;
    else if (s == "EXON") r = ce::EXON;
    else if (s == "PARTIAL_EXON") r = ce::PARTIAL_EXON;
    else if (s == "INTRON") r = ce::INTRON;
    //     else if (s == "HEADS") r = ce::HEADS;
    //     else if (s == "TAILS") r = ce::TAILS;
    //     else if (s == "TANDEM") r = ce::TANDEM;
    //     else if (s == "INTRONIC") r = ce::INTRONIC;
    //     else if (s == "EXONIC") r = ce::EXONIC;
    //     else if (s == "MIXED_INTRONIC") r = ce::MIXED_INTRONIC;
    //     else if (s == "INTERGENIC_INTRONIC") r = ce::INTERGENIC_INTRONIC;
    //     else if (s == "UNKNOWN") r = ce::UNKNOWN;
    //     else if (s == "NA") r = ce::NA;
    //     else if (s == "UPSTREAM") r = ce::UPSTREAM;
    //     else if (s == "IN") r = ce::IN;
    //     else if (s == "DOWNSTREAM") r = ce::DOWNSTREAM;

    return r;
  }


  std::string print(ce::rel_position d){
    std::string r;
    switch(d){
    case ce::HEAD: r = "HEAD"; break;
    case ce::TAIL: r = "TAIL"; break;
    case ce::EXON: r = "EXON"; break;
    case ce::PARTIAL_EXON: r = "PARTIAL_EXON"; break;
    case ce::INTRON: r = "INTRON"; break;
    }
    return r;
  }

//   std::string print(ce::syn_type d){
//     std::string r;
//     switch(d){
//     case ce::HEADS: r = "HEADS"; break;
//     case ce::TAILS: r = "TAILS"; break;
//     case ce::TANDEM: r = "TANDEM"; break;
//     case ce::INTRONIC: r = "INTRONIC"; break;
//     case ce::EXONIC: r = "EXONIC"; break;
//     case ce::MIXED_INTRONIC: r = "MIXED_INTRONIC"; break;
//     case ce::INTERGENIC_INTRONIC: r = "INTERGENIC_INTRONIC"; break;
//     case ce::UNKNOWN: r = "UNKNOWN"; break;
//     case ce::NA: r = "NA"; break;
//     }
//     return r;
//   }

//   std::string print(ce::offset_type d){
//     std::string r;
//     switch(d){
//     case ce::UPSTREAM: r = "UPSTREAM"; break;
//     case ce::IN: r = "IN"; break;
//     case ce::DOWNSTREAM: r = "DOWNSTREAM"; break;
//     }
//     return r;
//   }


} // namespace cis_enum
