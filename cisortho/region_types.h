#ifndef _REGION_TYPES_H
#define _REGION_TYPES_H


namespace cis {

  class region;

  typedef region * REG_P;
  typedef region const* REG_PC;

  typedef region & REG_R;
  typedef region const& REG_CR;

  struct less_region_pos;

  struct index_tag {};

} // namespace cis

#endif // _REGIONS_TYPES_H
