#ifndef _DNA_TYPES_H
#define _DNA_TYPES_H

#include <string>

namespace cis {

  class dna_t;
  typedef unsigned char NUC;
  typedef std::basic_string<NUC> dnastring;
  struct less_dna { bool operator()(dna_t const&, dna_t const&) const; }; 
  struct less_dna_ptr { bool operator()(dna_t const*, dna_t const*) const; };

  enum dna_strand { POS , NEG, POS_NEG, NON_STRANDED };

  const NUC X(0);
  const NUC A(1);
  const NUC C(2);
  const NUC M(3);
  const NUC G(4);
  const NUC R(5);
  const NUC S(6);
  const NUC V(7);
  const NUC T(8);
  const NUC W(9);
  const NUC Y(10);
  const NUC H(11);
  const NUC K(12);
  const NUC D(13);
  const NUC B(14);
  const NUC N(15);
  const NUC z(16);

  typedef unsigned short int CONS;

  const CONS mNULL(0);
  const CONS mX(1<<0);
  const CONS mA(1<<1);
  const CONS mC(1<<2);
  const CONS mM(1<<3);
  const CONS mG(1<<4);
  const CONS mR(1<<5);
  const CONS mS(1<<6);
  const CONS mV(1<<7);
  const CONS mT(1<<8);
  const CONS mW(1<<9);
  const CONS mY(1<<10);
  const CONS mH(1<<11);
  const CONS mK(1<<12);
  const CONS mD(1<<13);
  const CONS mB(1<<14);
  const CONS mN(1<<15);

} // namespace cis

#endif //_DNA_TYPES_H
