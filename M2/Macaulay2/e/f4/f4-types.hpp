/* Copyright 2005-2021, Michael E. Stillman */

#ifndef _F4types_h_
#define _F4types_h_


#include <climits>                   // for INT_MIN
#include "VectorArithmetic.hpp"      // for ElementArray
#include "f4/f4-monlookup.hpp"       // for F4MonomialLookupTableT
#include "f4/moninfo.hpp"            // for MonomialInfo, monomial_word, pac...
#include "f4/varpower-monomial.hpp"  // for varpower_monomials, varpower_mon...
#include "newdelete.hpp"             // for our_new_delete, VECTOR, (gc_allocator)
#include "style.hpp"                 // for LT

#define sizeofspair(s, len) \
  (sizeof(*s) - sizeof(s->lcm) + (len) * sizeof(s->lcm[0]))

enum gbelem_type {
  ELEM_IN_RING,          // These are ring elements
  ELEM_POSSIBLE_MINGEN,  // These are min GB elements which might
                         // also be min gens, in the graded case,
                         // they ARE minimal generators
  ELEM_MIN_GB,           // These are elements which are minimal GB elements
  ELEM_NON_MIN_GB        // These are elements which are not minimal GB elements
};

enum spair_type {
  F4_SPAIR_SPAIR,        // arising from an honest spair
  F4_SPAIR_GCD_ZZ,       // for gbs over the integers
  F4_SPAIR_RING,         // an spair between a generator and a gen of the defining ideal
  F4_SPAIR_SKEW,         // from exterior variables times a monomial
  F4_SPAIR_GEN,          // a generator of the defining ideal
  F4_SPAIR_ELEM          // 
};

enum class SPairType {
  SPair,
  Generator,
  Retired
  // later we would also like GCDZZ, Ring, Skew to handle those cases as well
};

struct GBF4Polynomial 
{
  int len;
  ElementArray coeffs;
  monomial_word *monoms;  // This is all of the monomials written contiguously
};

struct pre_spair 
{
  //enum spair_type type;
  SPairType type;
  int deg1;  // simple degree of quot
  varpower_monomial quot;
  int j;
  bool are_disjoint;
};

// old spair type (linked list-style container)
//struct spair
//{
//  spair *next;
//  spair_type type;
//  int deg; /* sugar degree of this spair */
//  int i;
//  int j;
//  monomial_word lcm[1];  // packed_monomial
//};

struct spair
{
public:
  // spair *next;
  SPairType type;
  int deg; /* sugar degree of this spair */
  int i;
  int j;
  // monomial_word lcm[1];  // packed_monomial
  monomial_word* lcm;  // pointer to a monomial space
  
  spair() : type(SPairType::Retired),deg(INT_MIN),i(-1),j(-1),lcm(nullptr) {}
  spair(SPairType t, int deg0, int i0, int j0, monomial_word* lcm0) :
    type(t), deg(deg0), i(i0), j(j0), lcm(lcm0) {}
};

struct gbelem
{
  GBF4Polynomial f;
  int deg;
  int alpha;  // the homogenizing degree
  gbelem_type minlevel;
};

// typedef VECTOR(gbelem *) gb_array;
typedef std::vector<gbelem *> gb_array;

struct sparse_row
{
  int len;
  ElementArray coeffs;
  int *comps; // of length len, allocated in a memory block.
};

struct row_elem
{
  // Header information
  packed_monomial monom; // pointer, allocated monomial in a memory block
  int elem;

  // The polynomial itself
  int len;
  ElementArray coeffs;
  int *comps; // of length len, no longer allocated in a memory block, but with new[]
};

struct column_elem
{
  packed_monomial monom; // pointer, allocated monomial in a memory block
  int head;  // which row is being used as a pivot for this column.
             // -1 means none, -2 means not set yet
};

struct coefficient_matrix
{
  //typedef VECTOR(row_elem) row_array;
  //typedef VECTOR(column_elem) column_array;
  typedef std::vector<row_elem> row_array;
  typedef std::vector<column_elem> column_array;

  row_array rows;
  column_array columns;
};

typedef int (MonomialInfo::*CompareFunction)(const monomial_word *,
                                             const monomial_word *) const;

class ColumnsSorter
{
 public:
  typedef MonomialInfo::value monomial;
  typedef int value;

 private:
  const MonomialInfo *M;
  //  const coefficient_matrix *mat;
  const coefficient_matrix::column_array &cols;

  static long ncmps;
  static long ncmps0;
  //  int (MonomialInfo::*compareFcn)(const monomial_word *, const monomial_word
  //  *) const;
 public:
  int compare(value a, value b)
  {
    //    ncmps ++;
    return M->compare_grevlex(cols[a].monom, cols[b].monom);
    // return (M->*(M->compare))(mat->columns[a].monom,mat->columns[b].monom);

    //    return (M->*compareFcn)(col[a].monom,col[b].monom);
    //    return (M->*compareFcn)(mat->columns[a].monom,mat->columns[b].monom);
    //    return
    //    M->compare_grevlex(mat->columns[a].monom,mat->columns[b].monom);
  }

  bool operator()(value a, value b)
  {
    // ncmps0 ++;
    return (M->compare_grevlex(cols[a].monom, cols[b].monom) == LT);
    //    return (M->*(M->compare))(mat->columns[a].monom,mat->columns[b].monom)
    //    == LT;
  }

  ColumnsSorter(const MonomialInfo *M0, const coefficient_matrix *mat0)
      : M(M0),
        /* mat(mat0), */ cols(
            mat0->columns) /* , compareFcn(&MonomialInfo::compare_grevlex) */
  {
  }

  long ncomparisons() const { return ncmps; }
  long ncomparisons0() const { return ncmps0; }
  void reset_ncomparisons()
  {
    ncmps0 = 0;
    ncmps = 0;
  }

  ~ColumnsSorter() {}
};

class PreSPairSorter
{
 public:
  typedef pre_spair *value;

 private:
  static long ncmps;

 public:
  int compare(value a, value b)
  {
    ncmps++;
    return varpower_monomials::compare(a->quot, b->quot);
  }

  bool operator()(value a, value b)
  {
    ncmps++;
    return varpower_monomials::compare(a->quot, b->quot) == LT;
  }

  PreSPairSorter() {}
  void reset_ncomparisons() { ncmps = 0; }
  long ncomparisons() const { return ncmps; }
  ~PreSPairSorter() {}
};

class SPairCompare
{
public:
  SPairCompare(const std::vector<spair> &spairs) : mSPairs(spairs) { }

public:
  bool operator()(size_t s, size_t t) const
  {
    const spair& a = mSPairs[s];
    const spair& b = mSPairs[t];
    if (a.deg > b.deg) return true;
    if (a.deg < b.deg) return false;
    if (a.i > b.i) return true;
    return false;
  }

private:
  const std::vector<spair>& mSPairs;
};

typedef F4MonomialLookupTableT<int32_t> MonomialLookupTable;

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
