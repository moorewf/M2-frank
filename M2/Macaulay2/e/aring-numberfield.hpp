// Copyright 2014 Michael E. Stillman

#ifndef _aring_gf_flint_hpp_
#define _aring_gf_flint_hpp_

#include <vector>

// The following needs to be included before any flint files are included.
#include <M2/gc-include.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/flint.h>      // for ???
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>       // for ???
#include <flint/fmpq_poly.h>  // for ???
#include <flint/nf.h>         // for ???
#include <flint/nf_elem.h>    // for ???
#pragma GCC diagnostic pop

#include "interface/random.h"
#include "aring.hpp"
#include "buffer.hpp"
#include "ringelem.hpp"
#include "exceptions.hpp"           // for exc::division_by_zero_error

class PolynomialRing;
class RingElement;

namespace M2 {

// types for translation between Flint poly types and ring_elems.
using ZZPolyStruct = std::vector<gmp_ZZ>;
struct QQPolyStruct
{
   std::vector<gmp_ZZ> coeffs;
   gmp_ZZ denom;
};
using ZZPoly = ZZPolyStruct&;
using QQPoly = QQPolyStruct&;

// Some helper functions to convert to/from arrays of gmp integers and flint's types

fmpz_poly_struct toFlintZZPoly(ZZPoly f)
{
   // we assume that result has already been initialized
   fmpz_poly_struct result;
   fmpz_poly_init2(&result, f.size());

   //   fmpz* resultCoeffs = fmpz_poly_get_coeff_ptr(&result);

   for (auto i = 0; i < f.size(); ++i)
   {
     fmpz_set_mpz(fmpz_poly_get_coeff_ptr(&result, i), f[i]);
     //      fmpz_set_mpz(&(resultCoeffs[i]), f[i]);
   }
   return result;
}

fmpq_poly_struct toFlintQQPoly(QQPoly f)
{
   // we assume that result has already been initialized
   fmpq_poly_struct result;
   fmpz_poly_init2(&result, f.coeffs.size());

   fmpz* resultCoeffs = fmpq_poly_numref(result);
   fmpz_t resultDen = fmpq_poly_denref(result);
   
   for (auto i = 0; i < f.coeffs.size(); ++i)
   {
      fmpz_set_mpz(&(resultCoeffs[i]), f.coeffs[i]);
   }
   fmpz_set_mpz(resultDen, f.denom);
   return result;
}

void fromFlintZZPoly(ZZPoly result, fmpz_poly_t f)
{
   // TODO
   result.resize(fmpz_poly_length(&f));
}

void fromFlintQQPoly(QQPoly result, fmpq_poly_t f)
{
   // TODO
   result.coeffs.resize(fmpz_poly_length(&f));   
}

/**
\ingroup rings
*/

class ARingNumberField : public RingInterface
{
 public:
  static const RingID ringID = ring_NumberField;
  typedef nf_elem_struct ElementType;
  typedef ElementType elem;
  typedef std::vector<elem> ElementContainerType;
  using ReadOnlyElement = ReadOnlyElementTempl<ARingNumberField>;

  /**
   * \brief A wrapper class for ElementType
   *
   * This keeps a pointer to the fq_zech_ctx_struct as it's needed to
   * implement the destructor
   */
  class Element : public ElementImpl<ElementType>
  {
   public:
    Element() = delete;
    Element(Element&& other) : mContext(other.mContext)
    {
      // figure out how to move the value without the context
      nf_elem_init(&mValue, mContext);
      nf_elem_set(&mValue, &other.mValue, mContext);
    }
    explicit Element(const ARingNumberField& R) : mContext(R.mContext)
    {
      nf_elem_init(&mValue, mContext);
    }
    Element(const ARingNumberField& R, const ElementType& value) : mContext(R.mContext)
    {
      R.init_set(mValue, value);
    }
    ~Element() { nf_elem_clear(&mValue, mContext); }

   private:
    const nf_struct* mContext;
  };

  class ElementArray
  {
    const nf_struct& mContext;
    const size_t mSize;
    std::unique_ptr<ElementType[]> mData;

   public:
    ElementArray(const ARingNumberField& R, size_t size)
        : mContext(R.mContext), mSize(size), mData(new ElementType[size])
    {
      for (size_t i = 0; i < mSize; i++) nf_elem_init(&mData[i], mContext);
    }
    ~ElementArray()
    {
      for (size_t i = 0; i < mSize; i++) nf_elem_clear(&mData[i], mContext);
    }
    ElementType& operator[](size_t idx) { return mData[idx]; }
    const ElementType& operator[](size_t idx) const { return mData[idx]; }
    ElementType *data() { return mData.get(); }
    const ElementType *data() const { return mData.get(); }
  };

  // TODO: Do we really want to take a ring_elem as input here?
  ARingNumberField(const PolynomialRing& R, const ring_elem a);

  ~ARingNumberField();

  const nf_struct& flintContext() const { return mContext; }

  long characteristic() const { return 0; }
  long dimension() const { return mDimension; }
  const PolynomialRing& originalRing() const { return mOriginalRing; }
  void getGenerator(ElementType& result_gen) const;

  void text_out(buffer& o) const;

 private:
  nf_struct& mContext;

  const PolynomialRing& mOriginalRing;   // This is a quotient ring k[a]/f(a).
  long mDimension;
  mutable flint_rand_t mRandomState;

  ////////////////////////////////
  /// Arithmetic functions ///////
  ////////////////////////////////

 public:
  unsigned int computeHashValue(const elem& a) const
  {
    // TODO: use entire data for the has value
    return static_cast<unsigned int>(a.value);
  }

  void to_ring_elem(ring_elem& result, const ElementType& a) const
  {
    // we only use this data type for GF's smaller than 32 bits, since
    // this class creates lookup tables that are way too big otherwise.

    // TODO
#if 0
    ElementType* b = getmemstructtype(ElementType*);
    init(*b);
    copy(*b, a);
    size_t coeffs_size = sizeof(mp_limb_t)*b->alloc;
    mp_ptr coeffs = reinterpret_cast<mp_ptr>(getmem_atomic(coeffs_size));
    memcpy(coeffs,b->coeffs,coeffs_size);
    flint_free(b->coeffs);
    b->coeffs = coeffs;
    result.poly_val = reinterpret_cast<Nterm*>(b);
#endif
  }

  void from_ring_elem(ElementType& result, const ring_elem& a) const
  {
    // TODO
  }

  ElementType from_ring_elem_const(const ring_elem& a) const
  {
    // TODO
    return *(a.get_number_field_elem());
  }

  void from_ring_elem_const_clear(ElementType a) const
  {
    // nothing needed to do or is that true? TODO
  }

  bool is_unit(const ElementType& f) const { return not is_zero(f); }
  bool is_zero(const ElementType& f) const
  {
    return nf_elem_is_zero(&f, mContext) != 0;
  }
  bool is_equal(const ElementType& f, const ElementType& g) const
  {
    return nf_elem_equal(&f, &g, mContext) != 0;
  }

  int compare_elems(const ElementType& f, const ElementType& g) const; // TODO

  void copy(ElementType& result, const ElementType& a) const
  {
    nf_elem_set(&result, &a, mContext);
  }
  void init(ElementType& result) const { nf_elem_init(&result, mContext); }
  void init_set(ElementType& result, const ElementType& a) const
  {
    init(result);
    copy(result, a);
  }
  void set(ElementType& result, const ElementType& a) const { copy(result, a); }
  void set_zero(ElementType& result) const { nf_elem_zero(&result, mContext); }
  void clear(ElementType& result) const { nf_elem_clear(&result, mContext); }
  void set_from_long(ElementType& result, long a) const
  {
    // TODO
    //fq_zech_set_ui(&result, a1, mContext);
  }

  void set_var(ElementType& result, int v) const
  {
    // TODO
    //if (v != 0) set_from_long(result, 1);
    //std::vector<long> poly = {0, 1};
    //fromSmallIntegerCoefficients(result, poly);
    // printf("variable is %lu\n", result.value);
  }

  void set_from_mpz(ElementType& result, mpz_srcptr a) const
  {
    // TODO
    //int b = static_cast<int>(mpz_fdiv_ui(a, characteristic()));
    //set_from_long(result, b);
  }

  bool set_from_mpq(ElementType& result, mpq_srcptr a) const
  {
    // TODO
#if 0
    ElementType n, d;
    init(n);
    init(d);
    set_from_mpz(n, mpq_numref(a));
    set_from_mpz(d, mpq_denref(a));
    if (is_zero(d)) return false;
    divide(result, n, d);
#endif
    return true;
  }

  bool set_from_BigReal(ElementType& result, gmp_RR a) const { return false; }
  void negate(ElementType& result, const ElementType& a) const
  {
    nf_elem_neg(&result, &a, mContext);
  }

  void invert(ElementType& result, const ElementType& a) const
  {
    if (is_zero(a))
      throw exc::division_by_zero_error();
    nf_elem_inv(&result, &a, mContext);
  }

  void add(ElementType& result,
           const ElementType& a,
           const ElementType& b) const
  {
    nf_elem_add(&result, &a, &b, mContext);
  }

  void subtract(ElementType& result,
                const ElementType& a,
                const ElementType& b) const
  {
    nf_elem_sub(&result, &a, &b, mContext);
  }

  void subtract_multiple(ElementType& result,
                         const ElementType& a,
                         const ElementType& b) const
  {
    ElementType c;
    init(c);
    mult(c, a, b);
    subtract(result, result, c);
    clear(c);
  }

  void mult(ElementType& result,
            const ElementType& a,
            const ElementType& b) const
  {
    nf_elem_mul(&result, &a, &b, mContext);
  }

  void divide(ElementType& result,
              const ElementType& a,
              const ElementType& b) const
  {
    // compute a/b
    if (is_zero(b))
      throw exc::division_by_zero_error();
    nf_elem_div(&result, &a, &b, mContext);
  }

  void power(ElementType& result, const ElementType& a, int n) const
  {
    // TODO: 0^0 == 1?  This is what flint does
    if (n < 0)
      {
        invert(result, a);
        nf_elem_pow(&result, &result, -n, mContext);
      }
    else
      nf_elem_pow(&result, &a, n, mContext);
  }

  void power_mpz(ElementType& result, const ElementType& a, mpz_srcptr n) const
  {
    // TODO: first see if n fits into an int, then call above power function.
  }

  void swap(ElementType& a, ElementType& b) const
  {
    nf_elem_swap(&a, &b, mContext);
  }

  void elem_text_out(buffer& o,
                     const ElementType& a,
                     bool p_one = true,
                     bool p_plus = false,
                     bool p_parens = false) const;

  void syzygy(const ElementType& a,
              const ElementType& b,
              ElementType& x,
              ElementType& y) const
  // returns x,y s.y. x*a + y*b == 0.
  // if possible, x is set to 1.
  // no need to consider the case a==0 or b==0.
  {
    assert(not is_zero(a)); // TODO: change to throw.
    assert(not is_zero(b)); // TODO: change to throw.
    set_from_long(x, 1);
    divide(y, a, b);
    negate(y, y);
  }

  void random(ElementType& result) const
  {
    mp_bitcnt_t bitcount = 5; // TODO: this would be a good parameter to have
    nf_elem_randtest(&result, mRandomState, bitcount, mContext)
  }

  bool promote(const Ring* Rf, const ring_elem f, ElementType& result) const; // TODO

  void lift_to_original_ring(ring_elem& result, const ElementType& f) const; // TODO

  bool lift(const Ring* Rg, const ElementType& f, ring_elem& result) const; // TODO

  // map : this --> target(map)
  //       primelem --> map->elem(first_var)
  // evaluate map(f)
  void eval(const RingMap* map,
            const ElementType& f,
            int first_var,
            ring_elem& result) const; // TODO
}; // end of ARingNumberField definition
}; // end namespace M2

#endif
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
