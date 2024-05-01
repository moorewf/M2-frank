// Copyright 2014 Michael E. Stillman

#include "aring-numberfield.hpp"
#include <iostream>
// #include "relem.hpp"
// #include "poly.hpp"
// #include "ringmap.hpp"

namespace M2 {

void toFlintZZPoly(fmpz_poly_struct* result, const std::vector<long>& coeffs)
{
   // result should be initialized beforehand
   for (auto i = 0; i < coeffs.size(); ++i)
   {
     fmpz_poly_set_coeff_si(result,i,coeffs[i]);
   }
}

fmpz_poly_struct toFlintZZPoly(ZZPoly f)
{
   fmpz_poly_struct result;
   fmpz tempCoeff;
   fmpz_poly_init2(&result, f.size());
   fmpz_init(&tempCoeff);

   for (auto i = 0; i < f.size(); ++i)
   {
     fmpz_set_mpz(&tempCoeff,f[i]);
     fmpz_poly_set_coeff_fmpz(&result, i, &tempCoeff);
   }

   fmpz_clear(&tempCoeff);

   return result;
}

fmpq_poly_struct toFlintQQPoly(QQPoly f)
{
   fmpq_poly_struct result;
   fmpq_poly_init2(&result, f.coeffs.size());

   fmpz* resultCoeffs = fmpq_poly_numref(&result);
   fmpz* resultDen = fmpq_poly_denref(&result);
   
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
   result.resize(fmpz_poly_length(f));
}

void fromFlintQQPoly(QQPoly result, fmpq_poly_t f)
{
   // TODO
   result.coeffs.resize(fmpq_poly_length(f));   
}

  ARingNumberField::ARingNumberField(const PolynomialRing& R, const ring_elem a)
    // Possible arguments: RingElement f, PolynomialRing + std::vector<gmp_ZZ>, QuotientRing QQ[a]/F(a)
    // Which of these is best?
    : mContext() // TODO
      //    , mOriginalRing(R)
      //    , mRandomState() // TODO
      //    , mPolynomial(a) // a standard vector or element of R?
    , mDimension() // TODO
  {
    // TODO
  }

  ARingNumberField::ARingNumberField(const std::vector<long>& coeffs)
    // Possible arguments: RingElement f, PolynomialRing + std::vector<gmp_ZZ>, QuotientRing QQ[a]/F(a)
    // Which of these is best?
  {
    fmpz_poly_struct temp;
    fmpz_poly_init2(&temp,coeffs.size());
    fmpq_poly_init2(&mFlintPolynomial,coeffs.size());
    toFlintZZPoly(&temp,coeffs);
    fmpq_poly_set_fmpz_poly(&mFlintPolynomial,&temp);
    nf_init(&mContext,&mFlintPolynomial);
    fmpz_poly_clear(&temp);
  }

  ARingNumberField::~ARingNumberField()
  {
    // clear mContext
    // TODO: anything else?
  }

  void ARingNumberField::elem_text_out(buffer& o,
                     const ElementType& a,
                     bool p_one,
                     bool p_plus,
                     bool p_parens) const
  {
    // TODO
  }

  bool ARingNumberField::promote(const Ring* Rf, const ring_elem f, ElementType& result) const
  {
    // TODO
    return false;
  }

  void ARingNumberField::lift_to_original_ring(ring_elem& result, const ElementType& f) const
  {
    // TODO
  }

  bool ARingNumberField::lift(const Ring* Rg, const ElementType& f, ring_elem& result) const
  {
    // TODO
    return false;
  }

  // map : this --> target(map)
  //       primelem --> map->elem(first_var)
  // evaluate map(f)
  void ARingNumberField::eval(const RingMap* map,
            const ElementType& f,
            int first_var,
            ring_elem& result) const
  {
    // TODO
  }
};

#if 0
-- M2 code

R = QQ[x]/(x^3-x-1)
K = numberField R
  F = (ideal R)_0
debug Core  
rawNumberField raw F -- returning the wrong thing.
#endif
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
