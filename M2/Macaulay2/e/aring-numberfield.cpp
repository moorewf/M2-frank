// Copyright 2014 Michael E. Stillman

#include "aring-numberfield.hpp"
// #include "relem.hpp"
// #include "poly.hpp"
// #include "ringmap.hpp"

namespace M2 {

  ARingNumberField::ARingNumberField(const PolynomialRing& R, const ring_elem a)
    // Possible arguments: RingElement f, PolynomialRing + std::vector<gmp_ZZ>, QuotientRing QQ[a]/F(a)
    // Which of these is best?
    : mContext() // TODO
    , mOriginalRing(R)
    , mRandomState() // TODO
    , mPolynomial(a) // a standard vector or element of R?
    , mDimension() // TODO
  {
    // TODO
  }

  ARingNumberField::~ARingNumberField()
  {
    // clear mContext
    // TODO: anything else?
  }

  void ARingNumberField::elem_text_out(buffer& o,
                     const ElementType& a,
                     bool p_one = true,
                     bool p_plus = false,
                     bool p_parens = false) const
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
