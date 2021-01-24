#ifndef __vector_arithmetic_2__
#define __vector_arithmetic_2__

#include "NCAlgebras/Range.hpp"  // for Range
#include "newdelete.hpp"         // for VECTOR
#include "ringelem.hpp"
#include "aring-glue.hpp"
#include <variant>
#include <type_traits>

#include <tbb/tbb.h>

class Ring;
union ring_elem;
using ComponentIndex = int;
class MemoryBlock;

// this is still deriving from our_new_delete, since this
// class doesn't know whether mValue points to a std::vector or a VECTOR.
// this issue persists in the Row struct in NCF4 as well.
class CoeffVector : public our_new_delete
{
  // disallow copy...
 public:
  CoeffVector() : mValue(nullptr) {}
  bool isNull() const { return mValue == nullptr; }
  void swap(CoeffVector& b) { std::swap(mValue, b.mValue); }
  void* getValue() { return mValue; }
  void setValue(void* newValue) { mValue = newValue; }
 private:
  void* mValue;
};

class DenseCoeffVector
{
  // disallow copy...
 public:
  DenseCoeffVector() : mValue(nullptr) {}
  bool isNull() const { return mValue == nullptr; }
  void swap(DenseCoeffVector& b) { std::swap(mValue, b.mValue); }
  void* getValue() { return mValue; }
  void setValue(void* newValue) { mValue = newValue; }
 private:
  void* mValue;
};

template<typename RingType,
         typename CoeffVectorType = CoeffVector,
         typename DenseCoeffVectorType = DenseCoeffVector>
class ConcreteVectorArithmetic
{
public:
  using FieldElement = typename RingType::ElementType;
  using DenseFieldElement = FieldElement;

  using CoeffVectorContainer = typename RingType::ElementContainerType;
  // we may want to change the container type of dense vectors in the future.
  using DenseCoeffVectorContainer = CoeffVectorContainer;

  friend class VectorArithmetic;

  ConcreteVectorArithmetic() : mRing(nullptr) {}

// private:
  const RingType* mRing;

  CoeffVectorContainer* coeffVector(CoeffVectorType f) const
  {
    return reinterpret_cast<CoeffVectorContainer*>(f.getValue());
  }

  DenseCoeffVectorContainer* denseCoeffVector(DenseCoeffVectorType f) const
  {
    return reinterpret_cast<DenseCoeffVectorContainer*>(f.getValue());
  }
  
  CoeffVectorType coeffVector(CoeffVectorContainer* vals) const
  {
    CoeffVectorType result;
    result.setValue(vals);
    vals = nullptr;
    return result;
  }

  DenseCoeffVectorType denseCoeffVector(DenseCoeffVectorContainer* vals) const
  {
    DenseCoeffVectorType result;
    result.setValue(vals);
    vals = nullptr;
    return result;
  }

public:
  ConcreteVectorArithmetic(const RingType* R) : mRing(R) {}
  
  /////////////////////////////
  // Allocation/Deallocation //
  /////////////////////////////
  CoeffVectorType allocateCoeffVector(ComponentIndex nelems) const
  {
    auto tempPtr = new CoeffVectorContainer(nelems);
    for (auto& d : *tempPtr)
      {
	mRing->init(d);
	mRing->set_zero(d);
      }
    return coeffVector(tempPtr);    
  }

  DenseCoeffVectorType allocateDenseCoeffVector(ComponentIndex nelems) const
  {
    auto tempPtr = new DenseCoeffVectorContainer(nelems);
    for (auto& d : *tempPtr)
      {
	mRing->init(d);
	mRing->set_zero(d);
      }
    return denseCoeffVector(tempPtr);    
  }
  // CoeffVector allocateCoeffVector() const
  // {
  //   return coeffVector(new std::vector<FieldElement>);
  // }
  void deallocateCoeffVector(CoeffVectorType& coeffs) const
  {
    delete coeffVector(coeffs);
  }

  void deallocateCoeffVector(DenseCoeffVectorType& coeffs) const
  {
    delete denseCoeffVector(coeffs);
  }
  
  ////////////////////////
  /// Linear Algebra /////
  ////////////////////////

  size_t size(const CoeffVectorType& coeffs) const
  {
    return coeffVector(coeffs)->size();
  }

  size_t size(const DenseCoeffVectorType& coeffs) const
  {
    return denseCoeffVector(coeffs)->size();
  }
  
  // better name: fillDenseRow.
  void sparseRowToDenseRow(DenseCoeffVectorType& dense,
                           const CoeffVectorType& sparse,
                           const Range<int>& comps) const
  {
    // Note: this function simply fills in the values coming from '(sparse, comps)'.
    //   Other values are not touched.
    // In our intended uses, the input `dense` is the vector consisting of all zeros`.
    
    auto& dvec = * denseCoeffVector(dense);
    auto& svec = * coeffVector(sparse);
    
    assert(comps.size() == svec.size());
    assert(comps[comps.size()-1] < dvec.size());
    
    auto len = comps.size();
    for (ComponentIndex i = 0; i < len; i++) mRing->set(dvec[comps[i]],svec[i]);
  }

  void denseRowCancelFromSparse(DenseCoeffVectorType& dense,
                                const CoeffVectorType& sparse,
                                const Range<int>& comps) const
  {
    // ASSUMPTION: svec[0] == 1.
    auto& dvec = * denseCoeffVector(dense);
    auto& svec = * coeffVector(sparse);

    DenseFieldElement a;
    mRing->init_set(a,dvec[comps[0]]);
    for (int i=0; i < comps.size(); ++i)
      mRing->subtract_multiple(dvec[comps[i]], a, svec[i]);
  }

  int denseRowNextNonzero(DenseCoeffVectorType& dense,
                          int first,
                          int last) const
  {
    auto& dvec = * denseCoeffVector(dense);
    for (int i = first; i <= last; i++)
      {
	if (mRing->is_zero(dvec[i])) continue;
	// these lines give the gist of how to handle delayed modulus
	//if (dvec[i] > mCharacteristic) dvec[i] %= mCharacteristic;
	//else if (dvec[i] < 0) dvec[i] %= mCharacteristic;
	//if (dvec[i] == 0) continue;
	return i;
      }
    return last + 1;
  }

  void denseRowToSparseRow(DenseCoeffVectorType& dense,
                           CoeffVectorType& sparse, // output value: sets this value
                           Range<int>& comps, // output value: sets comps
                           int first,
                           int last,
                           MemoryBlock& monomialSpace) const
  { 
    tbb::null_mutex noLock;
    safeDenseRowToSparseRow<tbb::null_mutex>(dense,
					     sparse,
					     comps,
					     first,
					     last,
					     monomialSpace,
					     noLock);
  }

  template<typename LockType>
  void safeDenseRowToSparseRow(DenseCoeffVectorType& dense,
                               CoeffVectorType& sparse, // output value: sets this value
                               Range<int>& comps, // output value: sets comps
                               int first,
                               int last,
                               MemoryBlock& monomialSpace,
                               LockType& lock) const
  {
    auto& dvec = * denseCoeffVector(dense);
    
    int len = 0;
    
    // first can be -1 if the row is zero.  in this case, we should
    // not be accessing dense[i] for i negative.
    for (int i = first; i >= 0 and i <= last; i++)
      {
	if (not mRing->is_zero(dvec[i])) len++;
      }

    comps = monomialSpace.safeAllocateArray<int,LockType>(len,lock);
    sparse = allocateCoeffVector(len);
    auto& svec = * coeffVector(sparse);

    int next = 0;
    for (int i = first; i >= 0 and i <= last; i++)
      if (not mRing->is_zero(dvec[i]))
	{
	  mRing->init_set(svec[next],dvec[i]);
	  comps[next] = i;
	  ++next;
	  mRing->set_zero(dvec[i]);
	}
  }

  void sparseRowMakeMonic(CoeffVectorType& sparse) const
  {
    auto& svec = * coeffVector(sparse);

    FieldElement leadCoeffInv;

    mRing->init(leadCoeffInv);
    mRing->invert(leadCoeffInv,svec[0]);

    for (auto& c : svec) { mRing->mult(c, c, leadCoeffInv); }
  }
  
  template<typename Container>
  void appendSparseVectorToContainer(const CoeffVectorType& sparse,
                                     Container& c) const
  {
    auto& svec = * coeffVector(sparse);
    for (auto a : svec)
      {
	ring_elem tmp;
	mRing->to_ring_elem(tmp,a);
	c.push_back(tmp);
      }
  }

  template<typename Container>
  CoeffVectorType sparseVectorFromContainer(const Container& c) const
  {
    CoeffVectorType sparse = allocateCoeffVector(c.size());
    auto& svec = * coeffVector(sparse);
    for (auto i = 0; i < c.size(); ++i)
      {
	mRing->from_ring_elem(svec[i], c[i]);
      }
    return sparse;
  }

  ring_elem ringElemFromSparseVector(const CoeffVectorType& sparse,
                                     int index) const
  {
    auto& svec = * coeffVector(sparse);
    ring_elem tmp;
    mRing->to_ring_elem(tmp,svec[index]);
    return tmp;
  }

  std::ostream& displayVector(std::ostream& o, const CoeffVectorType& v) const
  {
    auto& svec = * coeffVector(v);

    bool first = true;
    o << "[(" << svec.size() << ")";
    for (const auto& a : svec)
      {
        buffer b;
        mRing->elem_text_out(b, a, true, true, true);
        if (first)
          first = false;
        else
          o << ",";
        o << b.str();
      }
    o << "]" << std::endl;
    return o;
  }

  
  ////////////////////////
  /// Append support /////
  ////////////////////////
  
  void pushBackOne(CoeffVectorType& coeffs) const
  {
    auto& svec = * coeffVector(coeffs);
    FieldElement one;
    mRing->init_set(one, 1);
    svec.emplace_back(one); // This grabs 'one' in cases where it is allocated...  I think...!
  }

  void pushBackMinusOne(CoeffVectorType& coeffs) const
  {
    auto& svec = * coeffVector(coeffs);
    FieldElement minus_one;
    mRing->init_set(minus_one, -1);
    svec.emplace_back(minus_one);
  }

  void pushBackElement(CoeffVectorType& coeffs,
                       const CoeffVectorType& take_from_here,
                       size_t loc) const
  {
    auto& svec = * coeffVector(coeffs);
    auto& svec2 = * coeffVector(take_from_here);
    assert(loc < svec2.size());
    FieldElement a;
    mRing->init_set(a, svec2[loc]);
    svec.emplace_back(a);
  }
  
  void pushBackNegatedElement(CoeffVectorType& coeffs,
                              const CoeffVectorType& take_from_here,
                              size_t loc) const
  {
    auto& svec = * coeffVector(coeffs);
    auto& svec2 = * coeffVector(take_from_here);
    assert(loc < svec2.size());
    FieldElement a;
    mRing->init_set(a, svec2[loc]);
    mRing->negate(a, a);
    svec.emplace_back(a);
  }
};

// `overloaded` construct (not standard until C++20)
template<class... Ts> 
struct overloaded : Ts... { using Ts::operator()...; };

template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

// this is the dispatching class using std::variant
class VectorArithmetic
{
  using CoeffVectorType = CoeffVector;
  using DenseCoeffVectorType = DenseCoeffVector;

  using CVA_Type = std::variant<ConcreteVectorArithmetic<M2::ARingZZpFlint>*,
                                ConcreteVectorArithmetic<M2::ARingGFFlintBig>*,
                                ConcreteVectorArithmetic<M2::ARingGFFlint>*,
                                ConcreteVectorArithmetic<M2::ARingQQGMP>*,
                                ConcreteVectorArithmetic<M2::ARingZZpFFPACK>*,
                                ConcreteVectorArithmetic<M2::ARingZZp>*,
                                ConcreteVectorArithmetic<M2::ARingGFM2>*,
                                ConcreteVectorArithmetic<M2::ARingGFGivaro>*,
                                ConcreteVectorArithmetic<CoefficientRingR>*,
                                ConcreteVectorArithmetic<M2::DummyRing>*>;

private:
  CVA_Type mConcreteVector;

public:
  VectorArithmetic(const Ring* R)
  {
    // if the rings coming in are defined in a similar way, we could avoid the switch
    // statement here as well and make it compile time.  It would make the class templated
    // on the variant class describing all the rings.
    switch (R->ringID())
      {
      case M2::ring_ZZpFlint:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingZZpFlint>* >(R)->ring());
        break;
      case M2::ring_GFFlintBig:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingGFFlintBig>* >(R)->ring());
        break;
      case M2::ring_GFFlintZech:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingGFFlint>* >(R)->ring());
        break;
      case M2::ring_QQ:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingQQGMP>* >(R)->ring());
	break;
      case M2::ring_ZZpFfpack:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingZZpFFPACK>* >(R)->ring());
	break;
      case M2::ring_ZZp:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingZZp>* >(R)->ring());
        break;
      case M2::ring_GFM2:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingGFM2>* >(R)->ring());
        break;
      case M2::ring_GFGivaro:
        mConcreteVector = new ConcreteVectorArithmetic
          (&dynamic_cast< const M2::ConcreteRing<M2::ARingGFGivaro>* >(R)->ring());
        break;
      case M2::ring_old:
        std::cout << "Using GC ring in VectorArithmetic." << std::endl;
        mConcreteVector = new ConcreteVectorArithmetic(R->getCoefficientRingR());
        break;
      default:
        // dummy ring type for default
	// throw an error here...
        std::cerr << "*** error! *** ring ID not found....!" << std::endl;
        mConcreteVector = new ConcreteVectorArithmetic<M2::DummyRing>();
      }
  }

  // provide simple visitor interface to underlying std::variant types
  size_t size(const CoeffVectorType& coeffs) const {
    return std::visit([&](auto& arg) -> size_t { return arg->size(coeffs); }, mConcreteVector);
  }
  
  size_t size(const DenseCoeffVectorType& coeffs) const {
    return std::visit([&](auto& arg) -> size_t { return arg->size(coeffs); }, mConcreteVector);
  }

  /////////////////////////////
  // Allocation/Deallocation //
  /////////////////////////////
  /// Create a coefficient vector with room for `nelems` coefficients
  CoeffVectorType allocateCoeffVector(ComponentIndex nelems) const {
    return std::visit([&](auto& arg) -> CoeffVectorType { return arg->allocateCoeffVector(nelems);}, mConcreteVector);
  }

  DenseCoeffVectorType allocateDenseCoeffVector(ComponentIndex nelems) const {
    return std::visit([&](auto& arg) -> DenseCoeffVectorType { return arg->allocateDenseCoeffVector(nelems);}, mConcreteVector);
  }

  /// Create a coefficient vector with 0 elements.  Can be increased in size.
  // CoeffVectorType allocateCoeffVector() const = 0;

  /// Deallocate all the coefficients, and the array itself.
  void deallocateCoeffVector(CoeffVectorType& coeffs) const {
    std::visit([&](auto& arg) { delete arg->coeffVector(coeffs); }, mConcreteVector);
  }
  
  void deallocateCoeffVector(DenseCoeffVectorType& coeffs) const {
    std::visit([&](auto& arg) { delete arg->denseCoeffVector(coeffs); }, mConcreteVector);
  }

  ////////////////////////
  /// Linear Algebra /////
  ////////////////////////
  void sparseRowToDenseRow(DenseCoeffVectorType& dense,
			   const CoeffVectorType& coeffs,
			   const Range<int>& comps) const {
    std::visit([&](auto& arg) { arg->sparseRowToDenseRow(dense,coeffs,comps); }, mConcreteVector);
  }

  void denseRowCancelFromSparse(DenseCoeffVectorType& dense,
                                const CoeffVectorType& coeffs,
                                const Range<int>& comps) const {
    std::visit([&](auto& arg) { arg->denseRowCancelFromSparse(dense,coeffs,comps); }, mConcreteVector);
  }

  int denseRowNextNonzero(DenseCoeffVectorType& dense,
                          int first,
                          int last) const {
    return std::visit([&](auto& arg) -> int { return arg->denseRowNextNonzero(dense,first,last); }, mConcreteVector);
  }

  void denseRowToSparseRow(DenseCoeffVectorType& dense,
                           CoeffVectorType& coeffs, // sets coeffs
                           Range<int>& comps, // sets comps
                           int first,
                           int last,
                           MemoryBlock& monomialSpace) const {
    std::visit([&](auto& arg) { arg->denseRowToSparseRow(dense,coeffs,comps,first,last,monomialSpace); }, mConcreteVector);  
  }

  template<typename LockType>
  void safeDenseRowToSparseRow(DenseCoeffVectorType& dense,
                               CoeffVectorType& coeffs, // sets coeffs
                               Range<int>& comps, // sets comps
                               int first,
                               int last,
                               MemoryBlock& monomialSpace,
                               LockType& lock) const {
    // the existence of 'template' here is a bit jarring to me, but it tells the compiler
    // that safeDenseRowToSparseRow is a template function
    std::visit([&](auto& arg) { arg->template safeDenseRowToSparseRow<LockType>(dense,coeffs,comps,first,last,monomialSpace,lock); }, mConcreteVector);      
  }
  
  void sparseRowMakeMonic(CoeffVectorType& coeffs) const {
    std::visit([&](auto& arg) { arg->sparseRowMakeMonic(coeffs); }, mConcreteVector);      
  }

  /////////////////////////////
  /// Translation    //////////
  /////////////////////////////

  template<class Container>
  void appendSparseVectorToContainer(const CoeffVectorType& coeffs,
                                     Container& c) const {
    std::visit([&](auto& arg) { arg->template appendSparseVectorToContainer<Container>(coeffs, c); }, mConcreteVector);      
  }

  template<class Container>
  CoeffVectorType sparseVectorFromContainer(const Container& c) const {
    return std::visit([&](auto& arg) -> CoeffVectorType { return arg->template sparseVectorFromContainer<Container>(c); }, mConcreteVector);
  }

  ring_elem ringElemFromSparseVector(const CoeffVectorType& coeffs,
                                     int index) const {
    return std::visit([&](auto& arg) -> ring_elem { return arg->ringElemFromSparseVector(coeffs,index); }, mConcreteVector);
  }

  std::ostream& displayVector(std::ostream& o, const CoeffVectorType& v) const
  {
    return std::visit([&](auto& arg) -> std::ostream& { return arg->displayVector(o, v); }, mConcreteVector);
  }
};

// don't need this right now, but may again later.
#if 0
template<typename T>
struct GCVector : public our_new_delete
{
  template<typename... Args>
  GCVector(Args&&... args) : mVector(args...) {}

  VECTOR(T)* value() { return &mVector; }

  VECTOR(T) mVector;
};

template<> 
class ConcreteVectorArithmetic<Ring>
{
  using FieldElement = ring_elem;
  using DenseFieldElement = FieldElement;

  friend class VectorArithmetic;

  ConcreteVectorArithmetic() : mRing(nullptr) {}

private:
  const Ring* mRing;

  VECTOR(FieldElement)* coeffVector(CoeffVector f) const
  {
    return (reinterpret_cast<GCVector<FieldElement>*>(f.mValue))->value();
  }

  VECTOR(DenseFieldElement)* denseCoeffVector(DenseCoeffVector f) const
  {
    return (reinterpret_cast<GCVector<DenseFieldElement>*>(f.mValue))->value();
  }
  
  CoeffVector coeffVector(GCVector<FieldElement>* vals) const
  {
    CoeffVector result;
    result.mValue = vals;
    vals = nullptr;
    return result;
  }

  DenseCoeffVector denseCoeffVector(GCVector<DenseFieldElement>* vals) const
  {
    DenseCoeffVector result;
    result.mValue = vals;
    vals = nullptr;
    return result;
  }

public:
  ConcreteVectorArithmetic(const Ring* R) : mRing(R) {}
  
  /////////////////////////////
  // Allocation/Deallocation //
  /////////////////////////////
  CoeffVector allocateCoeffVector(ComponentIndex nelems) const
  {
    auto tempPtr = new GCVector<FieldElement>(nelems);
    for (auto& d : *tempPtr->value())
      {
	d = mRing->zero();
      }
    return coeffVector(tempPtr);    
  }

  DenseCoeffVector allocateDenseCoeffVector(ComponentIndex nelems) const
  {
    auto tempPtr = new GCVector<DenseFieldElement>(nelems);
    for (auto& d : *tempPtr->value())
      {
        d = mRing->zero();
      }
    return denseCoeffVector(tempPtr);    
  }
  // CoeffVector allocateCoeffVector() const
  // {
  //   return coeffVector(new VECTOR(FieldElement));
  // }
  void deallocateCoeffVector(CoeffVector& coeffs) const
  {
    delete coeffVector(coeffs);
  }

  void deallocateCoeffVector(DenseCoeffVector& coeffs) const
  {
    delete denseCoeffVector(coeffs);
  }
  
  ////////////////////////
  /// Linear Algebra /////
  ////////////////////////

  size_t size(const CoeffVector& coeffs) const
  {
    return coeffVector(coeffs)->size();
  }

  size_t size(const DenseCoeffVector& coeffs) const
  {
    return denseCoeffVector(coeffs)->size();
  }
  
  // better name: fillDenseRow.
  void sparseRowToDenseRow(DenseCoeffVector& dense,
                           const CoeffVector& sparse,
                           const Range<int>& comps) const
  {
    // Note: this function simply fills in the values coming from '(sparse, comps)'.
    //   Other values are not touched.
    // In our intended uses, the input `dense` is the vector consisting of all zeros`.
    
    auto& dvec = * denseCoeffVector(dense);
    auto& svec = * coeffVector(sparse);
    
    assert(comps.size() == svec.size());
    assert(comps[comps.size()-1] < dense.size());
    
    auto len = comps.size();
    for (ComponentIndex i = 0; i < len; i++) 
      dvec[comps[i]] = svec[i];
  }

  void denseRowCancelFromSparse(DenseCoeffVector& dense,
                                const CoeffVector& sparse,
                                const Range<int>& comps) const
  {
    // ASSUMPTION: svec[0] == 1.
    auto& dvec = * denseCoeffVector(dense);
    auto& svec = * coeffVector(sparse);

    ring_elem a = dvec[comps[0]];
    for (int i=0; i < comps.size(); ++i)
    {
      ring_elem b = mRing->mult(a,svec[i]);
      dvec[comps[i]] = mRing->subtract(dvec[comps[i]],b);
    }
  }

  int denseRowNextNonzero(DenseCoeffVector& dense,
                          int first,
                          int last) const
  {
    auto& dvec = * denseCoeffVector(dense);
    for (int i = first; i <= last; i++)
      {
	if (mRing->is_zero(dvec[i])) continue;
	// these lines give the gist of how to handle delayed modulus
	//if (dvec[i] > mCharacteristic) dvec[i] %= mCharacteristic;
	//else if (dvec[i] < 0) dvec[i] %= mCharacteristic;
	//if (dvec[i] == 0) continue;
	return i;
      }
    return last + 1;
  }

  void denseRowToSparseRow(DenseCoeffVector& dense,
                           CoeffVector& sparse, // output value: sets this value
                           Range<int>& comps, // output value: sets comps
                           int first,
                           int last,
                           MemoryBlock& monomialSpace) const
  { 
    tbb::null_mutex noLock;
    safeDenseRowToSparseRow<tbb::null_mutex>(dense,
					     sparse,
					     comps,
					     first,
					     last,
					     monomialSpace,
					     noLock);
  }

  template<typename LockType>
  void safeDenseRowToSparseRow(DenseCoeffVector& dense,
                               CoeffVector& sparse, // output value: sets this value
                               Range<int>& comps, // output value: sets comps
                               int first,
                               int last,
                               MemoryBlock& monomialSpace,
                               LockType& lock) const
  {
    auto& dvec = * denseCoeffVector(dense);
    
    int len = 0;
    
    // first can be -1 if the row is zero.  in this case, we should
    // not be accessing dense[i] for i negative.
    for (int i = first; i >= 0 and i <= last; i++)
      {
	if (not mRing->is_zero(dvec[i])) len++;
      }

    comps = monomialSpace.safeAllocateArray<int,LockType>(len,lock);
    sparse = allocateCoeffVector(len);
    auto& svec = * coeffVector(sparse);

    int next = 0;
    for (int i = first; i >= 0 and i <= last; i++)
      if (not mRing->is_zero(dvec[i]))
	{
          svec[next] = dvec[i];
	  comps[next] = i;
	  ++next;
	  dvec[i] = mRing->zero();
	}
  }

  void sparseRowMakeMonic(CoeffVector& sparse) const
  {
    auto& svec = * coeffVector(sparse);

    FieldElement leadCoeffInv = mRing->invert(svec[0]);

    for (auto& c : svec) { c = mRing->mult(c, leadCoeffInv); }
  }
  
  template<typename Container>
  void appendSparseVectorToContainer(const CoeffVector& sparse,
                                     Container& c) const
  {
    auto& svec = * coeffVector(sparse);
    for (auto a : svec)
      {
	c.push_back(a);
      }
  }

  template<typename Container>
  CoeffVector sparseVectorFromContainer(const Container& c) const
  {
    CoeffVector sparse = allocateCoeffVector(c.size());
    auto& svec = * coeffVector(sparse);
    for (auto i = 0; i < c.size(); ++i)
      {
        svec[i] = c[i];
      }
    return sparse;
  }

  ring_elem ringElemFromSparseVector(const CoeffVector& sparse,
                                     int index) const
  {
    auto& svec = * coeffVector(sparse);
    return svec[index];
  }

  std::ostream& displayVector(std::ostream& o, const CoeffVector& v) const
  {
    auto& svec = * coeffVector(v);

    bool first = true;
    o << "[(" << svec.size() << ")";
    for (const auto& a : svec)
      {
        buffer b;
        mRing->elem_text_out(b, a, true, true, true);
        if (first)
          first = false;
        else
          o << ",";
        o << b.str();
      }
    o << "]" << std::endl;
    return o;
  }
  
};

#endif

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e  "
// indent-tabs-mode: nil
// End: