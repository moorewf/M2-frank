---------------------------------------------------------------------------
-- PURPOSE : Computation of truncations M_{>=d} for modules
--
-- UPDATE HISTORY : created Oct 2018
--                  updated Jun 2021
--
-- TODO :
-- 1. improve speed and caching
---------------------------------------------------------------------------
newPackage(
    "Truncations",
    Version => "0.7",
    Date => "31 Dec 2018",
    Headline => "truncation of a module",
    Authors => {
        { Name => "David Eisenbud", Email => "de@msri.org",           HomePage => "https://www.msri.org/~de" },
        { Name => "Mike Stillman",  Email => "mike@math.cornell.edu", HomePage => "https://www.math.cornell.edu/~mike" }
        },
    Keywords => { "Commutative Algebra" },
    PackageImports => { "Polyhedra", "NormalToricVarieties" },
    AuxiliaryFiles => true,
    DebuggingMode => true
    )

-- "truncate" is exported by Core

protect Exterior
protect Nef

--------------------------------------------------------------------
-- Helpers
--------------------------------------------------------------------

-- check whether truncation is implemented for this ring type.
truncateImplemented = method()
truncateImplemented Ring := Boolean => R -> (
    ZZ === (A := ambient first flattenRing R)
    or isAffineRing A
    or isPolynomialRing A and isAffineRing coefficientRing A and A.?SkewCommutative
    or isPolynomialRing A and ZZ === coefficientRing A
    )

-- checkOrMakeDegreeList: takes a degree, and degree rank:ZZ
-- output: a list of degrees, all of the correct length (degree rank), otherwise an error
-- in the following n represents an integer, and d represents a list of integers:
--  n          --> {{n}}          if degree rank is 1
-- {n0,...,ns} --> {{n0,...,ns}}  if length is the degree rank
-- {d0,...,ds} --> { d0,...,ds }  if length of each di is the degree rank
checkOrMakeDegreeList = method()
checkOrMakeDegreeList(ZZ,   ZZ) := (n, degrank) -> (
    if degrank === 1 then {{n}} else error("expected a degree of length " | degrank))
checkOrMakeDegreeList(List, ZZ) := (L, degrank) -> (
    if #L === 0 then error "expected nonempty list of degrees";
    if all(L, d -> instance(d, ZZ))
    then if #L === degrank then {L} else error("expected a multidegree of length " | degrank)
    else ( -- If L is a list of lists of integers, all the same length, L will be returned.
        if all(L, deg -> instance(deg, VisibleList) and #deg === degrank and all(deg, d -> instance(d, ZZ)))
        then L else error("expected a list of multidegrees, each of length " | degrank))
    )

-- Helpers for toric varieties

-- Generators of effective cone
effGenerators = method()
effGenerators Ring := (cacheValue symbol effGenerators) (R -> transpose matrix degrees R)
effGenerators NormalToricVariety := X -> effGenerators ring X

-- Nef cone of X as a polyhedral object
nefCone = method()
nefCone NormalToricVariety := (cacheValue symbol nefCone) (X -> convexHull(matrix{0_(picardGroup X)}, nefGenerators X))
nefCone Ring := R -> if R.?variety and instance(R.variety, NormalToricVariety) then nefCone R.variety

-- Effective cone of X as a polyhedral object
effCone = method()
effCone NormalToricVariety := (cacheValue symbol effCone) (X -> convexHull(matrix{0_(picardGroup X)}, effGenerators X))
effCone Ring := R -> if R.?variety and instance(R.variety, NormalToricVariety) then effCone R.variety

--------------------------------------------------------------------
-- Polyhedral algorithms
--------------------------------------------------------------------

-- If columns of A span the effective cone of a simplicial toric variety X
-- (e.g. when columns are the degrees of the variables of the Cox ring S)
-- and b is an element of the Picard group of X (e.g. a multidegree in S)
-- then truncationPolyhedron returns a polyhedron in the lattice of
-- exponent vectors of monomials of S whose degree is in b + nef cone of X
truncationPolyhedron = method(Options => { Exterior => {}, Nef => null })
-- Exterior should be a list of variable indices which are skew commutative; i.e. have max degree 1.
truncationPolyhedron(Matrix, List)   := Polyhedron => opts -> (A, b) -> truncationPolyhedron(A, transpose matrix {b}, opts)
truncationPolyhedron(Matrix, Matrix) := Polyhedron => opts -> (A, b) -> (
    -- assumption: A is m x n. b is m x 1.
    -- returns the polyhedron {Ax >= b, x_i >= 0}
    -- or if Nef is present   {Ax - b in Nef cone, x_i >= 0}
    -- TODO: why is it better to be over QQ?
    if ring A === ZZ then A = A ** QQ;
    if ring A =!= QQ then error "expected matrix over ZZ or QQ";
    -- added to ensure x_i >= 0
    I := id_(source A);
    z := map(source A, QQ^1, 0);
    -- data for the polyhedron inequalities
    hdataLHS := -(A || I);
    hdataRHS := -(b || z);
    -- added to ensure x_i <= 1 for skew commutating variables
    if #opts.Exterior > 0 then (
        -- also need to add in the conditions that each variable in the list has degree <= 1.
        hdataLHS = hdataLHS || (I ^ (opts.Exterior));
        hdataRHS = hdataRHS || matrix toList(#opts.Exterior : {1_QQ});
        );
    -- this is correct when the Nef cone equals the positive quadrant
    P := polyhedronFromHData(hdataLHS, hdataRHS);
    if opts.Nef === null or rays opts.Nef == id_(target A) then return P;
    -- otherwise, intersect with the preimage of the Nef cone
    intersection(P, affinePreimage(A, opts.Nef, -b)))

-- basisPolyhedron: this function is not used nor tested here.
-- It can be used for a perhaps better implementation of 'basis'.
-- TODO: it should have exterior variables added too... what would this mean?
basisPolyhedron = method()
basisPolyhedron(Matrix, List)   := (A, b) -> basisPolyhedron(A, transpose matrix {b})
basisPolyhedron(Matrix, Matrix) := (A, b) -> (
    -- assumption: A is m x n. b is m x 1.
    -- returns the polyhedron {Ax = b, x_i >= 0}
    if ring A === ZZ then A = A ** QQ;
    if ring A =!= QQ then error "expected matrix over ZZ or QQ";
    -- added to ensure x_i >= 0
    I := id_(source A);
    z := map(source A, QQ^1, 0);
    polyhedronFromHData(-I, -z, A, b))

--------------------------------------------------------------------
-- Algorithms for truncations of a polynomial ring
--------------------------------------------------------------------

truncationMonomials = method()
truncationMonomials(List, Ring) := (degs, R) -> (
    -- valid for total coordinate ring of any simplicial toric variety
    -- or any polynomial ring with skew commutating variables.
    degs = checkOrMakeDegreeList(degs, degreeLength R);
    -- we could call findMins on degs, but only with respect to the Nef cone
    if #degs > 1 then return sum for d in degs list truncationMonomials(d, R);
    d := degs#0;
    if  R#?(symbol truncate, d)
    then R#(symbol truncate, d)
    else R#(symbol truncate, d) = (
        (R1, phi1) := flattenRing R;
        -- generates the effective cone
        A := effGenerators R1;
        P := truncationPolyhedron(A, transpose matrix{d},
            Exterior => (options R1).SkewCommutative,
            Nef => nefCone R1);
        H := hilbertBasis cone P;
        H = for h in H list flatten entries h;
        J := ideal leadTerm ideal R1;
        ambR := ring J;
        -- generates the Nef cone
        --nefgens := for h in H list if h#0 === 0 then ambR_(drop(h, 1)) else continue;
        mongens := for h in H list if h#0 === 1 then ambR_(drop(h, 1)) else continue;
        result := mingens ideal (matrix(ambR, {mongens}) % J);
        if R1 =!= ambR then result = result ** R1;
        if R =!= R1 then result = phi1^-1 result;
        ideal result))

truncation0 = (deg, M) -> (
    -- WARNING: valid for a polynomial ring with degree length = 1.
    -- uses the engine routines for basis with Truncate => true
    -- deg: a List of integers
    -- M: Module
    -- returns a submodule of M
    if M.?generators then (
        F := cover basis(deg, deg, cokernel presentation M, Truncate => true);
        subquotient(M.generators * F, if M.?relations then M.relations))
    else image basis(deg, deg, M, Truncate => true))

truncation1 = (deg, M) -> (
    -- WARNING: valid for towers of rings with degree length = 1.
    -- uses the engine routines for basis with Truncate => true
    -- deg: a List of integers
    -- M: Module
    -- returns a submodule of M
    R := ring M;
    (R1, phi1) := flattenRing R;
    if R1 === R then return truncation0(deg, M);
    -- why not just do phi1? or M ** phi1?
    M1 := if isFreeModule M then phi1 M else subquotient(
        if M.?generators then phi1 M.generators,
        if M.?relations  then phi1 M.relations);
    result1 := truncation0(deg, M1);
    gensM := if not result1.?generators then null else phi1^-1 result1.generators;
    relnsM := if not result1.?relations then null else phi1^-1 result1.relations;
    if gensM === null and relnsM === null then phi1^-1 result1
    else subquotient(gensM, relnsM))

--------------------------------------------------------------------
-- truncate
--------------------------------------------------------------------

truncateModuleOpts := {
    MinimalGenerators => true -- whether to trim the output
    }

truncate(ZZ, Ring)   :=
truncate(ZZ, Ideal)  :=
truncate(ZZ, Matrix) :=
truncate(ZZ, Module) := truncateModuleOpts >> opts -> (d, m) -> truncate({d}, m, opts)

truncate(List, Ring) := Ideal => truncateModuleOpts >> opts -> (degs, R) -> (
    -- truncate the graded ring R in degrees >= d, for d in degs
    if not truncateImplemented R then error "cannot use truncate with this ring type";
    degs = checkOrMakeDegreeList(degs, degreeLength R);
    doTrim := if opts.MinimalGenerators then trim else identity;
    doTrim if degreeLength R === 1 and any(degrees R, d -> d =!= {0})
    then ideal truncation1(min degs, R^1)
    else ideal gens truncationMonomials(degs, R))

truncate(List, Ideal)  := Ideal  => truncateModuleOpts >> opts -> (degs, I) -> ideal truncate(degs, module I, opts)
truncate(List, Module) := Module => truncateModuleOpts >> opts -> (degs, M) -> (
    if M == 0 then return M;
    if not truncateImplemented(R := ring M) then error "cannot use truncate with this ring type";
    degs = checkOrMakeDegreeList(degs, degreeLength R);
    doTrim := if opts.MinimalGenerators then trim else identity;
    doTrim if degreeLength R === 1 and any(degrees R, d -> d =!= {0})
    then truncation1(min degs, M)
    else if isFreeModule M then (
        image map(M, , directSum for a in degrees M list
            gens truncationMonomials(for d in degs list (d - a), R)))
    else if M.?generators and not M.?relations then (
        image map(M, , gens truncate(degs, target presentation M, opts)))
    else subquotient(
        gens truncate(degs, image generators M, MinimalGenerators => false),
        gens truncate(degs, image  relations M, MinimalGenerators => false))
    )

truncate(List, Matrix) := Matrix => truncateModuleOpts >> opts -> (degs, f) -> (
    if not truncateImplemented(R := ring f) then error "cannot use truncate with this ring type";
    degs = checkOrMakeDegreeList(degs, degreeLength R);
    F := truncate(degs, source f, opts);
    G := truncate(degs, target f, opts);
    map(G, F, (f * gens F) // gens G))

--------------------------------------------------------------------
----- Tests section
--------------------------------------------------------------------

load "./Truncations/tests.m2"

--------------------------------------------------------------------
----- Documentation section
--------------------------------------------------------------------

beginDocumentation()
load "./Truncations/docs.m2"

--------------------------------------------------------------------
----- Development section
--------------------------------------------------------------------

end--

restart
uninstallPackage "Truncations"
restart
loadPackage "Truncations"
debug needsPackage "Truncations"
restart
installPackage "Truncations"
check "Truncations"
