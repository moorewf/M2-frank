newPackage(
    "FactorTower",
    Version => "0.1",
    Date => "",
    Headline => "Factoring over Towers and Adjoining Roots",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    AuxiliaryFiles => false,
    DebuggingMode => false
    )

export {"factorOverTower",
        "adjoinRoot"}

inducedPolynomialMap = method()
inducedPolynomialMap (RingMap, Ring, Ring) := (psi, tar, src) -> (
   -- psi: map between coefficient rings of src and tar
   if numgens tar != numgens src then error "Expected number of generators to match";
   if coefficientRing tar =!= target psi or coefficientRing src =!= source psi then error "Expected map between coefficient rings.";
   map(tar,src, vars tar | psi.matrix)
)

underlyingRing = method()
underlyingRing Ring := F -> (
    
)

-- TODO: update toField documentation to only have one 'toField' in a tower
--       Or: fix flattenRing to work for towers involving toField

primitiveElementRing = method(Options => {Variable => null, CoefficientRing => null})
primitiveElementRing Ring := opts -> F -> (
    -- assumes that F is toField of a flattened finite extension of QQ or ZZ/p
    (newF,phi) := flattenRing(F,CoefficientRing => opts.CoefficientRing);
    k := coefficientRing newF;
    tryPrim := sum gens newF;
    X := if opts.Variable === null then getSymbol("X") else opts.Variable;
    G := k(monoid[X]);
    psi := map(newF,G,{tryPrim});
    degExtn := numcols basis newF;
    minPoly := first flatten entries gens ker psi;
    if first degree minPoly != degExtn then error "Sum of variables not primitive.";
    -- TODO: Adapt our primitive element if the sum is not primitive
    --       allow for user-provided primitive element, or random
    quotG := G/ideal(minPoly);
    psiQuotG := map(newF,quotG,{tryPrim});
    tfQuotG := toField quotG;
    tfPsi := map(F,tfQuotG,{phi^(-1) tryPrim});
    tfPsiInv := map(tfQuotG,F,(psiQuotG^(-1)).matrix);
    assert(tfPsi * tfPsiInv == 1);
    (tfPsi,tfPsiInv)
)

factorOverTower = method()
factorOverTower RingElement := f -> (
   R := ring f;
   F := coefficientRing R;
   e := getSymbol "e";
   (psi,psiInv) := primitiveElementRing(F,Variable => e);
   G := source psi;
   y := getSymbol "y";
   S := G(monoid[y]);
   psiInd := inducedPolynomialMap(psi,R,S);
   psiInvInd := inducedPolynomialMap(psiInv,S,R);
   facList := (factor psiInvInd f) // toList / toList;
   Product apply(facList, p -> Power{psiInd(p#0),p#1})
)

adjoinRoot = method()
adjoinRoot RingElement := f -> (
   -- 
)

end--

restart
debug needsPackage "FactorTower"
check "FactorTower"

uninstallPackage "FactorTower"
restart
installPackage "FactorTower"
viewHelp "FactorTower"

restart
load "factorTower.m2"
F0 = QQ
R0 = F0[x]
f0 = x^3 - 2
F1num = F0[a]
phi1 = map(F1num,R0,{F1num_0})
F1 = toField(F1num/ideal (phi1 f0))
R1 = F1[x]
facs = factorOverTower sub(f0,R1)
f1 = first first select(1,facs // toList / toList, p -> degree(x,p#0) > 1)
F2num = F1[b]
phi2 = map(F2num,R1,{F2num_0})
F2 = toField(F2num/ideal (phi2 f1))
R2 = F2[x]
facs = factorOverTower 
sub(f0,R2)

restart
load "factorTower.m2"
F = toField(QQ[c,d]/(c^3 - 2, d^2 + d + 1))
R = F[x]
f = x^3 - 2
factorOverTower f
value factorOverTower f == f

restart
R = QQ[x]
factor (x^2 - 4)
peek toList oo

Product {Power {x-2,1}, Power{x+2,1}}

(psi,psiInv) = primitiveElementRing(F,Variable => e)
G = source psi
S = G[y]
value oo


-- we would like to:
-- 1. Flatten a field constructed as a tower
-- 2. Construct a primitive element for it
-- 3. Define the map to and from the original tower (perhaps when constructing primitive element)
-- 4. Also find change of basis matrix between two extensions
-- 5. 

restart
F = toField(QQ[c]/(c^3 - 2))
R = F[x]
factor (x^3 - 2) -- works

restart
F = toField(QQ[c]/(c^3 - 2))
G = toField(F[d]/(d^2 + c*d + c^2))
flattenRing(G,CoefficientRing => QQ)
R = G[x]
factor (x^3 - 2) -- doesn't work

restart
gbTrace = 2
F = toField(QQ[c,d]/(c^3 - 2, d^2 + d + 1))
R = F[x]
factor (x^3 - 2) -- doesn't work

S = QQ[e]
phi = map(R,S,{c + d})
f = (ker phi)_0

G = toField(S/ideal(f))
R' = G[y]
factor (y^3 - 2)

-- preimage doesn't work on elements
use R'
psi = map(F,G,{c+d})
psi e
preimage(psi,ideal c)


restart
load "factorTower.m2"
F = toField(QQ[c,d]/(c^3 - 2, d^2 + d + 1))
psi = map(F,G,{c+d})
psiInv = psi^(-1)
assert(psi(psiInv(c)) == c)
assert(psi(psiInv(d)) == d)
assert(psiInv(psi(e)) == e)
R = F[x]
f = x^3 - 2
S = G[y]
psiInvInd = inducedPolynomialMap(psiInv,S,R)
psiInd = inducedPolynomialMap(psi,R,S)
factor (psiInvInd f)

restart
load "factorTower.m2"
F = toField(QQ[c,d]/(c^3 - 2, d^2 + d + 1))
(psi,psiInv) = primitiveElementRing(F,Variable => e)
G = source psi
R = F[x]
S = G[y]
psiInd = inducedPolynomialMap(psi,R,S)
psiInvInd = inducedPolynomialMap(psiInv,S,R)
f = x^3 - 2
facList = (factor psiInvInd f) // toList / toList
product apply(facList, p -> (expression psiInd(p#0))^(p#1))
value oo

restart
load "factorTower.m2"
F = toField(QQ[c,d]/(c^3 - 2, d^2 + d + 1))
(psi,psiInv) = primitiveElementRing(F,Variable => e)
G = source psi
R = F[x]
S = G[y]
psiInd = inducedPolynomialMap(psi,R,S)
psiInvInd = inducedPolynomialMap(psiInv,S,R)
f = x^3 - 2

restart
H = QQ[c,d,e, MonomialOrder => Lex]
I = ideal (c^3 - 2, d^2 + d + 1, e^6+3*e^5+6*e^4+3*e^3+9*e+9, e - (c + d))
Igb = gb I
preimageC = c % Igb
(preimageC^3 - 2) % Igb
preimageD = d % Igb
(preimageD^2 + preimageD + 1) % Igb

-- fake matrix ring
restart
needsPackage "AssociativeAlgebras"
R = QQ <| e11,e12,e21,e22 |>
I = ideal (e11*e11 - e11, e11*e12 - e12, e11 * e21, e11*e22,
   e12*e11, e12*e12, e12 * e21 - e11, e12*e22 - e12,
   e21*e11 - e21, e21*e12 - e22, e21 * e21, e21*e22,
   e22*e11, e22*e12, e22 * e21 - e21, e22*e22 - e22)
Igb = NCGB(I, 10)
S = R/I

R' = QQ <| X |>
phi = map(S,R',{2*e11 + 4*e21 + 3*e22})
ncKernel phi   -- take a look at this example, looks wrong.

-- finite dimensional algebras
-- (commuting) matrices for variables representing multiplication map
-- can easily determine invertibility and compute inverse
