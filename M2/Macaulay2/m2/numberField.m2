needs "enginering.m2"
needs "quotring.m2"
needs "polyrings.m2"

NumberField = new Type of EngineRing
NumberField.synonym = "Number field"

numberField = method (
     Options => { 
	  -- PrimitiveElement => FindOne,
	  Variable => null    -- use the variable in the corresponding ring QQ[x]/(f)
	  }
     )

unpack = (S) -> (			  -- a quotient ring
     if class S =!= QuotientRing
     then error("expected ",toString S," to be a quotient ring");
     R := ambient S;
     if class R =!= PolynomialRing
     or numgens R =!= 1
     then error("expected ",toString R," to be a polynomial ring in one variable");
     if degreeLength R =!= 1
     then error("expected ",toString R," to be a simply graded polynomial ring");
     if coefficientRing R =!= QQ
     then error("expected ",toString R," to be a polynomial ring over QQ");
     I := flatten entries presentation S;
     if #I =!= 1
     then error("expected ",toString S," to be a quotient ring by a principal ideal");
     f := I_0;
     n := first degree f;
     (R,n,f)
)

numberField(Ring) := NumberField => opts -> (S) -> (
     (R,n,f) := unpack S;
     if not isPrime f
     then error("expected ",toString S," to be a quotient ring by an irreducible polynomial");
     var := if opts.Variable === null then S.generatorSymbols#0 else baseName opts.Variable;

     rawF := rawNumberField raw f;
     F := new NumberField from rawF;
     F.degreeLength = 0;
     F.rawNumberField = true;

     F.isBasic = true;

     toString F := h -> toString expression h;
     net F := h -> net expression h;

     F.baseRings = append(S.baseRings,S);

     F.promoteDegree = makepromoter 0;			    -- do this before commonEngineRingInitializations
     F.liftDegree = makepromoter degreeLength S;	    -- do this before commonEngineRingInitializations
     commonEngineRingInitializations F;

     F.isCommutative = true;
     expression F := t -> expression lift(t, S);

     F.char = 0;
     F.frac = F;
     F.generators = apply(generators S, m -> promote(m,F)); -- this will be wrong if S is a tower
     if S.?generatorSymbols then F.generatorSymbols = S.generatorSymbols;
     if S.?generatorExpressions then F.generatorExpressions = S.generatorExpressions;
     if S.?generators then F.generators = apply(S.generators, r -> promote(r,F));
     if S.?indexSymbols then F.indexSymbols = applyValues(S.indexSymbols, r -> promote(r,F));
     if S.?indexStrings then F.indexStrings = applyValues(S.indexStrings, r -> promote(r,F));
     baseName F := y -> (
	  if F_0 == y then var else error "expected a generator"
	  );
     F.use = F -> var <- F_0;
     F.use F;
     F / F := (x,y) -> if y == 0 then error "division by zero" else x // y;
     F % F := (x,y) -> if y == 0 then x else 0_F;
     F
)
