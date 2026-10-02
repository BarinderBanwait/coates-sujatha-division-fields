// Finite subfields in the division towers of Coates--Sujatha Examples 4.7--4.8.
// Run: magma -n -b case_id:=5a descent.m
// Optional: grh:=1 (conditional upper bounds), effort:=2, twist:=-7.
SetSeed(20261002);
SetNthreads(1);
SetMemoryLimit(4*1024^3);
if not assigned case_id then case_id := "5a"; end if;
if not assigned grh then grh := "0"; end if;
if not assigned effort then effort := "2"; end if;
if not assigned twist then twist := "1"; end if;
if not assigned diagnostic then diagnostic := "0"; end if;
if not assigned local_hints then local_hints := "0"; end if;
if not assigned selmer_only then selmer_only := "0"; end if;
if diagnostic eq "1" then
    SetVerbose("TwoDescent", 1);
    SetVerbose("ClassGroup", 1);
    SetVerbose("Selmer", 2);
end if;
R<x> := PolynomialRing(Rationals());
polys := AssociativeArray();
polys["5a"] := x^5-x^4+4*x^3+4*x^2-x+13;
polys["5b"] := x^5+5*x^3+5*x-1;
polys["5c"] := x^10-18*x^5-324;
polys["5d"] := x^10-54*x^5+324;
polys["7a"] := x^7-18;
polys["7b"] := x^7-14*x^4+7*x^3+28*x^2+21*x+8;
assert IsDefined(polys, case_id);
if case_id[1] eq "5" then
    aq := [1,0,0,-3,-3];
else
    aq := [1,0,0,-1,-1];
end if;
K<a> := NumberField(polys[case_id]);
E0 := EllipticCurve([Rationals() | c : c in aq]);
if StringToInteger(twist) ne 1 then
    E0 := QuadraticTwist(E0, StringToInteger(twist));
end if;
E := BaseChange(E0, K);
v1,v2,v3 := GetVersion();
printf "VERSION %o.%o-%o\n", v1,v2,v3;
printf "CASE %o TWIST %o GRH %o\n", case_id, twist, grh;
printf "FIELD %o\nCURVE %o\n", DefiningPolynomial(K), aInvariants(E);
if grh eq "1" then
    SetClassGroupBounds("GRH");
    print "WARNING: rank upper bounds assume GRH for class-group computations.";
else
    print "Class groups: default unconditional bounds; clean startup (-n).";
end if;
t := Cputime();
printf "DISCRIMINANT %o\n", Discriminant(Integers(K));
printf "ROOT_NUMBER %o\n", RootNumber(E);
printf "TORSION_BOUND %o\n", TorsionBound(E,10);
if diagnostic eq "1" or local_hints eq "1" or selmer_only eq "1" then
    hints := {K|};
    if local_hints eq "1" then
        assert case_id eq "7a" and twist eq "1";
        // x(t) for the minimal model, truncated at t^20. Here t=-x/y.
        // a is a uniformizer at the prime over 2. Multiplication by 4
        // matches the integral cubic X^3+X^2-16X-64 used in TwoSelmerGroup.
        // Hints are suggestions only: Magma verifies each local point.
        cc := [1,1,2,3,3,2,-3,-15,-37,-72,-113,-133,-65,228,965,2413,4709,7370,8295];
        fx := func< t | t^-2-t^-1 + &+[cc[i]*t^(i+1):i in [1..#cc]] >;
        hints := {4*fx(s*a^k): k in [1..14], s in [-1,1]};
        printf "FORMAL_LOCAL_HINTS %o\n", #hints;
    end if;
    ss := TwoSelmerGroup(E : Hints:=hints);
    printf "TWO_SELMER_INVARIANTS %o\n", Invariants(ss);
end if;
if selmer_only eq "1" then
    printf "CPU_SECONDS %o\nMEMORY_BYTES %o\n", Cputime(t), GetMaximumMemoryUsage();
    print "COMPLETED: Selmer computation only; no claim of a rank lower bound.";
    quit;
end if;
r, pts, sha := MordellWeilShaInformation(E : Effort:=StringToInteger(effort));
printf "RANK_BOUNDS %o\n", r;
printf "POINTS %o\n", pts;
printf "SHA_INFORMATION %o\n", sha;
for P in pts do
    printf "POINT_ORDER %o\n", Order(P);
end for;
printf "CPU_SECONDS %o\nMEMORY_BYTES %o\n", Cputime(t), GetMaximumMemoryUsage();
print "COMPLETED";
quit;
