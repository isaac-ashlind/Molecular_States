# Independent GAP cross-check of the group-theoretic claims used in the figures.
# Run:  gap -q -b --quitonbreak checks/verify.g   (prints a report; a failed Assert exits nonzero)
# Orders and cosets; [H,B] (10 subgroups, 17 covers) and [H,S] (36, 73); the three order-12 extensions of G6; the
# class sizes of G6; the proton-spin weights by parity; C[G6/H] = A1 + E and the sign of H induced to A2 + E; the KRb
# groups.
#
# Model: the five protons are points 1..5; inversion E* is modeled as the
# transposition (8,9) of two extra points, so S = S_5 x C_2 is a permutation
# group of order 240.  A starred operation sigma* is sigma*(8,9).

b := (2,3)(4,5)(8,9);;  t := (1,2,3);;  u := (4,5);;  Estar := (8,9);;
S := Group((1,2),(1,2,3,4,5),Estar);;
H := Group(b);;  G6 := Group(b,t);;  G12 := Group(b,t,u);;  B := Group(b,t,u,Estar);;

Print("orders H,G6,G12,B,S = ", List([H,G6,G12,B,S],Size), "\n");
Assert(0, List([H,G6,G12,B,S],Size) = [2,6,12,24,240]);
Assert(0, IsSubgroup(G6,H) and IsSubgroup(G12,G6) and IsSubgroup(B,G12) and IsSubgroup(S,B));
Assert(0, G12 = Group(b,t,t*u));            # G12 = <G6, tu>
Assert(0, b*t*b = t^-1 and b*t <> t*b and b*t = t^2*b);
Print("relations btb=t^-1, bt=t^2b<>tb: ok\n");

Print("cosets |G6/H|,|G12/H|,|S/G12|,|S/H| = ",
      [Index(G6,H),Index(G12,H),Index(S,G12),Index(S,H)], "\n");
Assert(0, [Index(G6,H),Index(G12,H),Index(S,G12),Index(S,H)] = [3,6,20,120]);

# bond-preserving interval [H,B]: 10 subgroups, 17 covering relations
intB := IntermediateSubgroups(B,H);;
nB := Length(intB.subgroups)+2;;
Print("interval [H,B]: subgroups = ", nB, ", covers = ", Length(intB.inclusions), "\n");
Assert(0, nB = 10 and Length(intB.inclusions) = 17);
Print("orders in [H,B]: ", SortedList(Concatenation([2,24],List(intB.subgroups,Size))), "\n");
# the three order-12 candidates above G6 and that any two generate B
cand := Filtered(intB.subgroups, K -> Size(K)=12 and IsSubgroup(K,G6));;
Assert(0, Length(cand) = 3);
Assert(0, ForAll(Combinations(cand,2), p -> ClosureGroup(p[1],p[2]) = B));
Print("three order-12 extensions of G6; any two generate B: ok\n");

# full interval [H,S]: manuscript states 36 subgroups
intS := IntermediateSubgroups(S,H);;
Print("interval [H,S]: subgroups = ", Length(intS.subgroups)+2, ", covers = ", Length(intS.inclusions), "\n");
Assert(0, Length(intS.subgroups)+2 = 36 and Length(intS.inclusions) = 73);

# character table of G6 (= S_3 as abstract group) and the five-proton spin weights
tbl := CharacterTable(G6);;
irr := Irr(tbl);;
cls := ConjugacyClasses(tbl);;
reps := List(cls, Representative);;
Assert(0, SortedList(List(cls, Size)) = [1, 2, 3]);
# class functions: spin character 2^(cycles on 1..5), statistics sgn, parity (+/- on E*)
cyclesOn5 := g -> Length(Cycles(g,[1..5]));;
spinchar := List(reps, g -> 2^cyclesOn5(g));;
statchar := List(reps, g -> SignPerm(RestrictedPerm(g,[1..5])));;
parchar  := List(reps, g -> SignPerm(RestrictedPerm(g,[8,9])));;
# name the irreducibles by their values on t and b
nameOf := function(chi)
  local vb;
  vb := chi[Position(reps, First(reps, g -> g in ConjugacyClass(G6,b)))];
  if chi[1] = 2 then return "E"; elif vb = 1 then return "A1"; else return "A2"; fi;
end;;
weights := function(parity)
  local w, chi, cf, i;
  w := [];
  for chi in irr do
    cf := List([1..Length(reps)], i -> chi[i]*spinchar[i]*statchar[i]*parchar[i]^((1-parity)/2));
    Add(w, [nameOf(chi), ScalarProduct(tbl, ClassFunction(tbl,cf), TrivialCharacter(tbl))]);
  od;
  return w;
end;;
we := weights(1);; wo := weights(-1);;
Print("spin weights even parity: ", we, "\n");
Print("spin weights odd parity:  ", wo, "\n");
lookup := function(w,n) return First(w, p -> p[1]=n)[2]; end;;
Assert(0, [lookup(we,"A1"),lookup(we,"A2"),lookup(we,"E")] = [12,4,8]);
Assert(0, [lookup(wo,"A1"),lookup(wo,"A2"),lookup(wo,"E")] = [4,12,8]);
Assert(0, lookup(we,"A1")+lookup(we,"A2")+2*lookup(we,"E") = 32);
# the full spin space: carbon-12 (C^1) and nitrogen-14 (C^3) are not permuted, so every weight triples
Assert(0, List(["A1","A2","E"], n -> 1*3*lookup(we,n)) = [36,12,24] and List(["A1","A2","E"], n -> 1*3*lookup(wo,n)) = [12,36,24]);
# permutation character of G6 on G6/H: multiplicities (A1,A2,E) = (1,0,1)
permchar := PermutationCharacter(G6,H);;
mult := List(irr, chi -> [nameOf(chi), ScalarProduct(tbl, ClassFunction(tbl, List(reps, g -> permchar[Position(reps,g)])), chi)]);;
Print("C[G6/H] multiplicities: ", mult, "\n");
Assert(0, [lookup(mult,"A1"),lookup(mult,"A2"),lookup(mult,"E")] = [1,0,1]);
# induced from the sign of H: A2 + E, so C[G6] = C[G6/H] + A2 + E
sgnH := First(Irr(H), chi -> chi[1] = 1 and chi <> TrivialCharacter(H));;
indsgn := InducedClassFunction(sgnH, G6);;
msgn := List(irr, chi -> [nameOf(chi), ScalarProduct(tbl, chi, indsgn)]);;
Assert(0, [lookup(msgn,"A1"),lookup(msgn,"A2"),lookup(msgn,"E")] = [0,1,1]);
Assert(0, List(irr, chi -> ScalarProduct(tbl, chi, PermutationCharacter(G6, TrivialSubgroup(G6)))) =
          List(irr, chi -> ScalarProduct(tbl, chi, permchar) + ScalarProduct(tbl, chi, indsgn)));
Print("induced from the sign of H: ", msgn, "\n");

# KRb + KRb: S = <(1,2),(3,4),E*> = C2^3, 16 subgroups; channel groups of order 4, any two generate S
SK := Group((1,2),(3,4),(8,9));;
Gin := Group((1,2)(3,4),(8,9));;  GK := Group((1,2),(8,9));;  GRb := Group((3,4),(8,9));;
Print("KRb: |S| = ", Size(SK), ", subgroups = ", Length(AllSubgroups(SK)),
      ", channel orders = ", List([Gin,GK,GRb],Size), "\n");
Assert(0, Size(SK)=8 and Length(AllSubgroups(SK))=16 and List([Gin,GK,GRb],Size)=[4,4,4]);
Assert(0, ForAll(Combinations([Gin,GK,GRb],2), p -> ClosureGroup(p[1],p[2]) = SK));

# rigid methane in exact cyclotomic arithmetic: Td(M) is S4 on the four protons with the odd relabelings starred,
# T the even ones (A4); the species of Td(M) named by degree and value on a transposition (A2 and T1 odd there)
Td := SymmetricGroup(4);;  T := AlternatingGroup(4);;
spin4 := ClassFunction(Td, List(ConjugacyClasses(Td), c -> 2^Length(Cycles(Representative(c), [1..4]))));;
tdname := chi -> [chi[1], chi[Position(ConjugacyClasses(Td), ConjugacyClass(Td, (1,2)))]];;
Assert(0, Set(List(Irr(Td), chi -> [tdname(chi), ScalarProduct(chi, spin4)])) =
  Set([[[1,1],5], [[1,-1],0], [[2,0],1], [[3,1],3], [[3,-1],0]]));          # 5 A1 + E + 3 T2
irrT := Irr(T);;  spinT := RestrictedClassFunction(spin4, T);;
Assert(0, SortedList(List(irrT, chi -> [chi[1], ScalarProduct(chi, spinT)])) = [[1,1],[1,1],[1,5],[3,3]]);   # 5 A + 1E + 2E + 3 T
Assert(0, SortedList(List(irrT, chi -> ScalarProduct(chi, PermutationCharacter(T, TrivialSubgroup(T))))) = [1,1,1,3]);              # C[T]
Assert(0, SortedList(List(Irr(Td), chi -> SortedList(List(irrT, psi -> ScalarProduct(RestrictedClassFunction(chi, T), psi))))) =
  SortedList([[0,0,0,1],[0,0,0,1],[0,0,1,1],[0,0,0,1],[0,0,0,1]]));             # A1, A2 -> A; E -> 1E + 2E; T1, T2 -> T
isomers := Filtered(Cartesian(irrT, irrT), p -> ScalarProduct(p[1]*p[2], TrivialCharacter(T)) <> 0);;
Assert(0, SortedList(List(isomers, p -> [p[1][1], ScalarProduct(p[2], spinT)])) = [[1,1],[1,1],[1,5],[3,3]]);
Assert(0, SortedList(List(irrT, chi -> Index(T, KernelOfCharacter(chi)))) = [1,3,3,12]);                   # monodromy
Assert(0, ScalarProduct(First(irrT, chi -> chi[1] = 3)^2, TrivialCharacter(T)) = 1);   # one invariant in T x T
Print("methane: 5 A1 + E + 3 T2, 5 A + 1E + 2E + 3 T, C[T], the species on T, the isomers, |T/ker| = 1, 3, 3, 12: ok\n");
Print("GAP cross-check: all assertions passed\n");
QUIT;
