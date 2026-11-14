LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
roots:=[];;
for pair in Combinations([1..8],2) do
 for u in [-2,2] do for v in [-2,2] do
  rr:=List([1..8],i->0);;rr[pair[1]]:=u;;rr[pair[2]]:=v;;Add(roots,rr);
 od;od;
od;
for signs in Tuples([-1,1],8) do
 if Product(signs)=1 then Add(roots,signs);fi;
od;
Sort(roots);;Assert(0,Length(Set(roots))=240);
simple:=[[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],
 [0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],
 [0,0,0,0,-2,2,0,0],[0,0,0,0,0,-2,2,0]];;
refl:=List(simple,a->PermList(List(roots,v->Position(roots,v-(v*a)/4*a))));;
wy:=Group(refl);;Assert(0,Size(wy)=696729600);
wg:=DerivedSubgroup(wy);;Assert(0,Size(wg)=348364800 and IsPerfectGroup(wg));
zc:=Centre(wg);;Assert(0,Size(zc)=2);;zi:=First(Elements(zc),x->Order(x)=2);;
Print("WEYL|degree=",NrMovedPoints(wy),"|order=",Size(wy),"|derived=",Size(wg),"|centre=",Size(zc),"\n");
Rows2W:=function(lo,hi)
 local fu,cp,j,is,vals;
 fu:=FusionConjugacyClasses(lo,hi);
 Assert(0,fu<>fail);
 cp:=[];
 for j in [1..Length(fu)] do
  is:=Positions(fu,fu[j]);Assert(0,Length(is)<=2);
  if Length(is)=1 then cp[j]:=j;else cp[j]:=First(is,x->x<>j);fi;
 od;
 vals:=List(Irr(lo),ValuesOfClassFunction);
 return List(vals,v->Position(vals,v{cp}));
end;;
for prime in [3,5] do
 pp:=SylowSubgroup(wg,prime);;nn:=Normalizer(wg,pp);;hh:=Normalizer(wy,pp);;
 cc:=Centralizer(wg,pp);;
 Print("ACTUAL_NORMALIZER|p=",prime,"|P=",Size(pp),"|N=",Size(nn),"|H=",Size(hh),"|C=",Size(cc),"|ZP=",Size(Centre(pp)),"|NmodP=",[Size(nn/pp),Exponent(nn/pp),Size(Centre(nn/pp)),Size(DerivedSubgroup(nn/pp)),AbelianInvariants(nn/pp)],"\n");
 Assert(0,Size(cc)=2*Size(Centre(pp)) and cc=ClosureGroup(Centre(pp),zc));
 ct:=CharacterTable(nn);;ht:=CharacterTable(hh);;rp:=Rows2W(ct,ht);;
 zpos:=PositionProperty(ConjugacyClasses(nn),c->zi in c);;
 ids:=Filtered([1..Length(Irr(ct))],j->Irr(ct)[j][1] mod prime<>0 and Irr(ct)[j][zpos]=-Irr(ct)[j][1]);;
 bl:=PrimeBlocks(ct,prime);;
 Assert(0,Length(Set(ids,j->bl.block[j]))=1);
 Print("FAITHFUL_LOCAL|p=",prime,"|rows=",ids,"|degree_collection=",Collected(List(ids,j->Irr(ct)[j][1])),
 "|outer_orbits=",Collected(List(Orbits(Group(PermList(rp)),ids),Length)),"|images=",rp{ids},"\n");
od;
Print("FRESH_O8_COMPLETE\n");
QUIT;
