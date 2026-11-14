LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
Val63:=function(n,p) local v; v:=0; while n mod p=0 do n:=n/p;v:=v+1;od;return v;end;;
RowAction63:=function(t,u)
 local f,cp,j,x,vals,ans;
 f:=GetFusionMap(t,u); Assert(0,f<>fail);
 cp:=[];
 for j in [1..Length(f)] do
  x:=Positions(f,f[j]); Assert(0,Length(x)<=2);
  if Length(x)=1 then cp[j]:=j;else cp[j]:=First(x,a->a<>j);fi;
 od;
 vals:=List(Irr(t),ValuesOfClassFunction);
 ans:=List(vals,x->Position(vals,x{cp})); Assert(0,not fail in ans);
 Assert(0,ForAll([1..Length(ans)],j->ans[ans[j]]=j));
 return ans;
end;;
for nm in ["2.O8+(2)"] do
 t:=CharacterTable(nm);; Assert(0,t<>fail);
 zz:=ClassPositionsOfCentre(t);; vv:=List(Irr(t),ValuesOfClassFunction);;
 upper:=CharacterTable(Concatenation(nm,".2"));;
 rp:=fail;;
 if upper<>fail and GetFusionMap(t,upper)<>fail then rp:=RowAction63(t,upper);fi;
 Print("GROUP|",nm,"|size=",Size(t),"|center_orders=",OrdersClassRepresentatives(t){zz},"\n");
 for ell in [3,5,7] do
  if Size(t) mod ell=0 then
   bb:=PrimeBlocks(t,ell);; ppart:=Val63(Size(t),ell);;
   reg:=Filtered([1..Length(vv)],j->OrdersClassRepresentatives(t)[j] mod ell<>0);;
   for b in Set(bb.block) do
    if bb.defect[b]>0 then
     ids:=Positions(bb.block,b);;
     zero:=Filtered(ids,j->Val63(vv[j][1],ell)=ppart-bb.defect[b]);;
     colours:=Collected(List(zero,j->List(zz,z->vv[j][z]/vv[j][1])));;
     orb:="not computed";;
     if rp<>fail and Set(rp{zero})=Set(zero) then
      orb:=Collected(List(Orbits(Group(PermList(rp)),zero),Length));
     fi;
     Print("BLOCK|",nm,"|p=",ell,"|id=",b,"|def=",bb.defect[b],
      "|k=",Length(ids),"|k0=",Length(zero),"|l=",RankMat(List(ids,j->vv[j]{reg})),
      "|colors=",colours,"|outer_orbits=",orb,"|h0_rows=",zero,"\n");
    fi;
   od;
  fi;
 od;
od;
Print("FRESH_O8_GLOBAL_COMPLETE\n");
QUIT;
