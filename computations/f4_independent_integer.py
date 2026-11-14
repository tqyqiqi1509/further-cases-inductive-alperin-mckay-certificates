"""Independent finite F4 checks used in the odd-prime audit.

The roots are built in Euclidean coordinates, rather than imported from the
previous certificates.  All arithmetic is rational or integral.  This is a
certificate of the stated root/normalizer assertions only, not of an ordinary
character bijection or of the inductive Alperin--McKay condition.
"""
from fractions import Fraction as F
from itertools import product, combinations
from collections import deque
from pathlib import Path
import json

def dot(x, y): return sum(a*b for a,b in zip(x,y))
def mul(a,b): return tuple(tuple(dot(row,col) for col in zip(*b)) for row in a)
def add(a,b): return tuple(tuple(x+y for x,y in zip(r,s)) for r,s in zip(a,b))
def scale(k,a): return tuple(tuple(k*x for x in r) for r in a)
def sub(a,b): return add(a,scale(-1,b))
def vec(a,v): return tuple(dot(r,v) for r in a)
def ident(n): return tuple(tuple(int(i==j) for j in range(n)) for i in range(n))
I=ident(4)
def inverse(a):
    m=[list(map(F,r))+list(map(F,s)) for r,s in zip(a,I)]
    for j in range(4):
        k=next(i for i in range(j,4) if m[i][j]);m[j],m[k]=m[k],m[j]
        d=m[j][j];m[j]=[x/d for x in m[j]]
        for i in range(4):
            if i!=j:
                d=m[i][j];m[i]=[x-d*y for x,y in zip(m[i],m[j])]
    return tuple(tuple(r[4:]) for r in m)
def integral(a):
    assert all(F(x).denominator==1 for r in a for x in r)
    return tuple(tuple(int(x) for x in r) for r in a)
def reflect(a):
    return tuple(tuple(F(int(i==j))-F(2*a[i]*a[j],dot(a,a)) for j in range(4)) for i in range(4))

# Twice the usual Euclidean roots.
euclidean=set()
for i in range(4):
    for z in [-2,2]:
        x=[0]*4;x[i]=z;euclidean.add(tuple(x))
for i,j in combinations(range(4),2):
    for z,w in product([-2,2],repeat=2):
        x=[0]*4;x[i]=z;x[j]=w;euclidean.add(tuple(x))
euclidean.update(product([-1,1],repeat=4))
simple=((0,2,-2,0),(0,0,2,-2),(0,0,0,2),(1,-1,-1,-1))
P=tuple(zip(*simple));Pi=inverse(P)
roots={tuple(int(x) for x in vec(Pi,r)) for r in euclidean}
S=[integral(mul(mul(Pi,reflect(r)),P)) for r in simple]
W={I};queue=deque([I])
while queue:
    h=queue.popleft()
    for s in S:
        g=mul(s,h)
        if g not in W: W.add(g);queue.append(g)
assert len(roots)==48 and len(W)==1152
assert all(vec(s,r) in roots for s in S for r in roots)
assert all(F(x).denominator==1 for r in euclidean for x in vec(Pi,r))

def zero_mod(a,m): return all(x % m==0 for r in a for x in r)
def vzero_mod(x,m): return all(y % m==0 for y in x)
def max_abs(a): return max(abs(x) for r in a for x in r)
mixed=[]
for i,j in [(0,1),(2,3)]:
    w=mul(S[i],S[j]);w2=mul(w,w)
    assert mul(w2,w)==I
    N=add(add(I,w),w2);K=add(w,scale(2,w2))
    low=[]
    for a,u in product([1,2],[1,2]):
        Q=add(N,scale(3**a*u,K));mod=3**(a+1)
        r0=sum(vzero_mod(vec(Q,r),mod) for r in roots)
        w0=sum(zero_mod(mul(Q,sub(g,I)),mod) for g in W)
        assert (r0,w0)==(0,1)
        low.append([a,u,r0,w0])
    rb=max(abs(x) for r in roots for x in vec(N,r))
    wb=max(max_abs(mul(N,sub(g,I))) for g in W)
    assert max(rb,wb)<27
    assert all(not vzero_mod(vec(K,r),3) for r in roots if not any(vec(N,r)))
    assert sum(zero_mod(mul(K,sub(g,I)),3) for g in W if max_abs(mul(N,sub(g,I)))==0)==1
    perp=[r for r in euclidean if dot(r,simple[i])==dot(r,simple[j])==0]
    assert len(perp)==6
    b,c=next((b,c) for b in perp for c in perp if dot(b,b)==dot(c,c) and 2*dot(b,c)==-dot(b,b))
    v=integral(mul(mul(Pi,mul(reflect(b),reflect(c))),P))
    delta=mul(sub(v,I),sub(v,I))
    assert mul(v,w)==mul(w,v) and mul(mul(v,v),v)==I
    assert all(not zero_mod(mul(add(N,scale(3*u,K)),delta),9) for u in [1,2])
    ND=mul(N,delta)
    assert any(x%3==0 and x%9!=0 for r in ND for x in r)
    mixed.append({'word':f's{i+1}s{j+1}','low_tests_a_u_roots_weyl':low,
      'N_root_bound':rb,'N_weyl_bound':wb,'perpendicular_roots':len(perp),
      'perpendicular_A2_pair_euclidean_scaled2':[b,c], 'N_delta':ND})

assert len({tuple(tuple(x%3 for x in r) for r in g) for g in W})==1152
order3=[g for g in W if g!=I and mul(mul(g,g),g)==I]
# For order three over Q, rank(g-1)=2 iff trace(g)=1.
rank2=[g for g in order3 if sum(g[i][i] for i in range(4))==1]
assert len(rank2)==64
assert all(not zero_mod(mul(sub(g,I),sub(g,I)),3) for g in rank2)
remaining=set(rank2);classes=[]
while remaining:
    g=next(iter(remaining));orb={integral(mul(mul(h,g),inverse(h))) for h in W}
    classes.append(len(orb));remaining-=orb
assert sorted(classes)==[32,32]
rank1=[]
for j in [0,2]:
    alpha=tuple(I[i][j] for i in range(4))
    line={tuple(z*x%3 for x in alpha) for z in range(3)}
    rz=[r for r in roots if tuple(x%3 for x in r) in line]
    wz=[g for g in W if all(tuple(x%3 for x in col) in line for col in zip(*sub(g,I)))]
    assert set(rz)=={alpha,tuple(-x for x in alpha)} and set(wz)=={I,S[j]}
    rank1.append({'simple_root':j+1,'root_annihilator_count':len(rz),'pointwise_W_stabilizer':len(wz)})

# The connected centralizer roots of a B2 subsystem on e3,e4.
H=[r for r in euclidean if r[0]==r[1]==0]
orthogonal=[r for r in euclidean if all(dot(r,h)==0 for h in H)]
odd_centralizer=[r for r in orthogonal if all(tuple(x+y for x,y in zip(r,h)) not in euclidean for h in H)]
assert len(H)==8 and len(orthogonal)==8 and len(odd_centralizer)==4
result={'root_count':48,'W_order':1152,'mixed':mixed,'pure_torus_mod3_faithful':True,
  'order3_rank2_class_sizes':sorted(classes),'rank1':rank1,
  'B2_source_root_count':len(H),'B2_orthogonal_root_count':len(orthogonal),
  'odd_characteristic_centralizer_roots':sorted(odd_centralizer),
  'scope':'Exact finite root checks; character exhaustion, blocks and outer covariance are separate proofs.'}
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
