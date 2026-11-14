"""Cold-round small-field check for the intrinsic A in the mixed F4 rows.

This check addresses q=4,8 only.  It constructs the character lattice and
the dual finite 3-torus directly, without importing any previous certificate.
It proves finite lattice/action assertions, not block labels or iAM.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from pathlib import Path
import json

def tr(a):
    return tuple(zip(*a))

def mm(a, b):
    return tuple(tuple(sum(x*y for x,y in zip(r,c)) for c in tr(b)) for r in a)

def mv(a, x):
    return tuple(sum(u*v for u,v in zip(r,x)) for r in a)

I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))

def inv(a):
    z=[list(map(Q,r))+list(map(Q,s)) for r,s in zip(a,I)]
    for j in range(4):
        k=next(k for k in range(j,4) if z[k][j])
        z[j],z[k]=z[k],z[j]
        v=z[j][j]
        z[j]=[x/v for x in z[j]]
        for k in range(4):
            if k!=j:
                v=z[k][j]
                z[k]=[x-v*y for x,y in zip(z[k],z[j])]
    return tuple(tuple(r[4:]) for r in z)

def refl(x):
    d=sum(z*z for z in x)
    return tuple(tuple(Q(int(i==j))-2*x[i]*x[j]/d for j in range(4)) for i in range(4))

simple=tuple(tuple(map(Q,x)) for x in (
    (0,1,-1,0),(0,0,1,-1),(0,0,0,1),
    (Q(1,2),Q(-1,2),Q(-1,2),Q(-1,2))))
B=tr(simple)
Bi=inv(B)

def lattice(a):
    z=mm(mm(Bi,a),B)
    assert all(x.denominator==1 for r in z for x in r)
    return tuple(tuple(int(x) for x in r) for r in z)

S=[lattice(refl(x)) for x in simple]
weyl={I}
frontier=[I]
while frontier:
    h=frontier.pop()
    for s in S:
        v=mm(s,h)
        if v not in weyl:
            weyl.add(v)
            frontier.append(v)
assert len(weyl)==1152

eroots=[]
for j in range(4):
    for sign in (-1,1):
        r=[Q(0)]*4
        r[j]=Q(sign)
        eroots.append(tuple(r))
for i,j in combinations(range(4),2):
    for x,y in product((-1,1),repeat=2):
        r=[Q(0)]*4
        r[i],r[j]=Q(x),Q(y)
        eroots.append(tuple(r))
eroots.extend(tuple(Q(x,2) for x in r) for r in product((-1,1),repeat=4))
roots=[mv(Bi,r) for r in eroots]
assert len(set(roots))==48
assert all(x.denominator==1 for r in roots for x in r)
roots=[tuple(map(int,r)) for r in roots]

orientations=[
    ('long_A2_twist', mm(S[0],S[1]),
     (Q(1,2),Q(-1,2),Q(-1,2),Q(-1,2)),(-1,0,0,0)),
    ('short_A2_twist', mm(S[2],S[3]),(1,1,0,0),(-1,0,-1,0)),
]
out=[]
for name,w,r1,r2 in orientations:
    u=lattice(mm(refl(tuple(map(Q,r1))),refl(tuple(map(Q,r2)))))
    assert mm(u,w)==mm(w,u) and mm(mm(u,u),u)==I
    # Characters act contragrediently on the torus points.
    U=tr(mm(u,u))
    for q,epsilon,a in [(4,1,1),(8,-1,2)]:
        modulus=3**(a+1)
        M=tr(tuple(tuple(epsilon*q*w[i][j]-I[i][j] for j in range(4)) for i in range(4)))
        A=[y for y in product(range(modulus),repeat=4)
           if all(v%modulus==0 for v in mv(M,y))]
        Aset=set(A)
        assert len(A)==3**(2*a+1)
        omega1=sum(all(3*x%modulus==0 for x in y) for y in A)
        assert omega1==9
        assert any(any(3**a*x%modulus for x in y) for y in A)
        assert all(tuple(x%modulus for x in mv(U,y)) in Aset for y in A)
        commutators={tuple((x-z)%modulus for x,z in zip(mv(U,y),y)) for y in A}
        square={tuple((x-z)%modulus for x,z in zip(mv(U,y),y)) for y in commutators}
        assert len(square)>1
        annihilating_roots=[r for r in roots if all(sum(x*z for x,z in zip(r,y))%modulus==0 for y in A)]
        assert not annihilating_roots
        # For a pointwise stabilizer the choice h versus h^{-1} is immaterial.
        pointwise=sum(all(tuple(x%modulus for x in mv(tr(h),y))==y for y in A) for h in weyl)
        assert pointwise==1
        out.append(dict(q=q,epsilon=epsilon,a=a,orientation=name,
                        A_order=len(A),A_invariant_factors=[3**a,3**(a+1)],
                        D_order=3*len(A),commutator_order=len(commutators),
                        second_commutator_image_order=len(square),
                        annihilating_roots=len(annihilating_roots),
                        pointwise_W_stabilizer=pointwise))

result={'cases':out,'scope':
        'At q=4,8 the mixed torus has the required abelian A; the Weyl lift '
        'constructs A semidirect C3, its commutator is noncentral, and the '
        'algebraic centralizer of A is exactly the torus. Character theory '
        'and block identification are separate arguments.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
