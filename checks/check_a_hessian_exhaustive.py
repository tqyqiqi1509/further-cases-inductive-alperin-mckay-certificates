"""Independent exact finite-matrix and augmentation checks for the type-A audit.

No character-triple/iAM conclusion is asserted by these computations.  This
script enumerates all 4**9 matrices, rather than checking only a proposed
normalizer's generators.  GF(4) is F_2[w]/(w*w+w+1), encoded by 0,1,2,3.
"""
from itertools import permutations, product
from json import dumps
from pathlib import Path

MUL = ((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))
INV = (0,1,3,2)
I = (1,0,0,0,1,0,0,0,1)
d = (1,0,0,0,2,0,0,0,3)
p = (0,0,1,1,0,0,0,1,0)
u = (1,0,0,0,1,0,0,0,2)
F = (1,1,1,1,2,3,1,3,2)

def mul(a,b):
    m=MUL
    return tuple(m[a[3*i]][b[j]] ^ m[a[3*i+1]][b[3+j]] ^
                 m[a[3*i+2]][b[6+j]] for i in range(3) for j in range(3))

def det(a):
    m=MUL
    return (m[a[0]][m[a[4]][a[8]] ^ m[a[5]][a[7]]] ^
            m[a[1]][m[a[3]][a[8]] ^ m[a[5]][a[6]]] ^
            m[a[2]][m[a[3]][a[7]] ^ m[a[4]][a[6]]])

def inverse(a, determinant):
    # Adjugate, with minus equal to plus in characteristic two.
    m=MUL; r=INV[determinant]
    ans=[]
    for i in range(3):
        for j in range(3):
            rows=[x for x in range(3) if x != j]
            cols=[x for x in range(3) if x != i]
            z=m[a[3*rows[0]+cols[0]]][a[3*rows[1]+cols[1]]] ^ m[a[3*rows[0]+cols[1]]][a[3*rows[1]+cols[0]]]
            ans.append(m[r][z])
    return tuple(ans)

def adjoint(a):
    return tuple(MUL[a[3*j+i]][a[3*j+i]] for i in range(3) for j in range(3))

def closure(gens):
    seen={I}; queue=[I]
    while queue:
        x=queue.pop()
        for g in gens:
            y=mul(x,g)
            if y not in seen:
                seen.add(y); queue.append(y)
    return seen

E=closure([d,p])
generated=closure([d,p,F,u])
normalizer=set(); unitary=set(); general_order=0
for a in product(range(4), repeat=9):
    t=det(a)
    if not t:
        continue
    general_order += 1
    ai=inverse(a,t)
    if mul(a,ai) != I:
        raise AssertionError("incorrect inverse")
    if mul(mul(a,d),ai) in E and mul(mul(a,p),ai) in E:
        normalizer.add(a)
    if mul(adjoint(a),a)==I:
        unitary.add(a)

assert general_order == 181440
assert len(E)==27
assert generated==normalizer==unitary
assert len(normalizer)==648
special={a for a in normalizer if det(a)==1}
assert len(special)==216
assert all(mul(mul(a,a),a)==I for a in E)
assert all(adjoint(a)==inverse(a,det(a)) for a in normalizer)

augmentation=[]
for n in range(2,8):
    for ell,aexp in ((3,1),(3,2),(5,1),(5,2)):
        modulus=ell**aexp
        kernel=[]
        for perm in permutations(range(n)):
            works=True
            # e_i-e_{n-1} generate the augmentation module.
            for i in range(n-1):
                v=[0]*n
                v[perm[i]]+=1; v[perm[n-1]]-=1
                v[i]-=1; v[n-1]+=1
                v=[z % modulus for z in v]
                if len(set(v))!=1:
                    works=False; break
            if works:
                kernel.append(perm)
        expected=3 if (n,ell,aexp)==(3,3,1) else 1
        assert len(kernel)==expected, (n,ell,aexp,kernel)
        augmentation.append({"n":n,"ell":ell,"a":aexp,"kernel_order":len(kernel)})

out={
    "field":"GF(4), w^2+w+1=0",
    "enumerated_matrices":4**9,
    "GL3_order":general_order,
    "E_order":len(E),
    "exhaustive_normalizer_order":len(normalizer),
    "exhaustive_normalizer_equals_generated_and_unitary":True,
    "special_normalizer_order":len(special),
    "field_times_graph_centralizes_normalizer":True,
    "augmentation_results":augmentation,
    "scope":"Matrix normalizers and finite augmentation kernels only; not an iAM certificate."
}
destination=Path(__file__).with_suffix('.json')
destination.write_text(dumps(out,indent=2)+"\n")
print(dumps(out,indent=2))
