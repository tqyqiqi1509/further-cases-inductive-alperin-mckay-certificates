#!/usr/bin/env python3
"""Exact checks for the quartic D4 source and the E8 mod-2 cover model."""
from pathlib import Path
from itertools import product, combinations
import json
import sympy as s

I = s.I
eye = s.eye(2)
X = s.Matrix([[0, 1], [1, 0]])
Z = s.diag(1, -1)
P = s.diag(1, I)
H = (1 + I) * s.Matrix([[1, 1], [1, -1]]) / 2

def key(m):
    return tuple(s.expand(x) for x in m)

def closure(gens):
    found = {key(eye): eye}
    pending = [eye]
    while pending:
        x = pending.pop()
        for g in gens:
            y = x * g
            k = key(y)
            if k not in found:
                found[k] = s.Matrix(2, 2, k)
                pending.append(found[k])
    return found

W = closure([I * eye, X, Z])
R = closure([P, H])
assert len(W) == 16 and len(R) == 96
assert set(W) <= set(R)
assert P**2 == Z and s.simplify(H**2) == I*eye
assert s.simplify((H*P)**3) == -eye
reflections = [x for x in W.values() if x != eye and s.expand((x-eye).det()) == 0]
assert len(reflections) == 6
assert all(s.simplify(x*x) == eye for x in reflections)
cent_sizes = []
for x in reflections:
    c = [r for r in R.values() if key(r*x) == key(x*r)]
    assert all(key(a*b) == key(b*a) for a in c for b in c)
    cent_sizes.append(len(c))
assert cent_sizes == [16]*6

lines = [s.Matrix(v) for v in [(1,0),(0,1),(1,1),(1,-1),(1,I),(1,-I)]]
for v in lines:
    fix = [x for x in W.values() if key(x*v) == key(v)]
    assert len(fix) == 2

def norm(z):
    return int(s.expand(z*s.conjugate(z)))

line_norms = sorted({norm(s.Matrix.hstack(v,w).det()) for v,w in combinations(lines,2)})
assert line_norms == [1,2,4]
other_norms = sorted({norm((x-eye).det()) for x in W.values()
                      if x != eye and s.expand((x-eye).det()) != 0})
assert all(n > 0 and (n & (n-1)) == 0 for n in other_norms)

# The eight linear characters of W/[W,W], and the actual R-action.
forms = {}
for a,b,c in product(range(4),range(2),range(2)):
    forms[key(I**a * X**b * Z**c)] = (a,b,c)
chars = list(product([-1,1],repeat=3))
actions = []
for g in [P,H]:
    gi = g.inv()
    action = []
    for chi in chars:
        values = []
        for x in [I*eye,X,Z]:
            a,b,c = forms[key(gi*x*g)]
            values.append(chi[0]**a * chi[1]**b * chi[2]**c)
        action.append(chars.index(tuple(values)))
    actions.append(action)
unseen = set(range(8)); orbits=[]
while unseen:
    seed=min(unseen); orb={seed}; todo=[seed]
    while todo:
        x=todo.pop()
        for a in actions:
            y=a[x]
            if y not in orb: orb.add(y);todo.append(y)
    unseen-=orb;orbits.append(len(orb))
assert sorted(orbits)==[1,1,3,3]

# Realifications preserve the actual D4 root set.
roots4 = set()
for i,j in combinations(range(4),2):
    for a,b in product([-1,1],repeat=2):
        v=[0]*4;v[i]=a;v[j]=b;roots4.add(tuple(v))
def realify(m):
    out=s.zeros(4)
    for i,j in product(range(2),repeat=2):
        z=m[i,j]
        out[2*i,2*j]=s.re(z);out[2*i,2*j+1]=-s.im(z)
        out[2*i+1,2*j]=s.im(z);out[2*i+1,2*j+1]=s.re(z)
    return out
for g in [P,H]:
    rr=realify(g)
    assert {tuple(rr*s.Matrix(v)) for v in roots4}==roots4

simple = [[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],
          [-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],
          [0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],
          [0,0,0,0,-2,2,0,0],[0,0,0,0,0,-2,2,0]]
mat=s.Matrix(simple);cartan=mat*mat.T/4
zeros=sum(int((s.Matrix([v])*cartan*s.Matrix(v))[0]/2)%2==0
          for v in product([0,1],repeat=8))
assert cartan.det()==1 and zeros==136
result = {
    "sympy_version": s.__version__,
    "quartic": {"W_order":len(W),"R_order":len(R),
                "reflection_count":len(reflections),"reflection_centralizers":cent_sizes,
                "reflection_centralizers_abelian":True,
                "line_determinant_norms":line_norms,
                "other_determinant_norms":other_norms,
                "linear_character_orbits":sorted(orbits),
                "actual_D4_roots_preserved":True},
    "o8_cover": {"cartan_determinant":int(cartan.det()),"quadratic_zeros":zeros,
                 "quadratic_nonzeros":256-zeros,
                 "full_orthogonal_order":2*2**12*15*3*15*63}
}
target=Path(__file__).with_suffix(".json")
target.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))

