from itertools import combinations, product
from collections import deque
import json
from pathlib import Path
roots=sorted(tuple(s if k==i else t if k==j else 0 for k in range(4)) for i,j in combinations(range(4),2) for s,t in product((-1,1),repeat=2))
pos={r:i for i,r in enumerate(roots)}
identity=tuple(range(24))
def perm(M):
    values=[]
    for r in roots:
        v=tuple(sum(M[i][j]*r[j] for j in range(4)) for i in range(4))
        assert all(x%2==0 for x in v)
        values.append(pos[tuple(x//2 for x in v)])
    return tuple(values)
def compose(a,b):return tuple(a[b[i]] for i in range(24))
def closure(gens):
    seen={identity}; q=deque([identity])
    while q:
        x=q.popleft()
        for g in gens:
            y=compose(x,g)
            if y not in seen:seen.add(y);q.append(y)
    return seen
def diag(signs):return [[2*signs[i] if i==j else 0 for j in range(4)] for i in range(4)]
def swap(i,j):
    M=diag([1]*4);M[i],M[j]=M[j],M[i];return M
s=[perm(swap(i,i+1)) for i in range(3)]
sign=perm(diag([-1,1,1,1]))
wd=closure(s+[perm(diag([-1,-1,1,1]))])
rt=closure(s+[sign,perm([[2*(i==j)-1 for j in range(4)]for i in range(4)])])
assert len(wd)==192 and len(rt)==1152
v3=perm([[0,0,2,0],[2,0,0,0],[0,2,0,0],[0,0,0,2]])
v4=perm([[0,-2,0,0],[2,0,0,0],[0,0,0,-2],[0,0,2,0]])
v6=compose(perm(diag([-1]*4)),v3)
base=[(1,1,0,0),(1,-1,0,0),(0,0,1,1),(0,0,1,-1)]
def mat2(p):
    b=[roots[p[pos[r]]] for r in base]
    cols=[tuple(b[0][i]+b[1][i] for i in range(4)),tuple(b[0][i]-b[1][i] for i in range(4)),tuple(b[2][i]+b[3][i] for i in range(4)),tuple(b[2][i]-b[3][i] for i in range(4))]
    return tuple(tuple(cols[j][i] for j in range(4))for i in range(4))
def cyclo3restriction(p):
    M=mat2(p)
    a=tuple(M[i][0]-M[i][1] for i in range(4));b=tuple(M[i][1]-M[i][2] for i in range(4))
    assert a[3]==b[3]==0 and sum(a[:3])==sum(b[:3])==0
    return (a[0],b[0],-a[2],-b[2])
report={"D4WeylOrder":len(wd),"D4RootAutomorphismOrder":len(rt)}
for e,v in ((3,v3),(4,v4),(6,v6)):
    c=[g for g in rt if compose(g,v)==compose(v,g)]
    cw=[g for g in wd if compose(g,v)==compose(v,g)]
    d={"root_centralizer_order":len(c),"weyl_centralizer_order":len(cw)}
    if e in (3,6):
        images={cyclo3restriction(g)for g in c}; d["cyclotomic_image_order"]=len(images)
        scalar=closure([v3,perm(diag([-1]*4))])
        assert images=={cyclo3restriction(g)for g in scalar}
        d["image_is_C6"]=True
    else:
        assert len(c)==96 and len(cw)==16
    report[str(e)]=d
path=Path(__file__).with_suffix(".json");path.write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
