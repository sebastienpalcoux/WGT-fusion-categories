"""Necessary principal-block tests on a complete integral modular dimension type."""
from collections import Counter
from math import gcd, prod
from functools import lru_cache
from odd_blocks import mass_feasible
from families import factors


def projections(ds,p,e,screen=True):
    mult=Counter(ds)
    D=sum(x*x for x in ds)
    pe=p**e
    f={x:D//(x*x) for x in mult}
    free=[x for x in mult if x%p]
    def admissible(y):
        if f[y]%(p*p): return False
        if any(gcd(x,y)==1 for x in free if x>1):return False
        return sum(mult[x]*(x//gcd(x,y))**2 for x in free)<=f[y]
    keys=sorted((y for y in mult if y>1 and admissible(y)),reverse=True)
    # Unit always belongs to block; any field congruence over p applies to the entire Galois orbit.
    chosen={}
    def visit(i,residue):
        if i==len(keys):
            if residue%pe: return
            if any(sum(m*(y//gcd(y,x))**2 for y,m in chosen.items())+1>f[x] for x in free):return
            if any(not any(y%s for y in chosen) for s,_ in factors(D) if s!=p):return
            if p%2 and e%2 and not mass_feasible(chosen,p,D%2,True):return
            yield tuple(sorted(chosen.items()))
            return
        y=keys[i]
        for m in range(mult[y]+1):
            if m: chosen[y]=m
            yield from visit(i+1,residue+m*y*y)
        chosen.pop(y,None)
    yield from visit(0,1)


def possible(ds):
    return {p:sum(1 for _ in projections(ds,p,e)) for p,e in factors(sum(x*x for x in ds))}

if __name__=='__main__':
    import json,sys
    rows=json.load(open(sys.argv[1]));rej=Counter();allc=[]
    for i,ds in enumerate(rows):
        c=possible(ds);allc.append(c);rej[tuple(p for p in c if not c[p])]+=1
    print(rej);print(allc[:6]);json.dump(allc,open('r28new/blocks/counts.json','w'))
