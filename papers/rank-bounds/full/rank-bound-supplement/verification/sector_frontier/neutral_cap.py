"""Loss of coordinate capacity in simultaneous odd-primary sectors."""
from functools import lru_cache
from itertools import combinations
from math import prod

@lru_cache(maxsize=None)
def cap(n, parts):
    """Largest integer <= n representable by eligible orbit sizes."""
    reachable=1
    for h in parts:
        for _ in range(n//h):
            reachable |= reachable << h
    return (reachable & ((1<<(n+1))-1)).bit_length()-1

@lru_cache(maxsize=None)
def quotas(S,F,q,C,Bq):
    primes=tuple(p for p in S if p>2 and q%p)
    groups=[]
    for i in range(1,len(primes)+1):
        for J in combinations(primes,i):
            facets=C*sum(w-prod(t for t in face if t not in J) for face,w in F)
            disjoint=Bq*sum(p-1 for p in J)
            groups.append((J,max(disjoint,facets)))
    return tuple(groups)

def parity_placements(hist,H):
    # None selects a suffix label. The exact two-power part of M is at
    # least H. An odd-degree label's denominator has that full part.
    power={b:b & -b for b in hist}
    for odd in (None,*hist):
        if odd is not None and power[odd]<H:continue
        for even in (None,*hist):
            if odd is not None and even is not None and power[even]>=power[odd]:continue
            yield odd,even

def improve(S,F,q,hist,H,ell,Mfloor,baseline,C,Bq,multiple):
    h=len(baseline); easy=[]; strict=[]
    for J,quota in quotas(S,F,q,C,Bq):
        hit={b:tuple((p-1)//2 for p in J if b%p==0) for b in hist}
        used=sum(cap(m,hit[b]) for b,m in hist.items())
        floor=min(multiple(S,ell,p) for p in J)
        need=max(0,quota-used)
        if need>=h:return None
        easy.append((need,floor))
        if all(p%4==1 for p in J):strict.append((J,hit,quota,floor))
    for need,floor in easy:
        row=sorted([Mfloor]+[floor]*need+[ell]*(h-1-need))
        baseline=list(map(max,baseline,row))
    # A full two-power in the denominator is required at an odd-degree
    # fixed suffix label. With the actual dyadic height unknown, use H.
    floor_odd=multiple(S,ell,H)
    options=[]
    for odd,even in parity_placements(hist,H):
        booked=[Mfloor]
        if odd is None:booked.append(floor_odd)
        if even is None:booked.append(ell)
        if len(booked)>h:continue
        positions=baseline.copy()
        for J,hit,quota,floor in strict:
            used=sum(cap(m-int(b==odd)-int(b==even),hit[b]) for b,m in hist.items())
            need=max(0,quota-used)
            if need+len(booked)>h:break
            row=sorted(booked+[floor]*need+[ell]*(h-need-len(booked)))
            positions=list(map(max,positions,row))
        else:
            options.append(positions)
    if not options:return None
    return list(map(min,zip(*options)))
