"""Enumerate odd-square prime supports, prime classes, edges and nines."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil, isqrt, prod

from edge_budget import compatible, isolated_minimum


RANK = 209
SENTINEL = 501
ABSENT = (0, 0, 0)
PRIMES = tuple(p for p in range(3, (RANK + 1) // 2 + 1, 2)
               if all(p % d for d in range(3, isqrt(p) + 1, 2)))


def supported(n, support):
    for p in support:
        while n % p == 0:
            n //= p
    return n == 1


def profiles(p, support):
    yield ABSENT
    d = (p-1) & -(p-1)
    if len(support) < 4 or (p-1)//d < 5:
        return
    for j in range(d, p-1, d):
        if (p-1) % j == 0 and p*j <= RANK:
            yield (1, p*j, (p-1)*j)


@lru_cache(None)
def witnesses(p, m):
    if not m:
        return 0
    numer = 0
    for a in range(1,(p+1)//2):
        u = next(k for k in range(3,3*p,3) if (k+a)%p==0)
        l = u - 3*p
        best = min((m-j)*l*l+j*u*u + Q(((m-j)*l+j*u)**2,81-m)
                   for j in range(0,m+1,2))
        numer += 2*ceil(best/162)
    return 2*ceil(Q(numer,2*sum(a*a for a in range(1,(p+1)//2))))


def graph_profiles(support, classes):
    P=support
    positive={p:(v,m,a) for p,(v,m,a) in zip((p for p in P if p>=5),classes) if v}
    low={p for p,(v,m,a) in positive.items() if v==1}
    height={p:max(2 if len(P)==3 else 1, positive.get(p,ABSENT)[0]) for p in P}
    vacuum=Q(1,prod(p**v for p,v in height.items())**2)
    rest=[b for b in range(15,SENTINEL,2) if b not in P and supported(b,P)
          and all(b%(p*p) for p in low)]
    def first(test):
        return next((b for b in rest if test(b)),SENTINEL)
    floor=min(rest+[SENTINEL])
    pure={p:first(lambda b:all(b%t for t in P if p!=t)) for p in P}
    occupied=sum(m for v,m,a in positive.values())
    booked=sum((Q(m,p*p) for p,(v,m,a) in positive.items()),Q(0))
    # The rank and mass estimates permit only these two ternary edges.
    edges = [(3,t) for t in (5,7) if 3 in P and t in P]
    def quota(values):
        return sum(n for n,b in values),sum((Q(n,b*b) for n,b in values),Q(0))
    for size in range(len(edges)+1):
        for G in combinations(edges,size):
            solitary=[p for p in P if all(p not in e for e in G)]
            if sum(p*t for p,t in G)+sum(solitary)>(RANK+1)//2:
                continue
            for nines in (range(0,77,2) if 3 in P else (0,)):
                K=RANK-1-occupied-nines
                E=1-booked-Q(nines,81)
                def fits(n,v):
                    return K>=n and E<=vacuum+v+Q(K-n,floor*floor)
                if E<=0 or not fits(0,0):
                    continue
                other = []
                three = []
                for p,t in G:
                    b=first(lambda b:b%(p*t)==0)
                    other.append((2*p*(t-1) if p==3 else 2*(p*t-1),b))
                    if p==3:
                        three.append((4,b))
                for p in solitary:
                    a=positive.get(p,ABSENT)[2] if p!=3 else nines
                    count=max(0,isolated_minimum(p)-a, 2*p-a if p not in positive and p!=3 else 0)
                    if p in positive and (a//(p-1))%2:
                        count=max(count,ceil(Q((p-1)*pure[p]**2,p*p)))
                    (three if p==3 else other).append((count,pure[p]))
                no,vo=quota(other)
                n3,v3=quota(three)
                if not fits(no+n3,vo+v3):
                    continue
                hits=[]
                if len(P)==3 and nines:
                    factor=9 if nines>=9 else 3
                    hits=[(witnesses(p,nines),first(lambda b:b%(factor*p*p)==0))
                          for p in P if p!=3]
                    nh,vh=quota(hits)
                    if not fits(no+nh,vo+vh):
                        continue
                proper=[]
                for p in low:
                    cap=Q(6,(3*p)**2)
                    if len(P)==4:
                        cap=0
                        for t in P:
                            if t==p:
                                continue
                            e=max(height[t], 2 if ((t==3 and nines) or
                                      (t in solitary and t not in positive)) else 1)
                            cap+=Q(2,(p*t**e)**2)
                    proper.append(cap)
                if any(not fits(no+n3+6,vo+v3+cap) for cap in proper):
                    continue
                if G and not compatible(RANK,P,positive,G,nines,witnesses):
                    continue
                yield G,nines


def run():
    count = 0
    detail = []
    for size in range(3,len(PRIMES)+1):
        if 2*sum(PRIMES[:size])-1>RANK:
            break
        for P in combinations(PRIMES,size):
            if 2*sum(P)-1>RANK:
                continue
            for classes in product(*(tuple(profiles(p,P)) for p in P if p>=5)):
                count+=1
                prime={p:m for p,(v,m,a) in zip((p for p in P if p>=5),classes) if v}
                if sum(prime.values())>=RANK:
                    continue
                # Quick upper bound with all unassigned labels at the smallest
                # supported nonprime (>3). Vacuum is smaller than any such label.
                f=9 if 3 in P else min(P)**2
                if sum((Q(m,p*p) for p,m in prime.items()),Q(0))+Q(RANK-sum(prime.values()),f*f)<1:
                    continue
                for G,nines in graph_profiles(P,classes):
                    detail.append((P,classes,G,nines))
    return count,detail


