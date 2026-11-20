"""Exact odd-square support and multiplicity check. Run with Python 3.9+."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil, isqrt, prod

from quadratic import marked_possible, primitive_quota, two_low_possible, other_edge_possible
from isolated_gauss import gauss_possible
from edge_budget import compatible, isolated_minimum


RANK = 201
SENTINEL = 501
ABSENT = (0, 0, 0)
PRIMES = tuple(p for p in range(3, (RANK + 1) // 2 + 1, 2)
               if all(p % d for d in range(3, isqrt(p) + 1, 2)))


def supported(n, support):
    for p in support:
        while n % p == 0:
            n //= p
    return n == 1


def projection(p, m, a):
    return (m - a + Q((p+1)*a, (p-1)**2)) / p**2


def orbit_partition(p, m, j):
    """Possible inertia orbit lengths among the trivial-twist rows at b=p."""
    size = m - (p-1)*j
    step = (p-1) & -(p-1)
    reachable = {0}
    for h in range(step, size+1, step):
        if (p-1) % h:
            continue
        x = h * pow(j,-1,p) % p
        if h*h + m*min(x,p-x)**2 > h*p*p:
            continue
        for n in range(size-h+1):
            if n in reachable:
                reachable.add(n+h)
    return size in reachable


def fixed_inertia_possible(p, m, j):
    """Squared Weyl budget on proper labels fixed by p-inertia."""
    orbit = (p-1) & -(p-1)
    neutral = m-(p-1)*j
    # Work in orbit multiples to avoid an intermediate rational.
    left = orbit*(2*j-2)-neutral*(orbit-1)
    if left < 0:
        return False
    fixed = max(0,3 - left//(orbit*(orbit-1)))
    if not fixed:
        return True
    limit = isqrt((p*p-1)//(2*m))
    if not limit:
        return False
    smallest = min(min(j*t%p,(-j*t)%p) for t in range(1,limit+1))
    return 2*(fixed*smallest**2+3-fixed) <= j*(p-j)+p*j-m


def profiles(p, support):
    """(valuation class, total multiplicity, active multiplicity) at b=p."""
    yield ABSENT
    if p == 5 and RANK < 225:
        return  # odd:four-five-edge: centered quadratic projection and edge mass.
    step = (p-1) & -(p-1)
    cap = min(RANK, p*p-isqrt(p-1)-1)
    if len(support) >= 4 and p >= 11:
        n0 = prod(t for t in support if t != p)
        for j in range(2, p, 2):
            a = (p-1)*j
            for m in range(a, cap+1):
                if m % step or m == p*p-p:
                    continue
                if (p*abs(p*j-m)*n0 <= j*(p-j)*n0+(p+1)*j
                        and orbit_partition(p,m,j) and fixed_inertia_possible(p,m,j)):
                    yield 1,m,a
    if min(2+2*p*p,p-1+ceil(Q(p*p*(p+1),p-1))) <= RANK:
        for m in range(step, cap+1, step):
            if m != p*p-p:
                for a in range(0, m+1, p-1):
                    if m + max(4*(p-1)**2*(m-a), primitive_quota(p,a)) <= RANK:
                        yield 2,m,a


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
    edges=[(p,t) for p,t in combinations(P,2) if p not in low and t not in low
           and 2*p*t-1<=RANK]
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
                if len(P)==4 and nines>=72:
                    dimensions=sum(p*t for p,t in G)+sum(p for p in solitary if p!=3)
                    if nines+9+2*dimensions>RANK:
                        continue
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
                good=True
                for p,(v,m,a) in positive.items():
                    if v!=2:
                        continue
                    need=projection(p,m,a)-vacuum
                    # Reserve a lower number of hits with 25|b or 49|b,
                    # disjoint from each explicitly non-p-active reservation.
                    e=2*ceil(need*p**4/2)
                    wrong=[]
                    for x,y in G:
                        n=2*(x*y-1) if p not in (x,y) else 2*p*((y if x==p else x)-1)
                        wrong.append((n,first(lambda b:b%(x*y)==0)))
                    for t in solitary:
                        if t!=p:
                            cap=positive.get(t,ABSENT)[2] if t!=3 else nines
                            wrong.append((max(0,isolated_minimum(t)-cap),pure[t]))
                    nw,vw=quota(wrong)
                    if not fits(nw+e,vw+Q(e,first(lambda b:b%(p*p)==0)**2)):
                        good = False
                        break
                    if hits and (a//(p-1))%2:
                        # Compare the necessary opposite-square-class mass with
                        # all slots after removing hits and the other prime's activity.
                        outside=[(n,b) for t in solitary if t not in (p,3)
                                 for n,b in [(max(0,isolated_minimum(t)-positive.get(t,ABSENT)[2]),pure[t])]]
                        remaining=K-nh-quota(outside)[0]
                        top=max(0,ceil((need-vh)*p**4))
                        if remaining<top or Q(top,p**4)+Q(remaining-top,floor**2)<Q(p-1,p*p):
                            good = False
                            break
                if not good:
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
                if not marked_possible(RANK,P,positive,G,nines):
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


def check_case(case):
    from joint_gauss import possible
    support, classes, edges, nines=case
    marked = {p:v for p,v in zip((p for p in support if p>=5),classes) if v[0]}
    if edges:
        if edges != ((3,5),):
            if other_edge_possible(RANK,support,marked,edges,nines,witnesses):
                raise ArithmeticError(("unresolved twist graph",support,edges,nines))
            return 1
        if any(v[0]==2 for v in marked.values()):
            if RANK<225 and marked.get(5)==(2,4,4):
                return 1  # odd:four-five-edge combines phase counts and orbit mass.
            raise ArithmeticError(("unresolved high-valuation class",support,nines))
        profiles = set()
        compatible(RANK,support,marked,edges,nines,witnesses,profiles)
        if not profiles:
            raise ArithmeticError(("missing phase shapes",support,edges,nines))
        groups={}
        for shape in sorted(profiles):
            if len(marked)==2:
                if len(support) not in (4,5) or not nines or not all(v[0]==1 for v in marked.values()):
                    raise ArithmeticError(("unresolved two prime classes",support,nines))
                if len(support)==4 and not two_low_possible(RANK,nines,marked,shape):
                    continue
            groups.setdefault((shape[0],shape[3],shape[2]),[]).append(shape)
        for family in groups.values():
            if len(family)>1:
                envelope=tuple(min(s[i] for s in family) for i in range(6))
                try:
                    survives=possible(RANK,support,nines,envelope,witnesses,marked,profile_shapes=family)
                except (ArithmeticError,ValueError):
                    survives=True  # An unresolved relaxed hypothesis: retain every exact profile.
                if not survives:
                    continue
            for shape in family:
                if possible(RANK,support,nines,shape,witnesses,marked):
                    raise ArithmeticError(("surviving Gauss profile",support,nines,shape))
        return 1
    else:
        high={p for p,v in marked.items() if v[0]==2}
        low=set(marked)-high
        proper=10 if len(support)==4 and len(low)>=2 else (6 if low else 0)
        neutral=()
        if len(support)==3 and nines:
            # Nine-row q-free witnesses are neutral when3 is isolated.
            # The two omission sets are disjoint in three-prime support.
            neutral=tuple((witnesses(q,nines),(9 if nines>=9 else 3)*q*q)
                          for q in support if q!=3)
        if not proper:
            # odd:neutral-orbits; count already reserved witnesses only once.
            neutral+=((max(0,6-sum(n for n,b in neutral)),15),)
        if gauss_possible(RANK,support,{p:v[1] for p,v in marked.items()},
                          {p:v[2] for p,v in marked.items()},nines,proper,high,neutral):
            raise ArithmeticError(("surviving isolated profile",support,classes,nines))
        return 0

def verify(workers=1):
    count,cases=run()
    if workers==1:
        mixed=sum(map(check_case,cases))
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as pool:
            mixed=sum(pool.map(check_case,cases))
    isolated=len(cases)-mixed
    print(f"{count:,} prime profiles through rank {RANK}; {isolated} isolated "
          f"and {mixed} mixed cases reach the final checks; no survivors.")

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs",type=int,default=1)
    args=parser.parse_args()
    if args.jobs<1:
        parser.error("--jobs must be positive")
    verify(args.jobs)
