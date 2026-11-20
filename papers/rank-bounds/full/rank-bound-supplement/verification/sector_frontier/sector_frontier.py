"""Exact universal starting obstructions, no enumeration of degree lists."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm, prod
import argparse
import json


def primes_to(n):
    return tuple(p for p in range(2,n+1) if all(p % a for a in range(2,isqrt(p)+1)))


def costs(q, primes):
    ram = [(p+1)//2 for p in primes if p>2 and q%p==0]
    C = prod(ram)
    Bq = max([C*(1+(q%2==0))] + [C//h*(h+1) for h in ram])
    return C, Bq


def attainable_mass(rank, b, m, q, Mfloor, reservations):
    """Componentwise maximum of sorted lists is a necessary order-statistic floor."""
    base = [b]*m + [b+1]*(rank-m-1) + [Mfloor]
    for count_floors in reservations:
        other = [floor for count,floor in count_floors for _ in range(count)]
        if len(other)>rank:return False
        other = sorted([b]*(rank-len(other))+other)
        base = list(map(max,base,other))
    if sum(c <= b for c in base)<m:return False
    counts=Counter(base)
    common=lcm(*counts)
    return sum(n*(common//d)**2 for d,n in counts.items())>=q*common**2


@lru_cache(maxsize=2000)
def facet_plans(P,capacity):
    """All antichains covering P whose single summands fit capacity."""
    sets=[]
    for k in range(1,len(P)+1):
        for part in combinations(P,k):
            w=prod(part)
            if w<=capacity:
                sets.append((frozenset(part), w))
    plans=[]
    def visit(i, family, union):
        if union==set(P):
            plans.append(tuple(family))
        for j in range(i,len(sets)):
            a,w=sets[j]
            if all(not(a<=old or old<=a) for old,v in family):
                visit(j+1,family+[(a,w)],union|a)
    visit(0,[],set())
    return tuple(plans)


@lru_cache(maxsize=2000)
def supports(rank,b,q):
    primes=primes_to(rank)
    C,Bq=costs(q,primes)
    required={2}|{p for p in primes if (q*b)%p==0}
    # Every odd extra prime is allowed and contributes p-1 in the
    # packet reservation. All sums above the rank stay above it.
    extra=[p for p in primes if p not in required]
    out=[]
    def visit(i,S):
        P=tuple(p for p in sorted(S) if q%p)
        if Bq*(1+sum(p-1 for p in P if p>2))>rank or C*sum(P)>rank:
            return
        if len(S)>=3:
            out.append((tuple(sorted(S)), P))
        for j in range(i,len(extra)):
            visit(j+1,S|{extra[j]})
    visit(0,required)
    return tuple(out)


def _part(n,p):
    t=1
    while n%(p*t)==0:t*=p
    return t


@lru_cache(maxsize=None)
def _free_min(p):
    # D=0 mod p^2. At least two nonunits omit p (deletion);
    # all others contribute zero modulo p^2.
    q=p*p
    squares={i*i % q for i in range(q) if i%p}
    values={0}
    for n in range(1,46):
        values={(x+y)%q for x in values for y in squares}
        if n>=2 and q-1 in values:return n
    return 46


def structural_floors(S,b,powers,Mfloor):
    rows=[]
    def multiple_at_least(x, floor):
        return x*((floor+x-1)//x)
    for mask in range(1,2**len(S)-1,2):
        left=[p for j,p in enumerate(S) if mask>>j&1]
        right=[p for p in S if p not in left]
        if not right:continue
        products=sorted(x*y for x in left for y in right)
        first=set()
        for k in products:
            c=multiple_at_least(k,b)
            first.update((c,c+k,c+2*k))
        rows.append([(1,Mfloor)]+[(1,k) for k in sorted(first)[:3]])
    # Each t-free witness charge competes with at most |S|-2 charges.
    charged=[]
    for t in S:
        floor=powers[t]
        largest_other=max(p for p in S if p!=t)
        k=_free_min(t)
        items=[multiple_at_least(floor,b)]*(k-1)
        items.append(multiple_at_least(floor*largest_other,b))
        charged+=items
        rows.append([(1,Mfloor)]+[(1,k) for k in items])
    charged.sort(reverse=True)
    capacity=len(S)-2
    rows.append([(1,Mfloor)]+[(1,k) for k in charged[::capacity]])
    return rows


def frontier(rank=45, include_plans=False):
    if rank!=45:raise ValueError('This documented frontier is for rank 45')
    primes=primes_to(rank)
    by_start={}; open_rows=[]; cover=[]
    for b in range(2,isqrt(rank-2)+1):
        for q in range(1,(rank-2)//(b*b)+1):
            f=q*b*b
            if f<8 or f in (8,9) or any(q%(p*p)==0 for p in primes):continue
            C,Bq=costs(q,primes)
            s=Counter()
            for S,P in supports(rank,b,q):
                if len(S)==3 and b in primes and (q==1 or q==2 and b%2):
                    # Actual prime-square rows and equal dyadic blocks.
                    continue
                lower={p:_part(b,p) for p in S}
                for p in S:
                    lower[p]=max(lower[p],p**(1+(len(S)==3)))
                Mfloor=max(6*b,prod(lower.values()))
                geometric=structural_floors(S,b,lower,Mfloor)
                maxmult=min([f-1]+[f-1-(p-1)//2 for p in primes if p>2 and f%p==0 and f%(p*p)])
                pure=next((p for p in S if _part(f,p)==f),None)
                if pure:
                    maxmult=min(maxmult,f-isqrt(pure)-1)
                plans=facet_plans(P,rank//C)
                for m in range(1,maxmult+1):
                    if q%2==0 and b%2 and m%2:continue
                    if b==2 and q in primes and q%2:
                        h=(q-1)//2
                        if m%h or m>7*h:continue
                    if q==2 and b==3 and m>10:continue
                    if q==1 and b==5 and m!=5:continue
                    if pure and f>=pure*pure and m==f-pure:continue
                    baseline=[(Bq*(p-1),b+(-b)%p) for p in P if p>2]
                    if not attainable_mass(rank,b,m,q,Mfloor,(baseline,)):continue
                    s['base_starts']+=1
                    possible=[]
                    after_sector=False
                    for F in plans:
                        s['base_plans']+=1
                        slots=[(C*w,b+(-b)%w) for a,w in F]
                        if attainable_mass(rank,b,m,q,Mfloor,(baseline,slots)):
                            s['sector_plans']+=1
                            after_sector=True
                            if not attainable_mass(rank,b,m,q,Mfloor,(baseline,slots,*geometric)):
                                continue
                            s['with_witnesses_plans']+=1
                            if q==1:
                                isolated=[w for part,w in F if len(part)==1 and w%2]
                                if b in primes and b not in isolated:
                                    # Every prime denominator in square dimension has valuation one.
                                    # A present prime at valuation one is isolated.
                                    continue
                                pure_slots=[]
                                for p in isolated:
                                    pb=p
                                    while pb<b:pb*=p
                                    pure_slots.append((p,pb))
                                other_slots=[(w,b+(-b)%w) for part,w in F if w not in isolated]
                                # The pure slots are disjoint from the unit and from all
                                # the mixed denominators of any particular prime cut.
                                cuts=geometric[:2**(len(S)-1)-1]
                                augmented=[pure_slots+row for row in cuts]
                                if not attainable_mass(rank,b,m,q,Mfloor,
                                    (baseline,slots,*geometric,pure_slots+other_slots,*augmented)):
                                    continue
                            s['final_plans']+=1
                            possible.append(F)
                    if after_sector:
                        s['sector_starts']+=1
                    if possible:
                        s['final_starts']+=1
                        open_rows.append([q,b,m,list(S),len(possible)])
                        if include_plans:
                            for F in possible:
                                cover.append((q,b,m,S,F))
            if s:by_start[f'{q},{b}']=dict(s)
    total=Counter()
    for row in by_start.values():total.update(row)
    if include_plans:return cover
    return {'rank':rank,'counts':dict(total),'by_q_b':by_start,'surviving_starts':open_rows}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output')
    args=parser.parse_args();result=frontier()
    with open(args.output,'w') as f:json.dump(result,f,indent=2);f.write('\n')
    print(result['counts'],result['by_q_b'])
