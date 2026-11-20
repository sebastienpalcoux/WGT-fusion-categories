"""Independent finite checks of antichain generation, floor combination and prefix tiling."""
from collections import Counter
from fractions import Fraction
from itertools import chain, combinations, product
from math import gcd,isqrt,prod
import json
from pathlib import Path
import sector_frontier as sf
import sector_prefix as sp

def check_facet_plans(P,capacity):
    actual={frozenset(s for s,w in f) for f in sf.facet_plans(P,capacity)}
    subsets=[frozenset(t) for n in range(1,len(P)+1) for t in combinations(P,n) if prod(t)<=capacity]
    # Enumerate downward-closed simplicial complexes by their maximal faces
    # in a separate bit-mask search over all candidate maximal faces.
    found=set()
    for mask in range(1<<len(subsets)):
        selected=[s for j,s in enumerate(subsets) if mask>>j&1]
        if set().union(*selected)!=set(P):continue
        if any(s<t for s in selected for t in selected):continue
        found.add(frozenset(selected))
    if actual!=found:raise AssertionError((P,capacity))
    return len(actual)

def independent_mass(rank,b,m,q,v,reservations):
    quotas=[Counter({b:m,b+1:rank-1-m,v:1})]
    for groups in reservations:
        counts=Counter()
        for n,y in groups:counts[y]+=n
        length=sum(counts.values())
        if length>rank:return False
        counts[b]+=rank-length
        quotas.append(counts)
    levels=sorted(set(chain.from_iterable(quotas)),reverse=True)
    prev=0;mass=Fraction(0)
    for d in levels:
        necessary=max(sum(n for floor,n in row.items() if floor>=d) for row in quotas)
        if d>b and necessary>rank-m:return False
        mass+=Fraction(necessary-prev,d*d)
        prev=necessary
    return mass>=q

def brute_tile(t,hist,low,slots,forbidden_equal=frozenset()):
    shapes=sp._block_shapes(t)
    target={b:n for b,n in hist.items() if b%t}
    keys=sorted(target)
    used=Counter()
    def match(i,spare):
        if i==len(keys):return True
        a=keys[i]
        if used[a]>target[a]:return False
        if used[a]==target[a]:return match(i+1,spare)
        for take,k,n in shapes:
            if k==1 and a in forbidden_equal:continue
            c=k*a
            if n and c<low and c not in target:continue
            cost=n if c>=low else 0
            if cost>spare:continue
            used[a]+=take
            if n:used[c]+=n
            if match(i,spare-cost):return True
            if n:used[c]-=n
            used[a]-=take
        return False
    return match(0,slots)

def verify():
    verified=set()
    for B in range(2,7):
        for q in range(1,44//(B*B)+1):
            for S,P in sf.supports(45,B,q):
                C,Bq=sf.costs(q,sf.primes_to(45))
                key=P,45//C
                if key not in verified:
                    check_facet_plans(*key)
                    verified.add(key)
    old_mass=sf.attainable_mass
    def checked(*args):
        a=old_mass(*args);b=independent_mass(*args)
        if a!=b:raise AssertionError(args)
        return a
    sf.attainable_mass=checked
    result=sf.frontier()
    here=Path(__file__).resolve().parent
    if result!=json.loads((here/'sector-frontier45.json').read_text()):
        raise AssertionError('start cover differs')
    trials=0
    for t in (2,3,5,7):
        for n1,n2 in product(range(5),repeat=2):
            for b,c in ((2,3),(3,6),(4,8),(4,5)):
                for slots in range(5):
                    hist={b:n1,c:n2}
                    for veto in [frozenset(),frozenset({b}),frozenset({c}),frozenset({b,c})]:
                        a=sp.block_prefix(t,hist,c+1,slots,veto)
                        b2=brute_tile(t,hist,c+1,slots,veto)
                        if a!=b2:raise AssertionError((t,hist,slots,veto,a,b2))
                        trials+=1
    return {'antichain_domains':len(verified),'tiling_crosschecks':trials,
            'frontier_counts':result['counts']}
if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
