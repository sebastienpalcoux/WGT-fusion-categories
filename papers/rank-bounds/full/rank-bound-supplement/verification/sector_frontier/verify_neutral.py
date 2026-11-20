"""Independent exact finite checks of the simultaneous sector reservation."""
from functools import lru_cache
from itertools import combinations,product
from math import prod
from neutral_cap import cap, quotas, parity_placements
from sector_frontier import costs,facet_plans

def check():
    orbit_tests=0
    for sizes in [(1,), (2,), (3,5), (2,3,5), (6,8), (4,8,9), ()]:
        sums={sum(m*h for m,h in zip(coeff,sizes)) for coeff in
                product(*(range(45//h+1) for h in sizes))}
        for n in range(45):
            expected=max(v for v in sums if v<=n)
            assert cap(n,sizes)==expected
            orbit_tests+=1
    sector_tests=0
    for S,q in [((2,3,5,7),1),((2,3,5),3),((2,3,7,13),2),((2,3,5,11),6)]:
        C,Bq=costs(q,S)
        for F in facet_plans(tuple(p for p in S if q%p),45//C):
            for J,count in quotas(S,F,q,C,Bq):
                total=0
                for face,k in F:
                    for n in range(len(face)+1):
                        for K in combinations(face,n):
                            if set(J).intersection(K):total+=C*prod(t-1 for t in K)
                expected=max(Bq*sum(p-1 for p in J),total)
                assert expected==count
                sector_tests+=1
    parity_tests=0
    for keys in [(2,3,4),(4,5,12),(6,9,10),(7,11,13)]:
        hist=dict.fromkeys(keys,1)
        for lower in [1,2,4,8]:
            lower=max(lower,*(b&-b for b in keys))
            expected=set()
            for H in (lower,2*lower):
                for odd,even in product((None,*keys),repeat=2):
                    if (odd is None or odd&-odd==H) and (even is None or even&-even<H):
                        expected.add((odd,even))
            assert expected==set(parity_placements(hist,lower))
            parity_tests+=1
    return {'orbit_capacity_checks':orbit_tests,'exact_sector_counts':sector_tests,
            'parity_placement_checks':parity_tests}
if __name__=='__main__':print(check())
