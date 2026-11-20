"""Exact scalar checks for the rank-44 odd-edge isolation extension.
No degree-type enumeration is performed.
"""
from fractions import Fraction as F
from math import isqrt
import json

RANK = 44


def upper(groups):
    room = RANK-sum(n for n, _ in groups)
    if room < 0:
        return F(0)
    return F(room, 36)+sum(F(n, d*d) for n, d in groups)


def square_partitions(total, count, smallest=1):
    if not count:
        if not total:
            yield ()
        return
    for a in range(smallest, isqrt(total//count)+1):
        for tail in square_partitions(total-a*a, count-1, a):
            yield (a,)+tail


def surviving_profiles():
    for u in range(7, RANK, 2):
        if any(u % j == 0 for j in range(3, isqrt(u)+1, 2)):
            continue
        for m in range(14):
            odd = 3 if m >= 7 else 0
            M = 30*u
            base = [(1, M), (m, 4), (18, 15)]
            if upper(base+[(u, u*u), (odd, 8)]) >= 1:
                yield u, m, 0
            for j in range(1, u-1):
                if (u-1) % j or F(m, 16)+F(j, u) >= 1:
                    continue
                fixed = [*base, (u*j, u)]
                if upper(fixed+[(3, 2*u)]) < 1:
                    continue
                if j == 1 and m % 2 and 16*(1-F(m, 16)-F(j, u)+F(1, M)) < 1:
                    continue
                free = RANK-sum(n for n, _ in fixed)
                for size in range(3, free+1):
                    if u == 7 and j == 1 and m in (11,12) and size != 6:
                        continue
                    for shape in square_partitions(j*(u-j), size):
                        for mask in (range(1, 2**size-1) if odd else (0,)):
                            floors = []
                            for i, a in enumerate(shape):
                                dyadic = 8 if mask & (1 << i) else 1
                                k = next(k for k in range(dyadic, 8*u+1, dyadic)
                                         if k > 1 and (a*k) % u in (1, u-1))
                                floors.append((1, u*k))
                            extra = max(0, odd-mask.bit_count())
                            if upper([*fixed, *floors, (extra, 8)]) >= 1:
                                yield u, m, j


assert not list(surviving_profiles())

vac=F(1,210**2)
checks={
 'two_odd_edges_15_21':F(12,225)+F(18,441)+F(14,16),
 'edge21_and_dyadic5':F(16,100)+F(18,441)+F(10,16),
 'four_odd_primes_with_edge':F(11,16)+F(1,15)+F(1,7)+F(1,11),
 '21_absent5':F(13,16)+F(22,441)+F(5,625)+F(3,36)+vac,
 '21_present5_at_m11':F(11,16)+F(22,441)+F(1,5)+F(4,400)+F(1,36)+vac,
 '21_present5_at_m12_raised':F(12,16)+F(1,5)+F(21,441)+F(1,42**2)+F(4,80**2)+F(1,420**2),
 'dyadic7_with5':F(8,16)+F(24,196)+F(1,5)+F(4,400)+F(2,36)+vac,
 'dyadic7_without5_m10':F(10,16)+F(24,196)+F(9,36)+vac,
 'twenty_active5_high_four':F(20,625)+F(10,16)+F(3,64)+F(7,49)+F(3,36)+vac,
 'eighteen_active5_m11':F(18,625)+F(7,49)+F(3,36)+F(11,16)+F(3,125**2)+vac,
}
assert all(v<1 for v in checks.values()), checks
assert 1+24+6+11+3==45
out={k:str(v) for k,v in checks.items()}
out['mixed_fifteen_survivors']=0
print(json.dumps(out,indent=2))
