"""Exact terminal search with a bounded coprime-scale planning pass.

This changes enumeration order only.  A terminal solution has
E*k*k = 1+d*d+e*e with d,e dividing L and gcd(k,10*d*e)=1.
If J is the prime support of k and U is the part of L prime to J,
then d,e divide U and k <= sqrt((1+2*U*U)/E).  At least two allowed
primes lie outside J because both degrees have two prime factors.

For at most MAX_SCALE_CANDIDATES possible k, use this complete scale
list even when E is large.  Otherwise fall back to the unchanged
reference enumeration before yielding anything.  The bound controls
the optimization, not mathematical coverage.
"""
from bisect import bisect_left, bisect_right
from itertools import combinations
from math import isqrt, prod

from denominators import degree_divisors, smooth_interval
from denominators import terminal as reference_terminal

MAX_SCALE_CANDIDATES = 128
SQUARE_RESIDUES_64 = sum(1 << r for r in {x*x % 64 for x in range(64)})


def terminal(L, E, low, high, primes, *, fixed_support):
    if low > high or E % 8 not in (1,2,3,5,6):
        return
    if high-low <= 64:
        yield from reference_terminal(L,E,low,high,primes,fixed_support=fixed_support)
        return
    primes = tuple(sorted(primes))
    missing = tuple(p for p in primes if fixed_support and L % p)
    if 2 in missing or 5 in missing:
        return
    # Since k is odd, E == 3 (mod8) forces both d and e to be odd.
    # Their common divisor bound can then use the odd part of L.
    degree_bound = L
    if E % 8 == 3:
        while degree_bound % 2 == 0:
            degree_bound //= 2
    elif L % 2:
        return
    allowed = tuple(p for p in primes if p not in (2,5))
    plans = []
    for count in range(min(len(allowed),len(primes)-2)+1):
        for group in combinations(allowed,count):
            if not all(p in group for p in missing):
                continue
            U = degree_bound
            for p in group:
                while U % p == 0:
                    U //= p
            base = prod(group)
            maximum = isqrt((1+2*U*U)//E)
            if maximum < base:
                continue
            for extra in smooth_interval(1,maximum//base,group):
                plans.append((U,base*extra))
                if len(plans) > MAX_SCALE_CANDIDATES:
                    yield from reference_terminal(L,E,low,high,primes,fixed_support=fixed_support)
                    return
    for U,k in plans:
        residue = E*k*k-1
        if residue < 72:
            continue
        first = max((k*L+high-1)//high,isqrt((residue-1)//2)+1)
        last = min(k*L//low,isqrt(residue-36))
        if first > last:
            continue
        eligible = degree_divisors(U,primes)
        for i in range(bisect_left(eligible,first),bisect_right(eligible,last)):
            e = eligible[i]
            d2 = residue-e*e
            if not (SQUARE_RESIDUES_64 >> (d2 & 63)) & 1:
                continue
            d = isqrt(d2)
            if d*d == d2 and U % d == 0 and sum(d % p == 0 for p in primes) >= 2:
                yield k*(L//e),k,d,e
