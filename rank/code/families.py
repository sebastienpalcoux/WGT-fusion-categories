"""Sufficient divisor-closed WGT families used by the finite certificates."""
from functools import lru_cache
from math import prod


@lru_cache(maxsize=65536)
def factors(n):
    found = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            found.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        found.append((n, 1))
    return tuple(found)


def family(profile):
    """Upper caps can be infinite. Every returned condition is sufficient."""
    repeated = [(p,e) for p,e in profile if e >= 2]
    k = len(repeated)
    if k < 3:
        return 'at most two repeated primes'
    primes, e = zip(*repeated)
    if k == 3:
        if e[0] <= 5 or e[1] <= 4 or e[2] <= 3:
            return 'ordered three-prime family'
        if all(x <= y for x,y in zip(sorted(e), (4,5,7))):
            return 'three-prime box (7,5,4)'
    if k == 4:
        d,c,b,a = sorted(e)
        if b <= 3 or (e[0] <= 3 and e[1] <= 3):
            return 'four-prime cubic family'
        if d == 2 and c <= 3:
            if a < 12 and a//2 + b//2 <= 5:
                return 'four-prime height sum'
            if primes[0] == 2 and e[0] >= 4 and a <= 9 and b <= 5:
                return 'four-prime dyadic box'
    if k == 5:
        s = sorted(e)
        if s[-1] <= 3 or (s[-1] <= 4 and s[-2] <= 3 and s[0] == 2):
            return 'five-prime family'
    return None


@lru_cache(maxsize=65536)
def dimension_family(D):
    profile = factors(D)
    reason = family(profile)
    if reason:
        return reason
    squarefull = prod(p**e for p,e in profile if e >= 2)
    if squarefull < 29160000 or (squarefull % 2 and squarefull < 6087156075):
        return 'square-full part below threshold'
    return None

