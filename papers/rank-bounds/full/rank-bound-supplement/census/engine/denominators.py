"""Shared, arbitrary-precision arithmetic for the ordered denominator searches."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, prod


def smooth_over(value, primes):
    for p in primes:
        while value % p == 0:
            value //= p
    return value == 1


def smooth_interval(low, high, primes):
    if high - low <= 35:
        for value in range(low, high + 1):
            if smooth_over(value, primes):
                yield value
        return
    primes = sorted(primes, reverse=True)

    def visit(i, value):
        if i == len(primes):
            if value >= low:
                yield value
            return
        p = primes[i]
        while value <= high:
            yield from visit(i + 1, value)
            value *= p

    yield from visit(0, 1)


@lru_cache(maxsize=32)
def degree_divisors(value, primes):
    """Divisors with at least two distinct prime factors, in increasing order."""
    divisors = [(1, 0)]
    for p in primes:
        previous = tuple(divisors)
        power = 1
        while value % p == 0:
            value //= p
            power *= p
            divisors.extend((d * power, count + 1) for d, count in previous)
    if value != 1:
        raise ValueError("prefix lcm has a prime outside the declared support")
    return sorted(d for d, count in divisors if count >= 2)


def terminal(L, E, low, high, primes, *, fixed_support):
    """Yield (b,k,d,e) for b=kL/e, d<=e dividing L and Ek²=1+d²+e².

    The smooth scale k is coprime to 10de. If fixed_support, every allowed
    prime must occur in kL. Both d and e have two distinct prime factors.
    """
    if low > high or E % 8 not in (1, 2, 3, 5, 6):
        return
    primes = tuple(sorted(primes))
    missing = tuple(p for p in primes if fixed_support and L % p)
    if 2 in missing or 5 in missing:
        return

    if high-low <= 64:
        for b in range(low, high+1):
            g = gcd(b, L)
            k, e = b//g, L//g
            if gcd(k, 10) != 1 or not smooth_over(k, primes) or any(k % p for p in missing):
                continue
            d2 = E*k*k-e*e-1
            if not 36 <= d2 <= e*e:
                continue
            d = isqrt(d2)
            if d*d != d2 or L % d or gcd(k, d) != 1:
                continue
            if all(sum(x % p == 0 for p in primes) >= 2 for x in (d, e)):
                yield b, k, d, e
        return

    divisors = degree_divisors(L, primes)
    if not divisors:
        return
    if E >= len(divisors):
        # For large E the congruence pairs the divisors sparsely.
        buckets = defaultdict(list)
        for e in divisors:
            buckets[e*e % E].append(e)
        for d in divisors:
            group = buckets.get(-(1+d*d) % E, ())
            for i in range(bisect_left(group,d), len(group)):
                e = group[i]
                k2 = (1+d*d+e*e)//E
                k = isqrt(k2)
                if k*k != k2 or gcd(k,10*d*e) != 1 or not smooth_over(k,primes):
                    continue
                b = k*(L//e)
                if low <= b <= high and all(k % p == 0 for p in missing):
                    yield b, k, d, e
        return

    # For small E, enumerate supp(k). Both degrees must divide its coprime
    # part U in L. Thus E k² <= 1+2U² and at least two primes are outside it.
    allowed = [p for p in primes if p not in (2, 5)]
    for count in range(min(len(allowed), len(primes)-2) + 1):
        for group in combinations(allowed, count):
            if not all(p in group for p in missing):
                continue
            U = L
            for p in group:
                while U % p == 0:
                    U //= p
            base = prod(group)
            max_k = isqrt((1+2*U*U)//E)
            if max_k < base:
                continue
            eligible = degree_divisors(U,primes)
            for extra in smooth_interval(1, max_k//base, group):
                k = base*extra
                residue = E*k*k-1
                if residue < 72:
                    continue
                first = max((k*L+high-1)//high, isqrt((residue-1)//2)+1)
                last = min(k*L//low, isqrt(residue-36))
                for i in range(bisect_left(eligible,first), bisect_right(eligible,last)):
                    e = eligible[i]
                    d2 = residue-e*e
                    d = isqrt(d2)
                    if d*d == d2 and U % d == 0 and sum(d%p == 0 for p in primes) >= 2:
                        yield k*(L//e), k, d, e
