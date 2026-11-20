"""Exact necessary prime-cut conditions; no searches or factorization.

The caller supplies the distinct primes of M. For each nontrivial
bipartition, the nonunit denominators crossing it must have at least
three distinct values and lcm M. The theorem applies in the manuscript's
non-pointed, based-simple, rotation-admissible setting, D=q*M*M.

A supplied support is validated by division; no integer factorizer is used.
The number of cuts is 2**(len(primes)-1)-1. Above max_primes the checker
retains the input, with an explicit diagnostic, rather than excluding it.
"""
from math import lcm


def _data(degrees, M, primes):
    primes = tuple(primes)
    if M < 2 or not primes or tuple(sorted(set(primes))) != primes:
        raise ValueError('M and its distinct sorted prime support are required')
    remainder = M
    for p in primes:
        if p < 2 or M % p:
            raise ValueError('invalid supplied prime support')
        while remainder % p == 0:
            remainder //= p
    if remainder != 1:
        raise ValueError('supplied prime support does not exhaust M')
    denominators = []
    for d in sorted(set(degrees)):
        if d < 1 or M % d:
            raise ValueError('every degree must be a positive divisor of M')
        if d == 1:
            continue
        b = M // d
        mask = sum(1 << j for j, p in enumerate(primes) if b % p == 0)
        denominators.append((b, mask))
    return primes, sorted(denominators)


def prime_cut_diagnostics(degrees, M, primes, *, max_primes=12):
    """Return every failed cut, or an explicit safe-skip diagnostic.

    Primality of supplied support entries is the caller's contract; this
    validates divisibility and exhaustiveness. The census replay selects
    them from a fixed, independently checked prime list.
    """
    primes, denominators = _data(degrees, M, primes)
    if len(primes) > max_primes:
        return {'status': 'skipped', 'reason': 'prime-count resource limit',
                'failures': []}
    whole = (1 << len(primes)) - 1
    failures = []
    tested = 0
    for left in range(1, whole, 2):  # Put the first prime on the left.
        right = whole ^ left
        crossing = [b for b, mask in denominators
                    if mask & left and mask & right]
        tested += 1
        height = lcm(*crossing) if crossing else 1
        if len(crossing) < 3 or height != M:
            failures.append({
                'left_primes': [p for j, p in enumerate(primes) if left >> j & 1],
                'right_primes': [p for j, p in enumerate(primes) if right >> j & 1],
                'crossing_denominators': crossing,
                'distinct_count': len(crossing),
                'crossing_lcm': height,
                'crossing_degree_gcd': M // height,
                'fails_count': len(crossing) < 3,
                'fails_gcd': height != M,
            })
    return {'status': 'checked', 'cuts_tested': tested, 'failures': failures}


def prime_cut_failure(degrees, M, primes, *, max_primes=12):
    """First exact exclusion certificate, or None (including safe skip)."""
    result = prime_cut_diagnostics(degrees, M, primes, max_primes=max_primes)
    return result['failures'][0] if result['failures'] else None


def prefix_count_failure(prefix, remaining, primes):
    """Safe weak prefix test at a fixed final prime support.

    prefix contains only nonunit denominators; remaining counts all future
    nonunit labels. Each future label can add at most one new crossing
    denominator per cut. Deficits from different cuts are NOT added.
    A later rise of M is irrelevant to this count-only condition.
    """
    primes = tuple(primes)
    if remaining < 0:
        raise ValueError('remaining must be nonnegative')
    whole = (1 << len(primes)) - 1
    masks = {b: sum(1 << j for j, p in enumerate(primes) if b % p == 0)
             for b in prefix}
    for left in range(1, whole, 2):
        right = whole ^ left
        count = sum(bool(mask & left and mask & right) for mask in masks.values())
        if count + remaining < 3:
            return {'left_primes': [p for j, p in enumerate(primes) if left >> j & 1],
                    'observed_distinct': count, 'remaining': remaining}
    return None
