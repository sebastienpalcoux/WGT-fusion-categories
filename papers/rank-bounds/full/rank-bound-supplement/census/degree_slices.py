"""Degree slices of principal blocks, using integer sums of complete Galois orbits.

The input is a complete sorted degree list with unit first and its known
distinct dimension primes. No factorization or floating point is used.
A cap omits only the final cyclic subset-sum test at a large prime power.
"""
from collections import Counter
from math import gcd


def _prime_part(n, p):
    value = 1
    while n % p == 0:
        n //= p
        value *= p
    return value


def _cardinalities(classes, dimension, y, prime, pure_degrees):
    """Relax the exact integral sphere for a possible principal y-slice."""
    codegree = dimension // (y*y)
    multiplicity = dict(classes)[y]
    if codegree % (prime*prime):
        return [0]
    pure = y in pure_degrees
    forced = pure or (bool(pure_degrees) and y % prime != 0)
    allowed = [] if forced else [0]
    terms = []
    for x, count in classes:
        g = gcd(x, y)
        forced_zero = g == 1 or (pure and dimension // (x*x) % (prime*prime) != 0)
        if forced_zero:
            if x % prime:
                return allowed
        else:
            terms.append((x, count, g, x//g))

    candidates = [multiplicity] if forced else range(1, multiplicity+1)
    for a in candidates:
        budget = a*(codegree-a)
        baseline = anchor = 0
        shifts = []
        for x, count, g, u in terms:
            step = prime if u % prime else 1
            residue = (a*g) % step
            z = residue if 2*residue <= step else residue-step
            if abs(z) >= a*g:
                break
            cheapest = u*u*z*z
            baseline += count*cheapest
            if baseline > budget:
                break
            anchor += count*x*u*z
            extra_costs = [u*u*((z+delta)**2-z*z) for delta in (-step,step)
                           if abs(z+delta) < a*g]
            if extra_costs:
                shifts.append((min(extra_costs), x*u*step))
        else:
            slack = budget-baseline
            modulus = gcd(*(size for extra, size in shifts if extra <= slack))
            if ((a+anchor) % modulus == 0 if modulus else a+anchor == 0):
                allowed.append(a)
    return allowed


def slice_failure(degrees, dimension, primes, *, cyclic_limit=200_000):
    """Return a cheap, independently replayable obstruction or None."""
    if (len(degrees) < 2 or degrees[0] != 1 or min(degrees[1:]) <= 1
            or sum(x*x for x in degrees) != dimension):
        raise ValueError("A complete perfect nontrivial degree profile is required")
    classes = tuple(Counter(degrees[1:]).items())
    codegrees = {y:dimension//(y*y) for y, _ in classes}
    if any(dimension % (y*y) for y, _ in classes):
        raise ValueError("Nonintegral formal codegree")
    for p in primes:
        if p < 2 or dimension % p:
            raise ValueError("Prime support inconsistent with dimension")
        power = _prime_part(dimension, p)
        pure = {y for y, f in codegrees.items() if _prime_part(f,p) == f}
        slices = []
        for y, _ in classes:
            sizes = _cardinalities(classes, dimension, y, p, pure)
            if not sizes:
                return (p, 'forced slice empty', y)
            if any(sizes):
                slices.append((y,sizes))
        if gcd(*(y for y, _ in slices)) != 1:
            return (p, 'principal degree gcd', None)
        if power > cyclic_limit:
            continue
        mask = (1 << power)-1
        reachable = 1
        for y, sizes in slices:
            following = 0
            for a in sizes:
                shift = a*y*y % power
                following |= ((reachable << shift) & mask) | (reachable >> (power-shift))
            reachable = following
        if not (reachable >> (power-1)) & 1:
            return (p, 'principal mass', None)
    return None
