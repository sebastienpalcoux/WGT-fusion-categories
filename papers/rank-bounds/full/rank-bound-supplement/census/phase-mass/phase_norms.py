"""Necessary phase-column norm conditions for a perfect square-dimensional type."""
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, prod


@lru_cache(maxsize=32768)
def has_norm(entries):
    """An entry (degree, count, divisor) restricts each coefficient to divisor*Z.

    A zero divisor denotes a degree excluded from the principal p-block.
    True includes the explicitly capped large-coefficient case.
    """
    if not entries:
        return True
    mass = sum(count * degree**2 for degree, count, _ in entries)
    weights = sorted(degree*divisor for degree, count, divisor in entries for _ in range(count)
                     if 0 < degree*divisor <= isqrt(mass))
    if not weights:
        return False
    common = gcd(*weights)
    if mass % common**2:
        return False
    mass //= common**2
    weights = [weight//common for weight in weights]
    if weights.count(1) >= 4:  # Lagrange's four-square theorem
        return True
    suffix_gcd = [0] * (len(weights) + 1)
    for index in reversed(range(len(weights))):
        suffix_gcd[index] = gcd(weights[index], suffix_gcd[index+1])

    @lru_cache(maxsize=8192)
    def represented(index, residual):
        if residual == 0:
            return True
        if index == len(weights) or residual % suffix_gcd[index]**2:
            return False
        square = weights[index]**2
        largest = isqrt(residual//square)
        if square * largest**2 == residual:
            return True
        if index+1 == len(weights):
            return False
        if largest > 1000:
            return True  # Skip this potentially costly branch conservatively.
        return any(represented(index+1, residual-square*value**2)
                   for value in range(largest, -1, -1))

    return represented(0, mass)


@lru_cache(maxsize=32768)
def phase_columns_possible(groups, p, M, radical, least):
    """Combine opposite order-p classes with a relaxation of all higher orders.

    A higher-order square class at p^j uses p^(j-1)h labels per
    representative. Repeat that representative p^(j-2) times and gather all
    such classes in one higher-level list. Their norm identities add, with
    total capacity p*h times the size of this list. The relaxation allows
    list entries to have independent coefficients in the sum.
    """
    h = (p-1)//2
    degrees = [degree for degree, _, _, _ in groups]
    weights = [degree**2 for degree in degrees]
    slots = [size for _, size, _, _ in groups]
    divisors = [divisor for _, _, _, divisor in groups]
    modulus = M//p * radical

    def norm_ok(counts):
        return has_norm(tuple((d, count, divisor) for d, count, divisor in zip(degrees, counts, divisors)
                              if count))

    for higher in product(*(range(size//p+1 if allowed else 1) for _, size, allowed, _ in groups)):
        if not norm_ok(higher):
            continue
        higher_mass = sum(a*b for a, b in zip(higher, weights))
        limits = [size-p*count for size, count in zip(slots, higher)]
        matches = {}
        for vector in product(*(range(n+1) for n in limits)):
            mass = sum(a*b for a, b in zip(vector, weights))
            active = h*(2*mass+p*higher_mass)
            if ((mass+h*higher_mass) % modulus or 2*p*least*active < (p-1)*M**2
                    or not norm_ok(vector)):
                continue
            if all(2*a <= n for a, n in zip(vector, limits)):
                return True
            others = matches.setdefault(mass, [])
            if any(all(a+b <= n for a, b, n in zip(vector, other, limits)) for other in others):
                return True
            others.append(vector)
    return False
