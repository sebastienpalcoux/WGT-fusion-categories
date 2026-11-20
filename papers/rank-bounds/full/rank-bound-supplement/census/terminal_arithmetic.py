"""The five necessary conditions used after construction of a complete degree list."""
from collections import Counter
from math import gcd

from audit import factors_small

STAGES = ("primary_mass", "prime_square_gcd", "defect_one", "principal_block", "ternary_trace")


def tiles(counts, p):
    """Tile p-free denominators in increasing order, using the fewest mixed blocks."""
    owed = {}
    for b, count in sorted(counts.items()):
        if b % p == 0:
            continue
        left = count - owed.pop(b, 0)
        if left < 0 or (p == 2 and left % 2):
            return False
        if p != 2:
            mixed = (left if p == 3 else 3 * left) % p
            if p == 5 and 2 * mixed > left:
                return False
            if mixed:
                owed[2*b] = 2*mixed
    return not owed


def principal_degrees(classes, codegrees, p):
    """Necessary conditions on the degree of a principal p-block row."""
    free = [(x, count) for x, count in classes.items() if x % p]
    possible = set()
    for y, f in codegrees.items():
        if f % p:
            continue
        for x, count in free:
            g = gcd(x, y)
            if g == 1 or 2*count*g <= p or count*(x//g)**2 >= f:
                break
        else:
            possible.add(y)
    return possible


def valuation(n, p):
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def ternary_trace_failure(classes, M, possible):
    """The four proper five-free entries of a denominator-five row are signs."""
    if M % 15 or M % 25 == 0 or M//5 not in classes:
        return False
    proper = [d for d, count in classes.items() if d % 5 and d != M//5 for _ in range(count)]
    if len(proper) != 4 or any(d % 5 not in (1, 4) for d in proper):
        return False
    modulus = 3**valuation(M, 3)
    sums = {0}
    for d in proper:
        if d in possible:
            signed = d if d % 5 == 1 else -d
            sums |= {(value + signed) % modulus for value in sums}
    return modulus-1 not in sums


def prime_square_gcd(classes, q, M, primes):
    """The absence of partner degrees forces an equal defect-one block at q."""
    if q not in (1, 2, 3, 5):
        return False
    for p in primes:
        n = M//p
        if p == q or n not in classes:
            continue
        if q in (3, 5) and (2*n in classes or (n % 2 == 0 and n//2 in classes)):
            continue
        if gcd(*(x for x in classes if x % p)) != 1:
            return True
    return False


def primary_mass_failure(classes, codegrees, dimension, primes):
    least = min(classes)
    for p in primes:
        eligible_mass = sum(c*x*x for x, c in classes.items() if codegrees[x] % p == 0)
        if 2*p*least*eligible_mass < (p-1)*dimension:
            return True
    return False


def obstruction(row, task):
    classes = Counter(row[1:])
    q = task["q"]
    M = row[-1]*task["B"]
    dim = q*M*M
    primes = tuple(p for p, _ in factors_small(M))
    codegrees = {y: dim//(y*y) for y in classes}
    if primary_mass_failure(classes, codegrees, dim, primes):
        return STAGES[0]
    if prime_square_gcd(classes, q, M, primes):
        return STAGES[1]
    denominators = Counter(M//d for d in row)
    if any(q % p == 0 and not tiles(denominators, p) for p in (2, 3, 5)):
        return STAGES[2]
    ternary = set()
    for p in primes:
        possible = principal_degrees(classes, codegrees, p)
        if gcd(*possible) != 1:
            return STAGES[3]
        if p == 3:
            ternary = possible
    if q == 1 and ternary_trace_failure(classes, M, ternary):
        return STAGES[4]
    return None
