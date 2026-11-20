"""Principal-block character sums and the p-adic mass of a block."""
from collections import Counter
from functools import lru_cache
from math import gcd
from screen import factor
from degree_rows import row_possible


@lru_cache(maxsize=65536)
def _class_cost(x, count, y, p):
    """Minimum charge from the x-degree class in a possible principal y-row."""
    g = gcd(x, y)
    if x % p == 0:
        return 0
    if g == 1:
        return -1
    u = x // g
    residue = (count * g) % p
    least = min(residue, p - residue)
    if least >= count * g:
        return -1
    # The degree class is Galois stable; its normalized squared norm
    # is an integer congruent to count*g*g. All its entries are nonzero.
    budget = max(count, (least * least + count - 1) // count)
    s = g % p
    budget = max(budget, (count*s*(p-s)+p-2)//(p-1))
    budget += (count*g*g-budget) % p
    return u*u*budget


def _mass_possible(classes, fixed, possible, modulus):
    residual = -1 - sum(count*x*x for x, count in classes if x in fixed)
    if residual % modulus == 0:
        return True
    terms = [(x*x % modulus, count) for x, count in classes
             if x in possible and x not in fixed and x*x % modulus]
    common = gcd(modulus, *(w for w, _ in terms))
    if residual % common:
        return False
    modulus //= common
    if modulus > 16384:
        return True
    target = (residual//common) % modulus
    mask = (1 << modulus) - 1
    bits = 1
    for residue,count in terms:
        w = residue // common
        # Cyclic bounded subset sum by bundles of 1, 2, 4, ..., remainder.
        step = 1
        while count:
            copies = min(count, step)
            delta = (copies*w) % modulus
            bits |= ((bits << delta) & mask) | (bits >> (modulus-delta))
            if bits & (1 << target): return True
            count -= copies
            step *= 2
    return False


def block_failure(row, dimension):
    classes = list(Counter(row[1:]).items())
    by_y = {y: dimension // (y*y) for y,_ in classes}
    for p, exponent in factor(dimension):
        def pure(f):
            while f % p == 0:
                f //= p
            return f == 1
        powers = {y for y, f in by_y.items() if pure(f)}
        if any(not row_possible(classes, dimension, y, p) for y in powers):
            return True
        fixed = set(powers)
        if powers:
            fixed.update(y for y,_ in classes if y % p)
        possible = set()
        for y, _ in classes:
            f = by_y[y]
            if f % (p*p):
                continue
            budget = f-1
            for x,count in classes:
                cost = _class_cost(x,count,y,p)
                if cost < 0 or cost > budget:
                    break
                budget -= cost
            else:
                possible.add(y)
        if not fixed <= possible or gcd(*possible) != 1:
            return True
        if not _mass_possible(classes,fixed,possible,p**exponent):
            return True
    return False
