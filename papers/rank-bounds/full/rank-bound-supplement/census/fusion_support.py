"""Necessary coefficient bounds from balancing and cyclic fusion reciprocity.

Input is a validated sorted complete degree list of the manuscript's
nonpointed based-simple modular data. A skipped large product is retained.
The returned obstruction names a product with no possible degree partition.
"""
from collections import Counter
from math import gcd


def least_prime(n):
    if n < 2:
        raise ValueError('expected an integer at least two')
    p = 2
    while p*p <= n:
        if n % p == 0:
            return p
        p = 3 if p == 2 else p + 2
    return n


def bounded_sum(target, capacities):
    """Decide a bounded coin sum by binary subdivision of multiplicities."""
    bits = 1
    mask = (1 << (target + 1)) - 1
    for degree, copies in capacities:
        step = 1
        while copies:
            take = min(copies, step)
            bits = (bits | (bits << (degree * take))) & mask
            copies -= take
            step *= 2
        if bits >> target:
            return True
    return bool(bits >> target)


def fusion_failure(row, *, max_target=600):
    """Return the first failed product, or None if the bounded test passes."""
    if (len(row) < 2 or row[0] != 1 or row[1] <= 1 or
            any(type(x) is not int or x <= 0 for x in row) or
            any(a > b for a, b in zip(row, row[1:]))):
        raise ValueError('expected a sorted perfect degree list')
    counts = Counter(row[1:])
    values = sorted(counts)
    minimum = values[0]
    prime = least_prime(sum(x*x for x in row))

    def directional_bound(a, b):
        if gcd(a, b) == 1:
            return a*b // prime
        return a*b - (1 if a == b else minimum)

    for index, a in enumerate(values):
        for b in values[index:]:
            mass = a*b
            if mass > max_target:
                break
            # Among equal-degree inputs a dual pair always exists. A
            # nondual pair also exists exactly when that degree is repeated.
            vacuum_cases = (0,) if a != b else (1, 0) if counts[a] > 1 else (1,)
            for vacuum in vacuum_cases:
                target = mass - vacuum
                ceiling = directional_bound(a, b)
                if not vacuum:
                    ceiling = min(ceiling, mass - minimum)
                capacities = []
                for c in values:
                    if c > ceiling:
                        break
                    per_label = min(ceiling // c, directional_bound(a, c) // b,
                                    directional_bound(b, c) // a)
                    if per_label:
                        capacities.append((c, counts[c] * per_label))
                if not bounded_sum(target, capacities):
                    return {'a': a, 'b': b, 'vacuum': vacuum, 'target': target,
                            'capacities': [list(pair) for pair in capacities]}
    return None
