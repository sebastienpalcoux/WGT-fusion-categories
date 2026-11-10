"""Divisibility of a local Gauss sum, for the ordinary integral modular table."""
from collections import Counter
from math import gcd, lcm, prod
from fiber_packets import factors


def obstructions(degrees):
    """List prime obstructions, conservatively skipping possible higher-order twists."""
    counts = Counter(degrees)
    dimension = sum(d*d*m for d, m in counts.items())
    valuations = list(factors(dimension))
    square_class = prod(p for p, e in valuations if e % 2)
    root = prod(p**(e//2) for p, e in valuations)
    delta = 2 if dimension % 2 else 1
    results = []
    for p, e in valuations:
        if p == 2 or e % 2:
            continue
        h = delta * (p-1)//2
        eligible = [(d, m, dimension//(d*d)) for d, m in counts.items()
                    if d > 1 and m >= h and dimension//(d*d) % p == 0]
        if not eligible or any(m >= p*h and f % (p*p) == 0 for d, m, f in eligible):
            continue
        common_square = gcd(*(d*d for d, m, f in eligible))
        active_capacity = sum((m//h)*h*d*d for d, m, f in eligible)
        modulus = lcm(delta*common_square, square_class*root//p)
        if active_capacity < (p-1)*modulus:
            results.append((p, modulus, active_capacity))
    return results
