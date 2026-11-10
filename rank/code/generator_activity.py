"""Exact activity-mass obstruction in a simple perfect integral modular table."""
from collections import Counter
from functools import lru_cache
from math import isqrt


@lru_cache(maxsize=8192)
def prime_divisors(n):
    out = []
    t = 2
    while t*t <= n:
        if n % t == 0:
            out.append(t)
            while n % t == 0:
                n //= t
        t = 3 if t == 2 else t+2
    if n > 1:
        out.append(n)
    return tuple(out)


def mass_obstruction(degrees, packet=True):
    """Return the first violated prime mass inequality, or None.

    Input must be a valid integral candidate type with square divisibility.
    In odd dimension nonunit labels in a twist fiber have distinct duals.
    """
    counts = Counter(degrees)
    D = sum(d*d*m for d,m in counts.items())
    smallest = min(d for d in counts if d>1)
    for p in prime_divisors(D):
        step = 1
        if packet and p != 2:
            step = ((p-1)//2)*(2 if D%2 else 1)
        capacity = sum((m-m%step)*d*d for d,m in counts.items()
                       if d>1 and D//(d*d)%p==0)
        if 2*p*smallest*capacity < (p-1)*D:
            return p, capacity
    return None


def saturated_twelfth_floor(degrees, p):
    """Rank floor for square-free part 3, ten top codegree-12 labels, p>=5."""
    D = sum(d*d for d in degrees)
    M = isqrt(D//3)
    if 3*M*M != D or M%2 or degrees.count(M//2) != 10 or M%p or p<5:
        raise ValueError('inapplicable parameters')
    h = (p-1)//2
    dens = [M//d for d in degrees if d>1 and (M//d)%p==0]
    if not dens:
        return None  # Cauchy alone excludes the profile
    B = min(dens)
    g = min(d for d in degrees if d>1)
    den = 2*p*g*h
    units = (3*(p-1)*B*B+den-1)//den
    return 13+h*units
