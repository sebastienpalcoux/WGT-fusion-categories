"""Class-sum obstruction for a row in the principal block at a prime."""
from functools import lru_cache
from math import gcd, isqrt


@lru_cache(maxsize=4096)
def residue_cost(p, s, multiplicity):
    """(p-1) times the norm bound for an entry reducing to s at p."""
    if multiplicity >= p:
        return s * (p - s)
    best = p**3
    for e in range(1, multiplicity + 1):
        if (p - 1) % e:
            continue
        h = (p - 1) // e
        lower = s // h
        value = min((s - h*k)**2 + h*k for k in (lower, lower + 1))
        best = min(best, p*value - s*s)
    return best


def row_possible(classes, dimension, y, p, codegree_limit=64):
    """Necessary condition; untested rows above the limit pass.

    classes lists (degree,multiplicity) for the nonunit degree classes.
    The y-row is being tested for membership in the principal p-block.
    """
    f = dimension // (y*y)
    if f > codegree_limit:
        return True
    pure_part = f
    while pure_part % p == 0:
        pure_part //= p
    support_in_block = pure_part == 1
    mult_y = dict(classes)[y]
    terms = []
    for x, mult_x in classes:
        g = gcd(x, y)
        u = x // g
        outside = support_in_block and dimension // (x*x) % (p*p)
        if g == 1 or u*u >= f or outside:
            if x % p:
                return False
            continue
        constrained = u % p != 0
        s = g % p
        bound = min(isqrt(mult_x*(f-1)) // u, mult_x*g - 1)
        choices = []
        for a in range(-bound, bound + 1):
            if constrained and (a - mult_x*g) % p:
                continue
            t, rem = divmod(a, mult_x)
            energy = (mult_x-rem)*t*t + rem*(t+1)**2
            if constrained:
                for shift in (t-1, t, t+1):
                    cost = residue_cost(p, (s-shift) % p, min(mult_x, mult_y))
                    ceil = (mult_x*cost + p-2) // (p-1)
                    energy = max(energy, ceil + 2*shift*a - mult_x*shift**2)
                energy += (mult_x*g*g - energy) % p
            used = energy * u*u
            if used <= f-1:
                choices.append((x*u*a, used))
        if not choices:
            return False
        terms.append(choices)

    # Large coefficients first. Lower norms and ranges ignore correlations,
    # so they can only shorten the exact bounded search.
    terms.sort(key=lambda choices: choices[-1][0]-choices[0][0], reverse=True)
    n = len(terms)
    minimum = [0]*(n+1)
    maximum = [0]*(n+1)
    modulus = [0]*(n+1)
    budget = [0]*(n+1)
    for i in reversed(range(n)):
        choices = terms[i]
        first = choices[0][0]
        minimum[i] = minimum[i+1] + first
        maximum[i] = maximum[i+1] + choices[-1][0]
        modulus[i] = gcd(modulus[i+1], *(v-first for v, _ in choices))
        budget[i] = budget[i+1] + min(c for _, c in choices)

    def feasible(i, target, remaining):
        if budget[i] > remaining or not minimum[i] <= target <= maximum[i]:
            return False
        if modulus[i] and (target-minimum[i]) % modulus[i]:
            return False
        if i == n:
            return True
        return any(feasible(i+1, target-value, remaining-used)
                   for value, used in terms[i])

    return feasible(0, -1, f-1)
