"""Bounded exact product-group relaxations of the global Weyl sign.

Research only: retained rows are necessary-condition reservations, never
realizations. A resource cap returns unknown and cannot reject anything.
"""
from functools import lru_cache
from math import prod

CAP = 50000


def mul(a, b, primes):
    return tuple(x*y % p for x, y, p in zip(a, b, primes))


def cyc(q, primes, step=1):
    identity = (1,)*len(primes)
    g = tuple(pow(q, step, p) for p in primes)
    out, x = [], identity
    while not out or x != identity:
        out.append(x)
        x = mul(x, g, primes)
        if len(out) > CAP:
            return None
    return out


@lru_cache(None)
def possible_targets(support, marked_items, target_prime=None, height=None, minimum=1):
    marked = dict(marked_items)
    primes = tuple(marked)
    if not primes:
        return True
    if any(not 1 <= marked[p] < p for p in primes):
        return False  # The valuation-one Weyl theorem gives 1 <= j_p < p.
    # Each marked prime has exact valuation one and primitive multiplicity j.
    fixed = tuple(marked[p]*prod(q for q in primes if q != p) % p for p in primes)
    targets = {tuple(sign*pow(x, -1, p) % p for p, x in zip(primes, fixed))
               for sign in (1, -1)}
    values = {(1,)*len(primes)}
    for q in support:
        if q in marked:
            continue
        if q == target_prime:
            choices = [tuple(pow(q, h, p) for p in primes) for h in range(minimum, height+1)]
        else:
            choices = cyc(q, primes)
            if choices is None:
                return None
        if len(values)*len(choices) > CAP*10:
            return None
        values = {mul(a, b, primes) for a in values for b in choices}
        if len(values) > CAP:
            return None
    return bool(values & targets)


def global_constraints(row, maximum_height=8):
    """Return an exact failure, or certified lower bounds on unknown heights."""
    support, classes, edges, nines = row
    touching = {q for edge in edges for q in edge}
    marked = tuple((p, c[2]//(p-1)) for p,c in zip((q for q in support if q>=5),classes)
                   if c[0] == 1 and p not in touching)
    if not marked:
        return None, {}
    for p, j in marked:
        if not 1 <= j < p:
            return {'reason': 'valuation_one_primitive_multiplicity',
                    'prime': p, 'primitive_multiplicity': j}, {}
    if possible_targets(tuple(support), marked) is False:
        return {'reason':'global_weyl_sign', 'marked':marked}, {}
    heights = {}
    for q in support:
        if q in dict(marked):
            continue
        # A nine-row forces H3>=2. An isolated absent prime cannot have
        # height one: its active pure denominator would then be that prime.
        minimum = 2 if q == 3 and nines > 0 or q not in touching else 1
        bound = minimum
        for h in range(minimum, maximum_height+1):
            if possible_targets(tuple(support), marked, q, h, minimum) is False:
                bound = h+1
            else:
                break
        if bound > minimum:
            heights[q] = bound
    return None, heights

