"""Exact common-support and primary-twist obstructions in the simple integral case.

The input is a complete degree multiset with a unique unit, positive integral
formal codegrees and square-full dimension. Supplied prime factors must be the
certified prime factorization of the dimension (not of its square root).
A return value of None means that these necessary conditions pass.
"""
from collections import Counter
from itertools import combinations_with_replacement
from math import gcd, prod


def prime_factors(n):
    """Trial division, intended for the smooth dimensions in the reference census."""
    if type(n) is not int or n < 1:
        raise ValueError("expected a positive integer")
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append((n, 1))
    return factors


def _eligible_columns(classes, scale, codegree, p):
    """Superset of the nonunit degree classes in the principal p-block."""
    p_units = [a for a in classes if a != 1 and a % p]
    possible = []
    for c in classes:
        if c == 1 or codegree[c] % (p * p):
            continue
        if any(a not in scale[c] for a in p_units):
            continue
        energy = 0
        for a in p_units:
            count = classes[a]
            residue = gcd(a, c) % p
            # The normalized character entries reduce to residue at p.
            lower = max(count, (count * residue * (p - residue) + p - 2) // (p - 1))
            lower += (count * residue**2 - lower) % p
            energy += scale[c][a]**2 * lower
        if energy < codegree[c]:
            possible.append(c)
    return p_units, possible


def _common_hit_failure(classes, allowed, p_units, p, exponent):
    """Look for a tuple violating the p-adic common-hit inequality."""
    height = {}
    for a in classes:
        n = a
        height[a] = 0
        while n % p == 0:
            height[a] += 1
            n //= p
    bearing = [a for a in classes if height[a]]
    # Inserting every p-unit row costs nothing in the valuation inequality.
    pool = tuple(c for c in allowed[1] if all(c in allowed[a] for a in p_units))
    work = 0

    def search(start, candidates, cost, chosen):
        nonlocal work
        work += 1
        if work > 2000:
            return None  # Optional tuple tests are safely skipped.
        for c in candidates:
            if 2 * height[c] <= sum(min(height[c], height[a]) for a in chosen):
                break
        else:
            return {"prime": p, "p_unit_degrees": p_units,
                    "additional_degrees": list(chosen), "possible_columns": list(candidates)}
        for j in range(start, len(bearing)):
            a = bearing[j]
            if cost + height[a] >= exponent:
                continue
            hits = tuple(c for c in candidates if c in allowed[a])
            if len(hits) == len(candidates):
                # Same support and a weaker height inequality.
                continue
            failure = search(j + 1, hits, cost + height[a], chosen + (a,))
            if failure is not None:
                return failure
            if work > 2000:
                break
        return None

    return search(0, pool, 0, ())


def _phase_bounds(classes, allowed, codegree, p):
    """Upper bounds for active labels, distinct nontrivial phases, and their mass."""
    parent = {a: a for a in classes}
    same_phase = set()

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in combinations_with_replacement(classes, 2):
        if a == b and classes[a] == 1:
            continue
        if a in same_phase and root(a) == root(b):
            continue
        divisor = 0
        for c, ua in allowed[a].items():
            ub = allowed[b].get(c)
            if ub is not None:
                divisor = gcd(divisor, ua * ub)
                if divisor == 1:
                    break
        if divisor != 1:
            same_phase.update((a, b))
            parent[root(a)] = root(b)

    inactive = {root(a) for a in classes if a == 1 or codegree[a] % p}
    if p >= 5:
        inactive.update(root(a) for a in classes if a in same_phase or classes[a] == 1)
    labels = phases = mass = 0
    counted_groups = set()
    h = max(1, (p - 1) // 2)
    for a, count in classes.items():
        group = root(a)
        if group in inactive:
            continue
        if a in same_phase or count == 1:
            slots = count
            if group not in counted_groups:
                phases += 1
                counted_groups.add(group)
        else:
            slots = (count // h) * h
            phases += slots
        labels += slots
        mass += slots * a * a
    return labels, phases, mass


def phase_failure(row, dimension, factors=None):
    """Return a common-support or primary-activity obstruction when found."""
    classes = Counter(row)
    if (classes[1] != 1 or len(classes) < 2 or min(classes) < 1
            or dimension != sum(count * a * a for a, count in classes.items())
            or any(dimension % (a * a) for a in classes)):
        raise ValueError("expected a complete, nontrivial, perfect integral profile")
    if factors is None:
        factors = prime_factors(dimension)
    factors = tuple(factors)
    if (len(set(p for p, _ in factors)) != len(factors)
            or any(p < 2 or e < 2 for p, e in factors)
            or prod(p**e for p, e in factors) != dimension):
        raise ValueError("expected the square-full prime factorization of D")
    codegree = {a: dimension // (a * a) for a in classes}
    least = min(a for a in classes if a > 1)
    scale = {
        a: {c: c // gcd(a, c) for c in classes if c > 1 and
            (a == 1 or (gcd(a, c) > 1 and (a * (c // gcd(a, c)))**2 < dimension))}
        for a in classes
    }
    odd_kernel = [p for p, e in factors if p > 2 and e % 2]
    multiplier = prod((p + 1) // 2 for p in odd_kernel)
    dyadic = 2 if any(p == 2 and e % 2 for p, e in factors) else 1
    fiber = max([dyadic * multiplier] +
                [multiplier // ((t + 1) // 2) * ((t + 3) // 2) for t in odd_kernel])
    for p, exponent in factors:
        p_units, possible = _eligible_columns(classes, scale, codegree, p)
        allowed = {a: {c: scale[a][c] for c in possible if c in scale[a]} for a in classes}
        failure = _common_hit_failure(classes, allowed, p_units, p, exponent)
        if failure is not None:
            return dict(test="principal_product", **failure)
        labels, phases, mass = _phase_bounds(classes, allowed, codegree, p)
        phase_min = max(1, (p - 1) // 2)
        label_min = phase_min
        if p > 2 and exponent % 2 == 0:
            phase_min = p - 1
            label_min = fiber * phase_min
        if (labels < label_min or phases < phase_min
                or 2 * p * least * mass < (p - 1) * dimension):
            return {"test": "primary_phase", "prime": p, "possible_columns": possible,
                    "max_active": labels, "min_active": label_min,
                    "max_phases": phases, "min_phases": phase_min, "max_mass": mass}
    return None
