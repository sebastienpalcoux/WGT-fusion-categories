"""Prime-support certificates, BDP 4.17, 16.1--16.6.

Input profiles are lists of distinct (prime, upper exponent) pairs. Infinity
is allowed as math.inf for a symbolic exponent. No external dependencies.
A true divisor certificate covers every divisor and any coprime square-free
decoration. A false result is inconclusive. This module does not factor.
"""

from functools import lru_cache
from math import gcd, inf, isqrt


def b(weight):
    return max(5, (weight-1)*weight**2)


def bounded_product(powers, limit):
    """Exact product, or None as soon as the requested positive limit is exceeded."""
    if limit < 1:
        return None
    value = 1
    for p, exponent in powers:
        remaining = limit // value
        if exponent == inf or exponent >= remaining.bit_length():
            return None
        value *= p**exponent
        if value > limit:
            return None
    return value


def heights(profile, reduced):
    return {p: inf if a == inf else (a - int(reduced and p == 2)) // 2
            for p, a in profile}


def crossed_branch(profile, reduced):
    """One witness excluding all tables on this exact repeated support."""
    height = heights(profile, reduced)
    for r, _ in reversed(profile):
        if r == 2:
            continue
        for p, _ in profile:
            if p == r:
                continue
            for t, _ in reversed(profile):
                if t == r or (reduced and t == 2):
                    continue
                used = [(s, h - int(s == t)) for s, h in height.items()
                        if s != r and s != p]
                weight = bounded_product(used, r)
                if weight is not None and (weight == 1 or r > b(weight)):
                    return r, p, t, weight
    return None


def base_cones(profile):
    """The elementary terminal cones used in BDP 16.2."""
    n = len(profile)
    if n <= 2:
        return "Dong"
    if n <= 5 and sum(a > 3 for _, a in profile) <= 5 - n:
        return {3: "BDP 4.12", 4: "BDP 11.3", 5: "BDP 15.5"}[n]
    return None


def extension_certificate(profile, direct=None, symbolic_cap_limit=10000):
    """Try BDP 16.2, 16.4 or 16.6 recursively, returning a proof tree.

    A caller may supply additional proven divisor cones as direct(profile).
    Its truthy return must certify *all divisors* with square-free decoration.
    A profile is sorted and square-free factors are dropped before the call.
    For an infinite cap, symbolic_cap_limit bounds the intermediate finite
    cofactor promotion to try. Increase it to try a larger finite bound.
    """

    @lru_cache(None)
    def visit(items):
        initial = base_cones(items) or (direct(items) if direct else None)
        if initial:
            return {"rule": str(initial), "profile": items}

        prefix, (r, _) = items[:-1], items[-1]
        denominator = prefix[0][0] * prefix[-1][0]
        if len(prefix) == 2:
            denominator *= prefix[-1][0]
        t = bounded_product(heights(prefix, False).items(),
                            isqrt(r * denominator - 1))
        if t is not None:
            proof = visit(prefix)
            if proof:
                return {"rule": "BDP 16.4", "profile": items,
                        "T": t, "subproof": proof}

        ordinary = crossed_branch(items, False)
        if ordinary is not None:
            witnesses = {"ordinary": ordinary}
            has_reduced = dict(items).get(2, 0) >= 3
            if has_reduced:
                witnesses["reduced"] = crossed_branch(items, True)
            if not has_reduced or witnesses["reduced"] is not None:
                children = []
                for i in range(len(items)):
                    proof = visit(items[:i] + items[i+1:])
                    if not proof:
                        break
                    children.append(proof)
                else:
                    return {"rule": "BDP 16.2", "profile": items,
                            "witnesses": witnesses, "proper_supports": children}

        for i, (p, a) in enumerate(items):
            others = items[:i] + items[i + 1:]
            powers = tuple(heights(others, False).items())
            if p > 5:
                # Removing any component bounds U divided by its largest.
                for j in range(len(powers)):
                    w = bounded_product(powers[:j] + powers[j+1:], p)
                    if w is not None and p > b(w):
                        proof = visit(others)
                        if proof:
                            return {"rule": "BDP 16.6", "profile": items,
                                    "prime": p, "prime_bound": b(w),
                                    "subproof": proof}
                        break
            bound = symbolic_cap_limit if a == inf else a
            root = bounded_product(powers, bound)
            if root is None:
                continue
            cap = 2 * (2 if any(s == 3 for s, _ in others) else 1) * root
            cap += int(p == 2) - 1
            if cap >= a or (a == inf and cap > symbolic_cap_limit):
                continue
            subprofile = items[:i] + ((p, cap),) + items[i + 1:]
            proof = visit(subprofile)
            if proof:
                return {"rule": "BDP 16.6", "profile": items, "prime": p,
                        "source_cap": cap, "subproof": proof}
        return None

    return visit(tuple(sorted((p, a) for p, a in profile if a >= 2)))


def missing_pair(degrees, delta, primes, reduced=False):
    """Missing necessary BDP 16.1 witness at an actual complete table type.

    delta is the table dimension, i.e. half the category dimension if reduced.
    Input primes must be all distinct prime divisors of delta. The output is
    None or an excluded triple (r,p,t). Distinct labels are never required.
    """
    degrees = set(degrees) - {1}
    codegrees = {x: delta // (x*x) for x in degrees}
    for r in primes:
        if r == 2:
            continue
        for p in primes:
            if p == r:
                continue
            ys = [y for y in degrees if y % p and codegrees[y] % (r*r) == 0]
            for t in primes:
                if t == r or (reduced and t == 2):
                    continue
                xs = [x for x in degrees if x % r and codegrees[x] % (t*t) == 0]
                compatible = False
                for x in xs:
                    for y in ys:
                        g = gcd(x, y)
                        if g >= 2 and r <= b(g):
                            compatible = True
                            break
                    if compatible:
                        break
                if not compatible:
                    return r, p, t
    return None
