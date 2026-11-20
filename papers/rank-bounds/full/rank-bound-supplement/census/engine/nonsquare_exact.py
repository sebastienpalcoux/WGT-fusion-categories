"""Arbitrary-precision replay of one nonsquare search task."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, groupby
from math import gcd, isqrt, prod


from defect import defect_tiles as tile_counts
from denominators import smooth_interval, smooth_over, terminal
from primes import primes_up_to
from starts import start_pairs, top_multiplicities
from prime_cut_prefix import prime_cut_prefix_failure


def possible_packets(prefix, slots, primes, q, fiber):
    runs = [(b, sum(1 for _ in block)) for b, block in groupby(prefix)]
    debt = 0
    for p in primes:
        if p == 2 or q % p == 0:
            continue
        step = (p - 1) // 2
        eligible = [n for b, n in runs if b % p == 0]
        tail = runs[-1][1] if runs[-1][0] % p == 0 else 0
        single = sum(n // step for n in eligible)
        debt += max(0, fiber * (p - 1) - single * step - tail % step)
        if debt > slots:
            return False
        if q % 2 == 0:
            double = sum(n // (2 * step) for n in eligible)
            single += (slots + tail % step) // step
            double += (slots + tail % (2 * step)) // (2 * step)
            if single + double < 6:
                return False
    return True


def defect_tiles(denominators, slots, q):
    return tile_counts(Counter(denominators), q, last=denominators[-1], slots=slots)


def joint_packets(denominators, primes, q, fiber):
    active = [p for p in primes if p != 2 and q % p]
    for size in range(1, len(active)+1):
        for group in combinations(active, size):
            cost = fiber * sum(p-1 for p in group)
            if sum(gcd(z, prod(group)) > 1 for z in denominators) < cost:
                return False
    return True


def remaining_mass_possible(L, E, lo, slots, prefix, primes, q, fiber):
    """Round each mass upwards on an integer grid, as in the native search."""
    # The support is fixed in both users of this helper. The residual
    # E/L^2 includes the vacuum; M>=L and M>=6*lo bound its mass.
    if prime_cut_prefix_failure(prefix, slots, tuple(sorted(primes)), lo,
            Fraction(E, L*L), vacuum_upper=Fraction(1, max(L*L, 36*lo*lo))) is not None:
        return False
    scale = 1 << 32
    numerator = L*L*scale
    def weight(z):
        return (numerator-1)//(z*z)+1
    base = weight(lo)
    unit = (numerator-1)//(36*lo*lo)+1
    target = E*scale
    def affordable(missing, threshold):
        missing = max(0, missing)
        return missing <= slots and (slots-missing)*base + missing*weight(threshold) + unit >= target
    for p in primes:
        part, value = 1, L
        while value % p == 0:
            value //= p
            part *= p
        part = max(p, part)
        needed = 3 if p in (2,5) else 2
        missing = needed - sum(z%part == 0 for z in prefix)
        if not affordable(missing, ((lo-1)//part+1)*part):
            return False
    active = [p for p in primes if p != 2 and q % p]
    for size in range(1, len(active)+1):
        for group in combinations(active,size):
            product = prod(group)
            missing = fiber*sum(p-1 for p in group) - sum(gcd(z,product)>1 for z in prefix)
            threshold = min(((lo-1)//p+1)*p for p in group)
            if not affordable(missing, threshold):
                return False
    return True


def seed_mass_possible(q, b, m, rank, primes, fiber):
    """Exact optimistic upper bound using the unramified-prime union quotas."""
    active = {p: (b//p + 1)*p for p in primes if p!=2 and q%p}
    thresholds = sorted({b+1, *active.values()}, reverse=True)
    assigned = 0
    upper = Fraction(1,36) + m
    for threshold in thresholds:
        candidates = [p for p in active if active[p] >= threshold]
        required = assigned
        for size in range(1, len(candidates)+1):
            for group in combinations(candidates, size):
                demand = fiber * sum(p-1 for p in group)
                if any(b%p == 0 for p in group):
                    demand -= m
                required = max(required, demand)
        remaining = rank-1-m
        if required > remaining:
            return False
        if threshold == b+1:
            required = remaining
        upper += (required-assigned)*Fraction(b*b, threshold*threshold)
        assigned = required
    return q*b*b <= upper


def complete(q, b, m, rank):
    """Enumerate a complete nonsquare task with arbitrary-precision integers."""
    yield from _complete(q, b, m, rank, None, [b])


class InvalidPrefix(ValueError):
    pass


def complete_prefix(q, b, m, rank, primes, prefix):
    """Continue one native node on precisely its fixed support and prefix."""
    if (not primes or list(primes) != sorted(set(primes)) or any(p not in primes_up_to(rank) for p in primes)
            or not 1 <= len(prefix) <= rank-3 or prefix[0] != b
            or any(type(z) is not int or z < 1 for z in prefix) or list(prefix) != sorted(prefix)
            or any(not smooth_over(z, primes) for z in prefix)):
        raise InvalidPrefix('invalid exact-prefix handoff')
    yield from _complete(q, b, m, rank, tuple(primes), list(prefix))


def _complete(q, b, m, rank, fixed_primes, prefix):
    if (any(type(x) is not int for x in (q, b, m, rank)) or rank < 2
            or not 2 <= q < rank or b < 1 or not 0 <= m < rank):
        raise ValueError('expected rank >= 2, 2 <= q < rank, b >= 1 and 0 <= m < rank')
    if (b,q) not in start_pairs(rank):
        return
    multiplicities = top_multiplicities(rank, b, q)
    if m and m not in multiplicities:
        return
    top_limit = max(multiplicities)
    f = q*b*b
    all_primes = primes_up_to(rank)
    fiber = 2 if q == 2 else (q+3)//2
    must = tuple(p for p in all_primes if q*b % p == 0)
    L, E = b, f-1
    for z in prefix[1:]:
        g = gcd(z, L)
        k, e = z//g, L//g
        L, E = L*k, E*k*k-e*e
        if E <= 0:
            raise ValueError('nonpositive prefix residual')
    fixed_run = prefix.count(b)
    if fixed_run > top_limit or m and (fixed_run > m or prefix[-1] > b and fixed_run != m):
        return

    def visit(L, E, slots, primes):
        if not possible_packets(prefix, slots, primes, q, fiber):
            return
        if not defect_tiles(prefix, slots, q):
            return
        first_run = prefix[-1] == b
        if first_run and len(prefix) > top_limit:
            return
        lo = max(prefix[-1], isqrt(L*L//E) + 1)
        hi = isqrt((36*slots+1)*L*L//(36*E))
        if first_run and m:
            if len(prefix) < m:
                hi = min(hi, b)
            else:
                lo = max(lo, b+1)
        if lo > hi:
            return
        if not remaining_mass_possible(L, E, lo, slots, prefix, primes, q, fiber):
            return
        if slots == 2:
            for z, k, d, e in terminal(L, E, lo, hi, primes, fixed_support=True):
                M = k*L
                row = [1, d, e] + [M // t for t in reversed(prefix)]
                if m and row[1:].count(row[-1]) != m:
                    continue
                if any(sum(y % p == 0 for p in primes) < 2 for y in row[1:]):
                    continue
                if any(sum(y % p != 0 for y in row[1:]) < (3 if p in (2,5) else 2)
                       for p in primes):
                    continue
                denominators = prefix + [z, M//d]
                if not defect_tiles(denominators, 0, q):
                    continue
                if not joint_packets(denominators, primes, q, fiber):
                    continue
                if possible_packets(denominators, 0, primes, q, fiber):
                    yield row
            return
        for z in smooth_interval(lo, hi, primes):
            g = gcd(z, L)
            k, e = z//g, L//g
            prefix.append(z)
            yield from visit(k*L, E*k*k-e*e, slots-1, primes)
            prefix.pop()

    def supports(start, chosen, budget):
        if len(chosen) >= 3:
            yield tuple(sorted(chosen))
        for i in range(start, len(optional)):
            p = optional[i]
            cost = 0 if p == 2 else p-1
            if cost > budget:
                break
            yield from supports(i+1, chosen+(p,), budget-cost)

    required = must if fixed_primes is None else fixed_primes
    if not set(must).issubset(required):
        return
    # The jointly inactive fiber has at least four labels for q=3.
    inactive = 4 if q == 3 else fiber
    budget = (rank-inactive)//fiber - sum(p-1 for p in required if p != 2 and q % p)
    if budget < 0:
        return
    optional = [p for p in all_primes if p not in required] if fixed_primes is None else []
    for primes in supports(0, required, budget):
        if m and not seed_mass_possible(q,b,m,rank,primes,fiber):
            continue
        yield from visit(L, E, rank-1-len(prefix), primes)
