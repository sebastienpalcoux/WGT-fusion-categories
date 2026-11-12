"""Arbitrary-precision replay of one nonsquare search task."""

from bisect import bisect_left, bisect_right
from collections import Counter
from itertools import combinations, groupby
from math import gcd, isqrt, prod


from primes import primes_up_to


def smooth(lo, hi, primes, index=0, value=1):
    if index == len(primes):
        if lo <= value <= hi:
            yield value
        return
    p = primes[index]
    while value <= hi:
        yield from smooth(lo, hi, primes, index + 1, value)
        value *= p


def exponents(value, primes):
    """Return its exponent tuple, or None if there is a prime outside the set."""
    result = []
    for p in primes:
        n = 0
        while value % p == 0:
            value //= p
            n += 1
        result.append(n)
    return tuple(result) if value == 1 else None


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


def terminal(L, E, lo, hi, primes):
    """Solve E*k²=1+d²+e² with d≤e divisors of L and lo≤kL/e≤hi."""
    # k is odd and the only square residues modulo 8 are 0, 1 and 4.
    if E % 8 not in (1, 2, 3, 5, 6):
        return
    if hi - lo <= 64:
        matches = []
        for b in range(lo, hi+1):
            g = gcd(b, L)
            k, e = b // g, L // g
            if gcd(k, 10) != 1:
                continue
            d2 = E*k*k - e*e - 1
            if 36 <= d2 <= e*e:
                d = isqrt(d2)
                if d*d == d2 and L % d == 0:
                    powers = exponents(b, primes)
                    if powers is not None:
                        matches.append((powers, b, k, d, e))
    else:
        # For J = supp(k), both d and e divide the J-free part U of L.
        # Hence at least two target primes lie outside J, and E*k² ≤ 1+2U².
        # For fixed k, the residual and denominator bounds give an interval
        # of divisors e; d is then uniquely determined.
        Lpowers = exponents(L, primes)
        allowed = [p for p in primes if p not in (2, 5)]
        required = {p for p, a in zip(primes, Lpowers) if a == 0}
        matches = []
        for count in range(min(len(allowed), len(primes)-2) + 1):
            for support in combinations(allowed, count):
                if not required.issubset(support):
                    continue
                U = prod(p**a for p, a in zip(primes, Lpowers) if p not in support)
                base = prod(support)
                max_k = isqrt((1+2*U*U)//E)
                if max_k < base:
                    continue
                divisors = [1]
                for p, a in zip(primes, Lpowers):
                    if a and p not in support:
                        divisors = [d*p**i for d in divisors for i in range(a+1)]
                divisors = sorted(d for d in divisors if sum(d%p == 0 for p in primes) >= 2)
                for power in smooth(1, max_k//base, support):
                    k = base*power
                    residue = E*k*k - 1
                    if residue < 72:
                        continue
                    e_min = max((k*L+hi-1)//hi, isqrt((residue-1)//2)+1)
                    e_max = min(k*L//lo, isqrt(residue-36))
                    for index in range(bisect_left(divisors,e_min), bisect_right(divisors,e_max)):
                        e = divisors[index]
                        d2 = residue-e*e
                        d = isqrt(d2)
                        if d*d != d2 or U%d:
                            continue
                        b = k*(L//e)
                        matches.append((exponents(b, primes), b, k, d, e))
    for _, b, k, d, e in sorted(matches):
        yield b, k, d, e


def ternary_tiles(denominators, slots):
    # A closed multiplicity n left after prior obligations needs exactly
    # 2(n mod 3) occurrences at 2z. Three mixed tiles may be exchanged
    # for one equal tile at z and two equal tiles at 2z.
    counts = Counter(z for z in denominators if z % 3)
    last = denominators[-1]
    owed = Counter()
    for z, n in sorted(counts.items()):
        if z == last and slots:
            break
        remainder = n - owed.pop(z, 0)
        if remainder < 0:
            return False
        if remainder % 3:
            owed[2*z] += 2*(remainder % 3)
    if any(z < last for z in owed):
        return False
    return sum(owed.values()) - min(owed[last], counts[last]) <= slots if slots else not owed


def family_bound(prefix, primes, q, L, E, slots, high):
    """Use necessary heights from uncovered ordered three-prime dimensions."""
    if len(primes) != 3:
        return high
    choices = []
    for p, threshold in zip(primes, (6, 5, 4)):
        floor = p**((threshold-int(q % p == 0)+1)//2)
        if slots == 2 and p in (2, 5) and L % floor:
            return 0
        seen, rest = 1, L
        while rest % p == 0:
            seen *= p
            rest //= p
        options = []
        for growing in (False, True):
            if not growing and seen < floor:
                continue
            power = max(floor, seen*p) if growing else seen
            known = 0 if growing else sum(b % seen == 0 for b in prefix)
            need = max(0, (3 if p in (2, 5) else 2)-known)
            if p == 2:
                need += (3-known-need) % 4
            if p == 3:
                need += (2-known-need) % 3
            if need > slots:
                continue
            denominator = ((prefix[-1]-1)//power+1)*power
            mass = need*((L*L+denominator*denominator-1)//(denominator*denominator))
            options.append((need, mass))
        choices.append(options)
    best = 0
    for a in choices[0]:
        for b in choices[1]:
            for c in choices[2]:
                reserved = a[0]+b[0]+c[0]
                if reserved > slots:
                    continue
                mass = 1+a[1]+b[1]+c[1]
                bound = high if E <= mass else min(high, isqrt((slots-reserved)*L*L//(E-mass)))
                best = max(best, bound)
    return best


def complete(q, b, m, rank, families=False):
    """Enumerate a complete nonsquare task with arbitrary-precision integers."""
    yield from _complete(q, b, m, rank, None, [b], families)


class InvalidPrefix(ValueError):
    pass


def complete_prefix(q, b, m, rank, primes, prefix, families=False):
    """Continue one native node on precisely its fixed support and prefix."""
    if (not primes or list(primes) != sorted(set(primes)) or any(p not in primes_up_to(rank) for p in primes)
            or not 1 <= len(prefix) <= rank-3 or prefix[0] != b
            or any(type(z) is not int or z < 1 for z in prefix) or list(prefix) != sorted(prefix)
            or any(exponents(z, primes) is None for z in prefix)):
        raise InvalidPrefix('invalid exact-prefix handoff')
    yield from _complete(q, b, m, rank, tuple(primes), list(prefix), families)


def _complete(q, b, m, rank, fixed_primes, prefix, families):
    if (any(type(x) is not int for x in (q, b, m, rank)) or rank < 2
            or not 2 <= q < rank or b < 1 or not 0 <= m < rank):
        raise ValueError('expected rank >= 2, 2 <= q < rank, b >= 1 and 0 <= m < rank')
    f = q*b*b
    all_primes = primes_up_to(rank)
    ramified = [p for p in all_primes if q % p == 0]
    if (f < 8 or f >= rank or any(q % (p*p) == 0 for p in ramified)
            or (b == 1 and len(ramified) == 1)):
        return
    sizes = [(p+1)//2 for p in ramified if p != 2]
    product = prod(sizes)
    fiber = max([(2 if q % 2 == 0 else 1)*product]
                + [(a+1)*(product//a) for a in sizes])
    must = tuple(p for p in all_primes if q*b % p == 0)
    limits = [f-1-(p-1)//2 for p in all_primes if p % 2 and f % p == 0 and f % (p*p)]
    top_limit = min(limits, default=rank)
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
        if f == 8 and prefix.count(b) > 2 and L % 4 == 0:
            return
        if not possible_packets(prefix, slots, primes, q, fiber):
            return
        if q == 3 and not ternary_tiles(prefix, slots):
            return
        first_run = prefix[-1] == b
        if first_run and (len(prefix) > top_limit or f == 8 and len(primes) <= 4 and len(prefix) > 4):
            return
        lo = max(prefix[-1], isqrt(L*L//E) + 1)
        hi = isqrt((36*slots+1)*L*L//(36*E))
        if families:
            hi = family_bound(prefix, primes, q, L, E, slots, hi)
        if first_run and m:
            if len(prefix) < m:
                hi = min(hi, b)
            else:
                lo = max(lo, b+1)
        if lo > hi:
            return
        if slots == 2:
            for z, k, d, e in terminal(L, E, lo, hi, primes):
                if gcd(k, 10*d*e) != 1:
                    continue
                M = k*L
                if any(M % p for p in primes):
                    continue
                row = [1, d, e] + [M // t for t in reversed(prefix)]
                if m and row[1:].count(row[-1]) != m:
                    continue
                if any(sum(y % p == 0 for p in primes) < 2 for y in row[1:]):
                    continue
                if any(sum(y % p != 0 for y in row[1:]) < 2 for p in primes):
                    continue
                denominators = prefix + [z, M//d]
                if f == 8 and denominators.count(b) > 2 and (M % 4 == 0 or len(primes) <= 4 and denominators.count(b) > 4):
                    continue
                if q == 3 and not ternary_tiles(denominators, 0):
                    continue
                if possible_packets(denominators, 0, primes, q, fiber):
                    yield row
            return
        for z in smooth(lo, hi, primes):
            g = gcd(z, L)
            k, e = z//g, L//g
            prefix.append(z)
            yield from visit(k*L, E*k*k-e*e, slots-1, primes)
            prefix.pop()

    def supports(start, chosen, budget):
        if len(chosen) >= (4 if f == 8 else 3):
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
    budget = rank//fiber - 1 - sum(p-1 for p in required if p != 2 and q % p)
    if budget < 0:
        return
    optional = [p for p in all_primes if p not in required] if fixed_primes is None else []
    for primes in supports(0, required, budget):
        yield from visit(L, E, rank-1-len(prefix), primes)
