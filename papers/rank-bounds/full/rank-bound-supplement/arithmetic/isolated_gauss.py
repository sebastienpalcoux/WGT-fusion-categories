"""Gauss-sum check when no twist is active at two distinct primes."""
from fractions import Fraction
from functools import lru_cache
from math import ceil, gcd, lcm, prod

from edge_budget import isolated_minimum


def gauss_possible(rank, support, multiplicities, active, nines, proper=0, high=(), neutral=()):
    """Test whether any completion meets the separated Gauss and mass bounds."""
    support = tuple(support)
    other = tuple(p for p in support if p not in multiplicities or p in high)
    base = prod(p if p in multiplicities and p not in high else p * p for p in support)
    room = rank - 1 - sum(multiplicities.values()) - nines - proper - sum(n for n,b in neutral)
    charge = 1 - sum((Fraction(a, p * (p - 1)) for p, a in active.items()), Fraction())
    quota = tuple(isolated_minimum(p) for p in other)

    # Retain an upper bound on vacuum mass in integer units. Rounding
    # every positive term upward can only permit additional completions.
    scale = 1 << 75

    @lru_cache(None)
    def upper_weight(b):
        return (scale + b * b - 1) // (b * b)

    free_weight = upper_weight(15)
    initial_mass = (upper_weight(base) + nines * upper_weight(9)
                    + sum(m * upper_weight(p) for p, m in multiplicities.items())
                    + proper * upper_weight(33)
                    + sum(n*upper_weight(b) for n,b in neutral) + room * free_weight)
    high_min = {}
    for i, p in enumerate(other):
        if p in high:
            a = active[p]
            projected = (multiplicities[p] - a + Fraction((p + 1) * a, (p - 1)**2)) / p**2
            high_min[i] = a + ceil((projected - Fraction(1, base**2)) * p**4) - proper

    @lru_cache(None)
    def next_levels(last):
        result = []
        for p in other:
            b = p * p
            while b <= last:
                b *= p
            result.append(b)
        return tuple(result)

    def visit(last, slots, num, den, mass, multiple, credits):
        # num/den is the exact signed residual; reduction is unnecessary.
        if slots < 0 or mass < scale or num == 0:
            return False
        missing = [max(0, q - c) for q, c in zip(quota, credits)]
        total_missing = sum(missing)
        if total_missing > slots:
            return False
        for i, minimum in high_min.items():
            n = max(0, minimum - credits[i])
            p = other[i]
            if n + total_missing - missing[i] > slots or mass - n * (free_weight - upper_weight(p*p)) < scale:
                return False

        levels = next_levels(last)
        loss = sum(n * (free_weight - upper_weight(b)) for n, b in zip(missing, levels))
        if mass - loss < scale:
            return False
        if total_missing == 0 and den % abs(num) == 0:
            M = den // abs(num)
            if M % multiple == 0:
                return True

        for i, p in enumerate(other):
            limit = slots // (p - 1)
            if limit == 0:
                continue
            b = levels[i]
            while True:
                required = lcm(multiple, b)
                left = num * required
                square = b * b
                # -1/required <= num/den <= 1/required + 3*slots/(2*b^2).
                if left < -den or 2 * (left - den) * square > 3 * slots * den * required:
                    break

                # The current coefficient c decreases the residual by c/h.
                # At the following level, apply the same two inequalities
                # with the old slot count, which is an upper bound.
                h = square // p
                common = lcm(den, h)
                numerator = num * (common // den)
                unit = common // h
                next_square = min(next_levels(b))**2
                lower_num = (2 * (left - den) * next_square - 3 * slots * den * required) * h
                lower_den = 2 * den * required * next_square
                first = max(1, -(-lower_num // lower_den))
                last_c = min(limit, (left + den) * h // (den * required))

                # A largest primitive-square-class imbalance would need at
                # least p^2 packets. Higher phase orders contribute even c.
                step = 2 if limit < p*p else 1
                first += (-first) % step
                for c in range(first, last_c + 1, step):
                    used = (p - 1) * c
                    next_multiple = required
                    if limit < p*p:
                        next_multiple = lcm(required, h // gcd(h, c))
                    new_credits = list(credits)
                    new_credits[i] += p * c
                    new_mass = mass - used * (free_weight - upper_weight(b))
                    if visit(b, slots - used, numerator - c * unit, common,
                             new_mass, next_multiple, new_credits):
                        return True
                b *= p
        return False

    for c in range(nines // 2 + 1) if 3 in other else (0,):
        residual = charge - Fraction(c, 27)
        credits = [active.get(p, 0) if p in high else 0 for p in other]
        if 3 in other:
            credits[other.index(3)] = 3 * c
        if visit(9, room, residual.numerator, residual.denominator, initial_mass, base, credits):
            return True
    return False
