"""Exact mass and Gauss certificates for separated odd primes through rank44."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import prod, isqrt
import json

from phase_deficits import sparse_budget


# Each odd dimension prime needs at least p-1 active positions in addition to the vacuum.
PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def alternatives(p):
    """Reserved labels and denominator floor, then proper labels and their floor."""
    yield p, p*p, 0, 1
    if p in (3, 5):
        yield 4*(p - 1), 2*p, 0, 1
    if p == 5:
        yield 5, 5, 4, 20
    elif p > 5:
        yield p, p, 3, 2*p


def cases(rank):
    for primes in combinations(PRIMES, 4):
        for choices in product(*(alternatives(p) for p in primes)):
            room = rank - 1 - sum(count for count, _, _, _ in choices)
            if room < 0:
                continue
            fixed = Fraction(1, 2310**2) + sum(Fraction(n, b*b) for n, b, _, _ in choices)
            for fours in range(min(13, room) + 1):
                witnesses = [(n, b) for _, _, n, b in choices if n]
                if fours >= 9:
                    witnesses.append((3, 8))
                # The witness sets may overlap. Sort their union by decreasing denominator.
                size = max((n for n, _ in witnesses), default=0)
                floors = [max(b for n, b in witnesses if n >= j) for j in range(1, size + 1)]
                remaining = room - fours - size
                if remaining < 0:
                    continue
                mass = fixed + Fraction(fours, 16) + sum(Fraction(1, b*b) for b in floors)
                if mass + Fraction(remaining, 36) >= 1:
                    yield fours, primes, choices


def gauss_errors(rank, profiles):
    for fours, primes, choices in profiles:
        assert all(n == p or b == p for p, (n, b, _, _) in zip(primes, choices))
        present = tuple(p for p, (_, b, _, _) in zip(primes, choices) if p == b)
        assert present
        absent = tuple(p for p in primes if p not in present)
        required = 4 * prod(p if p in present else p*p for p in primes)
        proper = max(n for _, _, n, _ in choices)
        packets = [(1,) if p == 5 else range(1, p) for p in present]
        for counts in product(*packets):
            classes = [p if j == 1 else (p - 1)*j for p, j in zip(present, counts)]
            budget = rank - 1 - fours - sum(classes) - proper
            if budget < sum(p - 1 for p in absent):
                continue
            available_mass = 1 - Fraction(fours, 16) - sum(
                Fraction(n, p*p) for n, p in zip(classes, present))
            for trace in range(-fours, fours + 1):
                residual = 1 - Fraction(fours - trace, 16) - sum(
                    Fraction(j, p) for j, p in zip(counts, present))
                if -Fraction(1, required) < residual <= 2*available_mass + Fraction(1, required):
                    yield from sparse_budget(residual, required, absent, budget)



CUTOFF = 256


@lru_cache(None)
def square_rows(length, norm):
    """All nonnegative integer vectors with the given squared norm."""
    if length == 0:
        return [()] if norm == 0 else []
    return [(a,) + tail for a in range(isqrt(norm) + 1)
            for tail in square_rows(length - 1, norm - a*a)]


def coefficient_tables(present, slots):
    """Absolute proper Weyl coefficients, up to simultaneous slot permutations."""
    tables = {((),) * slots}
    for p in present:
        new_tables = set()
        for row in square_rows(slots, p - 1):
            if sum(a != 0 for a in row) < 3:
                continue
            for table in tables:
                if any(sum(a*b for a, b in zip(row, column)) % 2 == 0
                       for column in zip(*table)):
                    continue
                new_tables.add(tuple(sorted(old + (a,) for old, a in zip(table, row))))
        tables = new_tables
    return tables


def smooth(b, primes):
    for p in primes:
        while b % p == 0:
            b //= p
    return b == 1


def coefficient(b, p):
    if b % p:
        return 0
    value = pow(b // p, -1, p)
    return value - p if value > p // 2 else value


def mixed_counts(rank, fours, primes, choices):
    present = tuple(p for p, (_, b, _, _) in zip(primes, choices) if b == p)
    absent = tuple(p for p, (_, b, _, _) in zip(primes, choices) if b == p*p)
    assert primes[0] == 3 and choices[0][:2] == (8, 6)
    assert len(present) + len(absent) + 1 == len(primes)
    slots = rank - 1 - fours - sum(n for n, _, _, _ in choices)
    proper = max(n for _, _, n, _ in choices)
    # Five has j=1. For every other present prime, an additional packet
    # costs p places while the proper witnesses occupy at least `proper`.
    for p in present:
        if p > 5:
            for j in range(2, p - 1):
                if (p - 1) % j == 0:
                    assert slots - p*(j - 1) < proper

    minimum_M = 4 * prod(p*p if p in absent else p for p in primes)
    target = 1 - Fraction(fours, 16) - sum(Fraction(1, p) for p in present)
    zero = (0,) * len(present)
    small = range(6, CUTOFF + 1)
    mixed = tuple(b for b in small if b % 6 == 0 and smooth(b, (2, 3)))
    fixed = [(mixed, zero)] * 8
    for p in absent:
        pure = tuple(b for b in small if b >= p*p and smooth(b, (p,)))
        fixed.extend([(pure, zero)] * p)

    candidates = {}
    for b in small:
        if b in primes or any(b % (p*p) == 0 for p in present):
            continue
        if smooth(b, (2,) + primes):
            key = tuple(abs(coefficient(b, p)) for p in present)
            candidates.setdefault(key, []).append(b)

    mass_matches = survivors = 0
    tail = CUTOFF + 1
    for table in coefficient_tables(present, slots):
        pools = fixed + [(tuple(candidates.get(row, ())), row) for row in table]
        pools.sort(key=lambda item: (min(item[0], default=tail), item))
        # Every slot also admits a denominator above the cutoff. Its mass
        # interval includes zero, so unused slots cover all smaller ranks.
        upper = [Fraction(0)] * (len(pools) + 1)
        upper[-1] = Fraction(1, minimum_M**2)
        for i in reversed(range(len(pools))):
            upper[i] = upper[i+1] + Fraction(1, min(pools[i][0], default=tail)**2)
        chosen = []

        def visit(i, previous, mass, error):
            nonlocal mass_matches, survivors
            if mass >= target or mass + error + upper[i] < target:
                return
            if i == len(pools):
                mass_matches += 1
                for j, p in enumerate(present):
                    total = sum(Fraction(coefficient(b, p), b) for b, row in chosen
                                if b is not None and row[j])
                    unknown = sum(row[j] for b, row in chosen if b is None)
                    if abs(total) > Fraction(1, minimum_M) + Fraction(unknown, tail):
                        return
                survivors += 1
                return
            values, row = pools[i]
            start = previous if i and pools[i-1] == pools[i] else 0
            for k in range(start, len(values) + 1):
                b = values[k] if k < len(values) else None
                chosen.append((b, row))
                if b is None:
                    visit(i+1, k, mass, error + Fraction(1, tail**2))
                else:
                    visit(i+1, k, mass + Fraction(1, b*b), error)
                chosen.pop()

        visit(0, 0, Fraction(0), Fraction(0))
    return mass_matches, survivors


def verify():
    rank = 44
    rows = list(cases(rank))
    pure = [row for row in rows if all(n == p for p, (n, _, _, _) in zip(row[1], row[2]))]
    mixed = [row for row in rows if row not in pure]
    assert (len(rows), len(pure), len(mixed)) == (49, 38, 11)
    assert not list(gauss_errors(rank, pure))
    counts = [mixed_counts(rank, *profile) for profile in mixed]
    mass_matches = [mass for mass, _ in counts]
    survivors = [left for _, left in counts]
    assert mass_matches == [0, 0, 0, 0, 0, 0, 3, 1, 14, 0, 0]
    assert survivors == [0] * len(mixed)
    return dict(rank=rank, initial_profiles=len(rows), all_primary=len(pure),
                mixed=len(mixed), mixed_mass_matches=mass_matches,
                mixed_survivors=survivors, survivors=sum(survivors))


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
