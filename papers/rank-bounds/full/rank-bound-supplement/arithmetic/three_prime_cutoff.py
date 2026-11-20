"""Exclude square dimensions on three primes at rank at most 44."""
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import isqrt, lcm

from full_census import census
from row_quota import row_quota

RANK = 44
CUTOFF = 65536
ODD_PRIMES = [p for p in range(3, RANK, 2) if all(p % j for j in range(3, isqrt(p)+1, 2))]


@cache
def omissions(p, fours):
    hits = row_quota(p, fours) if fours >= 4 else 0
    if p == 3:
        return hits, hits + (2-hits) % 3
    if p in (5, 7) and fours >= 4:
        return hits, census(p, fours)
    return hits, max(hits, 3 if p == 5 else 2)


def enough_mass(groups, ordinary):
    available = RANK - sum(count for count, _ in groups)
    if available < 0:
        return False
    mass = sum(Fraction(count, floor*floor) for count, floor in groups)
    return mass + Fraction(available, ordinary*ordinary) >= 1


def surviving_modes(p, t, m):
    ordinary = 6 if p == 3 else 8
    # A=4 records exact dyadic valuation two, A=8 covers all greater valuations.
    for A in ([4, 8] if m < 4 else [8]):
        odd = 3 - (m if A == 4 else 0)
        if odd <= 0:
            odd += 4   # at least one principal-t witness is additional to the fours
        odd_floor = 4*p if A == 4 else A
        for mode, count, floor in [('dyadic', 4*(t-1), 2*t), ('odd', p*(t-1), p*t)]:
            if mode == 'odd' and m <= 3:
                groups = [(m, 4), (p*t, p*t)]      # may include the vacuum
            else:
                groups = [(1, A*p*p*t*t), (m, 4), (count, floor)]
                if m:
                    groups.extend([(1, A*t), (odd-1, odd_floor)])
                if m >= 4:
                    groups.extend((omissions(s, m)[0], 4*s*s) for s in (p, t))
            if enough_mass(groups, ordinary):
                yield mode, A

        if not enough_mass([(1, 900), (m, 4), (t, t*t)], ordinary):
            continue
        kp, np = omissions(p, m)
        kt, nt = omissions(t, m)
        for mode, pure in [('dense', t*(t-1)), ('nonminimal', 2*(t-1)), ('minimal', t)]:
            dense = mode == 'dense'
            T = t**(2 if dense else 3)
            groups = [(1, A*p*p*T), (m, 4), (1, A*(t if dense else t*t)),
                      (odd-1, odd_floor), (pure, t*t)]
            if m < 4:
                groups.extend([(1, p*p*t), (np-1, p*p), (1, p*T)])
                if not dense:
                    groups.append((nt-1, T))
            else:
                groups.extend([(kp, 4*p*p), (np-kp, p*p), (kt, 4*T)])
                if not dense:
                    groups.append((nt-kt, T))
            if enough_mass(groups, ordinary):
                yield mode, A


def first_large_exponent(p):
    exponent = 0
    while p**exponent <= CUTOFF:
        exponent += 1
    return exponent


A, B, C = map(first_large_exponent, (2, 3, 5))
SMOOTH = sorted((2**i*3**j*5**k, (i, j, k))
                for i, j, k in product(range(A), range(B), range(C))
                if 6 <= 2**i*3**j*5**k <= CUTOFF)
SCALE = lcm(4, *(d for d, _ in SMOOTH))**2


def scaled_upper(count, denominator):
    return -(-count*SCALE // (denominator*denominator))


@cache
def choices(a, b, c, omitted, divisor):
    exponents = (a, b, 2*c-1)
    floors = (2**a, 3**b, 5**(2*c-1))
    required = divisor if omitted < 0 else lcm(divisor, floors[omitted])
    intervals = []
    for d, powers in SMOOTH:
        if d % divisor or any(x > y for x, y in zip(powers, exponents)):
            continue
        missing = [i for i in range(3) if powers[i] == exponents[i]]
        if len(missing) > 1 or (omitted >= 0 and missing != [omitted]):
            continue
        mass = SCALE // (d*d)
        intervals.append((mass, mass))
    intervals.append((0, scaled_upper(1, max(CUTOFF+1, required))))
    return intervals


def interval_count(m, a, b, c):
    if a < A and b < B and pow(2, a, 125)*pow(3, b, 125) % 125 not in (1, 124):
        return 0, 0
    kp, np = omissions(3, m)
    kt, nt = omissions(5, m)
    vacant = RANK-m-1-5-3-np-nt
    requirements = [(vacant, -1, 1), (1, 0, 5), (2, 0, 1),
                    (kp, 1, 4), (np-kp, 1, 1), (kt, 2, 4), (nt-kt, 2, 1)]
    if vacant < 0:
        return 0, 0
    groups = [(n, choices(a, b, c, omitted, divisor))
              for n, omitted, divisor in requirements if n]
    groups.sort(key=lambda g: -g[1][0][1])
    upper = [0]*(len(groups)+1)
    for i in reversed(range(len(groups))):
        n, values = groups[i]
        upper[i] = upper[i+1]+n*values[0][1]
    low = m*SCALE//16 + (5*SCALE//5**(2*c) if c < C else 0)
    high = (m*SCALE//16 + scaled_upper(5, 5**c)
            + scaled_upper(1, 2**a*3**b*5**(2*c-1)))
    visited = survived = 0

    def visit(i, count, start, lower, higher):
        nonlocal visited, survived
        visited += 1
        if i == len(groups):
            survived += lower < SCALE <= higher
            return
        values = groups[i][1]
        for index in range(start, len(values)):
            lo, hi = values[index]
            if higher + count*hi + upper[i+1] < SCALE:
                break
            if lower + lo >= SCALE:
                continue
            if count > 1:
                visit(i, count-1, index, lower+lo, higher+hi)
            else:
                next_count = groups[i+1][0] if i+1 < len(groups) else 0
                visit(i+1, next_count, 0, lower+lo, higher+hi)

    if high+upper[0] >= SCALE:
        visit(0, groups[0][0], 0, low, high)
    return visited, survived


def main():
    profiles = [(p, t, m, mode, floor) for p, t in combinations(ODD_PRIMES, 2)
                for m in range(14) for mode, floor in surviving_modes(p, t, m)]
    assert all((p, t, mode, floor) == (3, 5, 'minimal', 8)
               for p, t, _, mode, floor in profiles)
    visits = residual = 0
    for _, _, m, _, _ in profiles:
        for a, b, c in product(range(3, A+1), range(2, B+1), range(2, C+1)):
            n, s = interval_count(m, a, b, c)
            visits += n
            residual += s
    assert residual == 0
    print(f'Three-prime supports: {len(profiles)} profiles, {visits} interval nodes, {residual} survivors.')


if __name__ == '__main__':
    main()
