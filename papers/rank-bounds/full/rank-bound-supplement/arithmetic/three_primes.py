"""Reciprocal-square bound for square dimension with three prime factors."""
from fractions import Fraction
from itertools import combinations
from math import isqrt

from row_quota import row_quota


def bounds(rank):
    odd_primes = [p for p in range(3, rank, 2)
                  if all(p % d for d in range(3, isqrt(p) + 1, 2))]
    for p, t in combinations(odd_primes, 2):
        modes = [('pure', t, t*t), ('dyadic', 4*(t - 1), 2*t),
                 ('joint', p*(t - 1), p*t)]
        for mode, active_count, active_floor in modes:
            for fours in range(14):
                groups = [(1, 900), (fours, 4)]
                if fours <= 3 and mode == 'joint':
                    # The full summand can contain the vacuum.
                    groups = [(fours, 4), (p*t, p*t)]
                else:
                    groups.append((active_count, active_floor))

                if fours <= 3 and mode == 'pure':
                    groups.extend([(1, p*p*t), (1, p*p), (1, t*t*p), (1, 4*t)])
                if fours >= 4:
                    p_free = row_quota(p, fours)
                    t_free = row_quota(t, fours)
                    groups.extend([(2, 8), (1, 8*t), (p_free, 4*p*p), (t_free, 4*t*t)])
                    if mode == 'pure':
                        if p == 3:
                            # The three-free nonunit count is two modulo three.
                            groups.append(((2 - p_free) % 3, p*p))
                        if fours == 10 and t == 5:
                            # A sixth five-free nonunit is forced by the row Gram bound.
                            groups.append((1, 25))
                remaining = rank - sum(count for count, _ in groups)
                if remaining >= 0:
                    mass = sum(Fraction(count, floor**2) for count, floor in groups)
                    yield (p, t, mode, fours), mass + Fraction(remaining, 36)


if __name__ == '__main__':
    key, maximum = max(bounds(41), key=lambda item: item[1])
    gram = min(Fraction((10-j)*36 + j*16, 16) + Fraction((-60 + 10*j)**2, 96)
               for j in range(11))
    assert gram == Fraction(355, 24)
    assert maximum == Fraction(215033, 216000) < 1
    print('Largest three-prime mass through rank 41:', maximum, key)
