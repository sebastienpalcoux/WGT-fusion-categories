"""Integral fusion categories without a braiding: ordinary rank 6 and odd rank 9.

Run next to families.py, which contains the divisor-closed dimension tests.
No list of codegrees, rings, or degree types is read.
"""
from collections import Counter
from functools import lru_cache
from math import gcd, lcm

from families import factors, family


@lru_cache(None)
def divisors(n):
    values = [1]
    for prime, exponent in factors(n):
        values = [d * prime**i for d in values for i in range(exponent + 1)]
    return values


def profiles(rank, odd):
    weight = 2 if odd else 1

    def search(p, q, left, common, previous, prefix):
        if left == 1:
            for k in divisors(q):
                if (1 + weight * k) % p:
                    continue
                dimension = q * ((1 + weight * k) // p)
                x = dimension // k
                if x >= previous and dimension % common == 0 and (not odd or dimension % 2):
                    yield (*prefix, x, dimension)
            return
        low = max(previous, weight * q // p + 1)
        high = (left * weight + 1) * q // p
        for x in range(low, high + 1):
            if odd and x % 2 == 0:
                continue
            a, b = p * x - weight * q, q * x
            g = gcd(a, b)
            yield from search(a // g, b // g, left - 1, lcm(common, x), x, (*prefix, x))

    yield from search(1, 1, (rank - 1) // weight, 1, weight + 1, ())


def center_covered(dimension):
    return (dimension**2 < (6087156075 if dimension % 2 else 29160000) or
            family(tuple((p, 2 * e) for p, e in factors(dimension))))


def main():
    for rank in (3, 5, 7, 9):
        tuples = list(profiles(rank, True))
        assert len(set(tuples)) == len(tuples)
        assert all(center_covered(fs[-1]) for fs in tuples)
        print("odd", rank, len(tuples), "all centers covered")
    for rank in range(3, 6):
        assert all(center_covered(fs[-1]) for fs in profiles(rank, False))
    counts = Counter()
    retained = []
    for fs in profiles(6, False):
        counts['all'] += 1
        if {2, 3, 7} <= set(fs):
            counts['greedy'] += 1
            continue
        if center_covered(fs[-1]):
            counts['center'] += 1
            continue
        if fs[0] == 2 and fs[-1] % 4 != 2:
            counts['stabilizer'] += 1
            continue
        retained.append(fs)
    assert not retained
    assert counts == {'all': 1568, 'greedy': 232, 'center': 1310, 'stabilizer': 26}
    print(dict(sorted(counts.items())))


if __name__ == '__main__':
    main()
