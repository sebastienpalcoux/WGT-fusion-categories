"""Verify the four-pair intersection counts in the quartic investigation.

The input describes only sets of exponents of roots of unity. No dimension,
category, or list of fusion rules is enumerated. Run with Python 3.10+.
"""
from itertools import combinations
from math import comb


def intersection_counts(prime, minimum):
    """Count normalized four-pair sets and their admissible difference masks."""
    h = (prime - 1) // 2

    # Enumerate F_prime^*/{+1,-1} in multiplicative order.
    for base in range(2, prime):
        representatives = [pow(base, j, prime) for j in range(h)]
        if len({min(x, prime - x) for x in representatives}) == h:
            break
    index = [0] * prime
    for i, x in enumerate(representatives):
        index[x] = index[-x] = 1 << i

    # Nonzero opposite pairs among the differences of {0,+/-x_i}.
    single = [index[x] | index[2 * x % prime] for x in representatives]
    pair = [[index[(x - y) % prime] | index[(x + y) % prime]
             for y in representatives] for x in representatives]
    mask = (1 << h) - 1
    large = survivors = 0
    for i, j, k in combinations(range(1, h), 3):
        # The first pair is +/-1, obtained by scaling.
        differences = (single[0] | single[i] | single[j] | single[k] |
                       pair[0][i] | pair[0][j] | pair[0][k] |
                       pair[i][j] | pair[i][k] | pair[j][k])
        if differences.bit_count() < minimum:
            continue
        large += 1
        for shift in range(1, h // 2 + 1):
            dilate = ((differences << shift) & mask) | (differences >> (h - shift))
            if (differences & dilate).bit_count() < 4:
                break
        else:
            survivors += 1
    return comb(h - 1, 3), large, survivors


if __name__ == "__main__":
    cases = {
        (173, 19): (98770, 52012, 0),
        (179, 20): (109736, 50376, 0),
        (181, 20): (113564, 52658, 0),
        (191, 20): (134044, 65108, 0),
        (193, 20): (138415, 67813, 0),
    }
    for key, expected in cases.items():
        result = intersection_counts(*key)
        assert result == expected, (key, result, expected)
        print(*key, *result)
