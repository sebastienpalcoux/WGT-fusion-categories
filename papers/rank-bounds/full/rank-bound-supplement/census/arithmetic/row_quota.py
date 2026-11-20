"""Gram lower bound for the number or weighted sum of prime-free four-row hits."""
from fractions import Fraction
from functools import lru_cache
from math import ceil


@lru_cache(None)
def row_quota(prime, rows, free_parity=False):
    """Use even row sums, or minimize over their two possible common parities."""
    assert prime > 2 and 0 < rows < 16
    period = 2 * prime
    half = (prime - 1) // 2
    total = 0
    for residue in range(1, half + 1):
        # The two adjacent lifts suffice for each permitted common parity.
        lifts = [prime - residue]
        if free_parity:
            lifts.append(2 * prime - residue)
        elif lifts[0] % 2:
            lifts[0] += prime
        minima = []
        for upper in lifts:
            lower = upper - period
            for count in range(rows + 1):
                trace = (rows - count) * lower + count * upper
                norm = (rows - count) * lower**2 + count * upper**2
                minima.append((16 - rows) * norm + trace**2)
        total += ceil(Fraction(min(minima), 16 * (16 - rows)))
    return ceil(Fraction(total, sum(a*a for a in range(1, half + 1))))
