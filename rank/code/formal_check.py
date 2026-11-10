"""Eliminate the ordinary formal-codegree tuples output by formal_sieve.cpp."""
import sys
from collections import Counter
from itertools import combinations_with_replacement
from math import isqrt, prod
from pathlib import Path

def divisors(n):
    small = [d for d in range(1, isqrt(n) + 1) if n % d == 0]
    return {k for d in small for k in (d, n // d)}

def square(n):
    return isqrt(n) ** 2 == n

def rational_case(dimension):
    # One unit, three consecutive nonunit degrees, and four more.
    valid = {d for d in divisors(dimension) if 1 < d and d * d < dimension}
    pairs = {x*x + y*y for x, y in combinations_with_replacement(valid, 2)}
    for x in sorted(valid):
        if x + 1 not in valid or x + 2 not in valid:
            continue
        target = dimension - 1 - x*x - (x+1)**2 - (x+2)**2
        if any(target - s in pairs for s in pairs):
            return True
    return False

def cubic_case(dimension, codegrees):
    # Three units, five nonunits of degrees 3a. At f=5 two a differ by one.
    if dimension % 9 != 3:
        return False
    if 5 not in codegrees:
        return True
    target = (dimension - 3) // 9
    valid = {d for d in divisors(dimension // 3) if d*d <= target}
    pairs = {x*x + y*y for x, y in combinations_with_replacement(valid, 2)}
    for a in valid:
        if a + 1 not in valid:
            continue
        remaining = target - a*a - (a+1)**2
        if any(remaining - b*b in pairs for b in valid):
            return True
    return False

def verify(path):
    counts = Counter()
    rows = set()
    survivors = []
    for line in Path(path).read_text().splitlines():
        rank, *fs = map(int, line.split())
        D = fs[-1]
        assert len(fs) == rank and fs == sorted(fs)
        assert all(D % x == 0 for x in fs) and sum(D // x for x in fs) == D
        assert tuple(fs) not in rows
        rows.add(tuple(fs))
        counts['profiles'] += 1
        multiplicities = Counter(fs)
        determinant = prod(fs)
        if len(multiplicities) == rank:
            if square(determinant):
                survivors.append(fs)
            else:
                counts['distinct and nonsquare'] += 1
            continue
        # There are no other multiplicity patterns in the computed file.
        assert fs[:2] == [3, 3] and len(multiplicities) == rank - 1 and rank == 8
        if square(determinant):
            counts['rational degree tests'] += 1
            if rational_case(D):
                survivors.append(fs)
        elif determinant % 3 == 0 and square(determinant // 3):
            counts['cubic degree tests'] += 1
            if cubic_case(D, fs):
                survivors.append(fs)
        else:
            counts['other square classes'] += 1
    print(dict(sorted(counts.items())))
    print('Distinct dimensions:', len({fs[-1] for fs in rows}))
    print('Survivors:', survivors)
    assert len(rows) == 52 and not survivors

if __name__ == '__main__':
    verify(sys.argv[1] if len(sys.argv) > 1 else 'formal_residual.txt')
