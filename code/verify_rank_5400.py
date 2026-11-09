"""Exact character-block check for ranks 10 and 12 in dimension 5400.

Run with Python 3.8 or later. No third-party packages are required.
Blocks (w,k) record central degree and multiplicity in induction of the unit.
The unit block is omitted throughout.
"""

from collections import Counter
from math import isqrt


DIM = 5400
DEGREES = tuple(
    d for d in range(10, 901)
    if DIM % d == 0 and sum(d % p == 0 for p in (2, 3, 5)) >= 2
)


def degree_profiles(left=9, total=DIM - 1, last=len(DEGREES) - 1):
    """All unordered lists of nine admissible central degrees summing to 5399."""
    if left == 0:
        if total == 0:
            yield ()
        return
    for j in range(last, -1, -1):
        d = DEGREES[j]
        if left * d < total:
            break
        if d + (left - 1) * DEGREES[0] <= total:
            for rest in degree_profiles(left - 1, total - d, j):
                yield (d,) + rest


def possible_rows(blocks):
    """Use the two canonical trace identities at X and XX* to list rows."""
    def entries(d, w, k, cap):
        if k == 1:
            limit = min(d - 1, isqrt(min(cap, DIM // w - 1)))
            for s in range(-limit, limit + 1):
                yield s, s*s
        else:
            for t in range(cap + 1):
                limit = min(k*d - 1, isqrt(k*t))
                for s in range(-limit, limit + 1):
                    yield s, t

    def complete(d, i, linear, square, chosen):
        w, k = blocks[i]
        if i == len(blocks) - 1:
            if linear % w or square % w:
                return
            s, t = linear // w, square // w
            if abs(s) >= k*d or s*s > k*t:
                return
            if k == 1 and (t != s*s or t >= DIM // w):
                return
            yield d, chosen + ((s, t),)
            return
        for s, t in entries(d, w, k, square // w):
            yield from complete(d, i+1, linear-w*s, square-w*t, chosen+((s, t),))

    for d in range(3, isqrt(DIM-1) + 1):
        yield from complete(d, 0, -d, DIM-d*d, ())


def moment_data(blocks, rows, double):
    """Orthogonality targets after removing the unit, and possible row vectors."""
    target = [DIM - 1]
    for i, (w, k) in enumerate(blocks):
        target.extend(((k + 2*(i == double))*(DIM//w) - k, -k))

    # Distinct character rows can have the same two kinds of moments.
    # One realizable singleton suffices when taking necessary conditions.
    vectors = {}
    for d, block_values in rows:
        moments = [d*d]
        for s, t in block_values:
            moments.extend((t, d*s))
        key = tuple(moments)
        group = 2 if any((s*s-t) % 2 for s, t in block_values) else 1
        vectors[key] = min(group, vectors.get(key, group))
    return tuple(target), sorted(vectors.items())


def attainable(vectors, target, slots):
    """Exact search for 'slots' rows whose vector sum is the target.

    For pruning, precompute the possible sums of each coordinate separately.
    A bit at position n represents the attainable integer n. Shifting every row
    coordinate by its minimum makes all the bit positions nonnegative.
    """
    if not vectors:
        return False
    minima = tuple(map(min, zip(*(v for v, group in vectors))))
    target = tuple(t - slots*m for t, m in zip(target, minima))
    if min(target) < 0:
        return False
    rows = [
        (tuple(v[j]-minima[j] for j in range(len(target))), group)
        for v, group in vectors
    ]

    # suffix[j][i][h] encodes the sums attainable in coordinate j with h
    # labels drawn from rows i, i+1, ... . A row may be reused.
    ranges = []
    for j, goal in enumerate(target):
        mask = (1 << (goal+1)) - 1
        suffix = [None] * (len(rows)+1)
        suffix[-1] = [1]+[0]*slots
        for i in reversed(range(len(rows))):
            row, group = rows[i]
            step = group*row[j]
            sums = suffix[i+1].copy()
            for h in range(group, slots+1):
                sums[h] |= (sums[h-group] << step) & mask
            suffix[i] = sums
        if not (suffix[0][slots] >> goal) & 1:
            return False
        ranges.append(suffix)

    def search(start, left, need):
        if left == 0:
            return not any(need)
        for j, n in enumerate(need):
            if n < 0 or not ((ranges[j][start][left] >> n) & 1):
                return False
        for i in range(start, len(rows)):
            row, group = rows[i]
            if left*row[0] > need[0]:
                break  # The rows are sorted by squared dimension.
            if group <= left:
                rest = tuple(n - group*a for n, a in zip(need, row))
                if search(i, left-group, rest):
                    return True
        return False

    return search(0, slots, target)


def main():
    profiles = placements = survivors = 0
    for degrees in degree_profiles():
        profiles += 1
        # Finish with a repeated block when possible, so its two trace entries
        # are determined by the remaining linear and squared sums.
        blocks = sorted(Counter(degrees).items(), key=lambda p: (p[1] > 1, -p[0]))
        rows = list(possible_rows(blocks))
        cases = [None] + [i for i, (w, k) in enumerate(blocks) if k >= 2]
        for double in cases:
            placements += double is not None
            rank = 10 if double is None else 12
            target, vectors = moment_data(blocks, rows, double)
            if attainable(vectors, target, rank-1):
                print("Survivor", degrees, "double block", double)
                survivors += 1
    print(f"Profiles {profiles}, noncommutative placements {placements}, survivors {survivors}")
    if (profiles, placements, survivors) != (26, 36, 0):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
