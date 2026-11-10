"""Prefix feasibility for the two tile sizes in eligible odd defect-one blocks.

Inputs are denominators z=M/d, not degrees.  The caller supplies an odd prime
p dividing the squarefree kernel of D=q M^2. The vacuum has defect >=2.
Use a new obligation dictionary for each p and copy it at a search fork.
"""
from math import isqrt


def eligible(p):
    """p is assumed prime."""
    if p <= 2 or (p-1) % 8 == 0:
        return False
    odd = p-1
    while odd % 2 == 0:
        odd //= 2
    return all(odd % (a*a) for a in range(3, isqrt(odd)+1, 2))


def close_run(p, pending, z, count, next_z, available):
    """Consume one entire run. Return updated obligations or None if impossible.

    pending maps future denominators to required label counts. next_z is
    the next denominator or None after the final run. available is the count
    of remaining nonunit slots, including the run that starts at next_z.
    """
    upcoming = pending.copy()
    if z % p:
        free = count - upcoming.pop(z, 0)
        if free < 0:
            return None
        outgoing = (-2*free) % p
        if (p-1)//2 * outgoing > free:
            return None
        if outgoing:
            upcoming[2*z] = upcoming.get(2*z, 0) + 2*outgoing
    if any(next_z is None or b < next_z for b in upcoming):
        return None
    return upcoming if sum(upcoming.values()) <= available else None


def can_partition(p, runs):
    """Exact for the indicated necessary tile relaxation on complete data."""
    if not eligible(p):
        raise ValueError('the two-tile rule does not apply to this prime')
    data = sorted((z, m) for z, m in runs.items() if m)
    if any(z <= 0 or m < 0 for z, m in data):
        raise ValueError('denominators must be positive and counts nonnegative')
    pending = {}
    available = sum(m for _, m in data)
    for i, (z, count) in enumerate(data):
        available -= count
        next_z = data[i+1][0] if i+1 < len(data) else None
        pending = close_run(p, pending, z, count, next_z, available)
        if pending is None:
            return False
    return True


if __name__ == '__main__':
    from functools import lru_cache
    from itertools import product

    def direct(p, counts):
        # All tilings of z=1,2,4,8, independent of the greedy recurrence.
        w = (p-1)//2
        @lru_cache(None)
        def solve(st):
            if not any(st):
                return True
            i = next(i for i,x in enumerate(st) if x)
            if st[i] >= p:
                nxt = list(st)
                nxt[i] -= p
                if solve(tuple(nxt)):
                    return True
            if i+1 < len(st) and st[i] >= w and st[i+1] >= 2:
                nxt = list(st)
                nxt[i] -= w
                nxt[i+1] -= 2
                return solve(tuple(nxt))
            return False
        return solve(tuple(counts))

    assert [p for p in (3,5,7,11,13,17,19,23,29,31) if eligible(p)] == [3,5,7,11,13,23,29,31]
    for p in (3,5,7,11,13):
        for counts in product(range(5), repeat=4):
            assert can_partition(p, dict(zip((1,2,4,8),counts))) == direct(p,counts)
        w = (p-1)//2
        for counts in ((p,w,2,0),(p*w,2*p,0,0),(w,2+w,2,0),(p,p,p,p)):
            assert can_partition(p, dict(zip((1,2,4,8),counts))) == direct(p,counts)
    assert can_partition(5,{1:2,2:2,3:5,5:1})
    assert not can_partition(5,{1:1,2:4})
    print('two-tile prefix checks passed')
