"""Exact necessary mass and multiplicity bounds at odd defect one."""
from collections import Counter
from math import lcm
from families import factors


def defect_cap(codegree):
    return codegree - 1 - sum((p - 1)//2 for p, e in factors(codegree) if p % 2 and e == 1)


def first_excess(degrees):
    total = sum(x*x for x in degrees)
    for x, m in sorted(Counter(degrees[1:]).items(), reverse=True):
        f = total // (x*x)
        # Even the strongest odd-prime saving leaves a cap >= (f-1)/2.
        if f > 2*len(degrees)+1:
            break
        if m > defect_cap(f):
            return f
    return None


def can_pack(q, p, runs):
    """Necessary test on *closed* multiplicities among nonunit denominators.

    More closed runs strengthen the test. Keeping the six smallest applicable
    denominators bounds the work and can only retain additional candidates.
    In a complete profile every nonunit run is closed.
    """
    selected = sorted((b, m) for b, m in runs if b % p)[:6]
    if not selected:
        return True
    scale = lcm(*(b for b, _ in selected))
    weight = [(scale//b)**2 for b, _ in selected]
    total = sum(u*m for u, (_, m) in zip(weight, selected))
    bound = (q//p)*scale**2
    low, high = total - bound + 1, (p-1)*bound - 1
    step = (p-1)//2
    choices = [(u*step, m//step) for u, (_, m) in zip(weight, selected)]
    suffix = [0]*(len(choices)+1)
    for i in reversed(range(len(choices))):
        suffix[i] = suffix[i+1] + choices[i][0]*choices[i][1]

    def feasible(i, value):
        if low <= value <= high:
            return True
        if value > high or value+suffix[i] < low or i == len(choices):
            return False
        u, count = choices[i]
        lower = max(0, -((value+suffix[i+1]-low)//u))
        upper = min(count, (high-value)//u)
        return any(feasible(i+1, value+j*u) for j in range(lower, upper+1))
    return feasible(0, 0)


def joint_obstruction(degrees, profile):
    q = 1
    M = 1
    for p, e in profile:
        M *= p**(e//2)
        if e % 2:
            q *= p
    runs = Counter(M//x for x in degrees[1:]).items()
    return next((p for p, e in profile if p % 2 and e % 2 and not can_pack(q, p, runs)), None)
