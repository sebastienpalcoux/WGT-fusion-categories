"""Necessary prime-power-row support test at a complete square-dimension type."""
from functools import reduce
from math import gcd, isqrt, lcm
from collections import Counter


def row_cover(degrees, top_base=2):
    """For D=M² and top denominator a a positive power of the prime top_base."""
    D = sum(d*d for d in degrees)
    M = isqrt(D)
    assert M*M == D and M % max(degrees) == 0
    a = M//max(degrees)
    check = a
    while check % top_base == 0:
        check //= top_base
    assert a > 1 and check == 1
    pool = set(M//d for d in degrees if d > 1 and M//d % top_base == 0)
    if reduce(lcm, pool, 1) != M:
        return 'lcm'
    # c_Y >= (a/gcd(a,b))^2 from lcm integrality.  One label of each
    # chosen degree is sufficient for this necessary minimum-cost test.
    bounds = {1:0}
    for b in pool:
        cost = (a//gcd(a,b))**2
        for old, spent in tuple(bounds.items()):
            if spent+cost >= a*a:
                continue
            new = lcm(old,b)
            bounds[new] = min(bounds.get(new,a*a),spent+cost)
    return 'passes' if M in bounds else 'cost'


def dyadic_trace_possible(degrees, *, use_modular=True):
    """Test the dyadic principal row at a square dimension and power-of-two top denominator."""
    from functools import lru_cache
    total = sum(y*y for y in degrees)
    M = isqrt(total)
    n = max(degrees)
    a = M//n
    assert M*M == total and M % n == 0 and a > 1 and a & (a-1) == 0
    principal = 2*degrees.count(n) > a*a
    budget = a*a-1
    terms = sorted(Counter(M//y for y in degrees if y>1 and (M//y)%2 == 0).items(),
                   reverse=True)
    groups = []
    scale = 256
    for b,m in terms:
        d = M//b
        step = a//gcd(a,b)
        parity = None
        if n%2 == 0 and (n//gcd(n,d))%2:
            parity = 0
        if principal and step%2:
            parity = d%2
        # For class trace t=sum trace(alpha), the smallest relaxed cost
        # is obtained by balanced coordinates because norm_floor is convex.
        def norm_floor(s):
            adjust = parity is not None and (s-parity)%2 != 0
            return s*s + int(adjust)
        bounds = []
        max_trace = isqrt(budget*m)//step
        for t in range(-max_trace,max_trace+1):
            if parity is not None and (t-m*parity)%2:
                continue
            j,r = divmod(abs(t),m)
            least = (m-r)*norm_floor(j)+r*norm_floor(j+1)
            for B in range(least,budget//(step*step)+1):
                if parity is not None and (B-m*parity)%2:
                    continue
                # An upper, rounded upward, for scale * d sum(mean |chi|).
                unit = scale*d*step
                cap = isqrt(unit*unit*m*B)+1
                if parity == 0:
                    cap = min(cap, isqrt((unit*B)**2//2)+1)
                else:
                    cap = min(cap, unit*B)
                bounds.append((step*d*t,step*step*B,cap))
        groups.append((step*d,bounds))
    suffix = [0]*(len(groups)+1)
    for i in reversed(range(len(groups))):
        suffix[i] = gcd(suffix[i+1],groups[i][0])
    @lru_cache(None)
    def best(i, target, remaining):
        if i == len(groups):
            return 0 if target==0 and remaining==0 else -1
        if target % suffix[i]:
            return -1
        maximum = -1
        for signed,cost,cap in groups[i][1]:
            if cost <= remaining:
                child = best(i+1,target-signed,remaining-cost)
                if child >= 0:
                    maximum = max(maximum,child+cap)
        return maximum
    bound = best(0,-1,budget)
    return bound >= (scale*(M-1) if use_modular else 0)


if __name__ == '__main__':
    import json,sys
    from pathlib import Path
    for path in map(Path,sys.argv[1:]):
        counts = Counter()
        for line in path.read_text().splitlines():
            degrees = json.loads(line)
            if row_cover(degrees) != 'passes':
                counts['cover'] += 1
            elif not dyadic_trace_possible(degrees,use_modular=False):
                counts['trace'] += 1
            elif not dyadic_trace_possible(degrees):
                counts['modular'] += 1
            else:
                counts['passes'] += 1
        print(path.name,dict(counts))
