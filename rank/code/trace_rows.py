"""Rational orbit-sum tests on a perfect simple integral character table."""
from families import factors
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm


def nonunit_trace_row(degrees, x, specified_h=None, ordinary=True, principal=False):
    """Test the orbit sum at one nonunit degree x by exact short-vector search."""
    counts = Counter(degrees)
    total = sum(d*d for d in degrees)
    f = total // (x*x)

    factorization = factors(f) if principal else ()
    p = factorization[0][0] if len(factorization) == 1 else None

    def for_orbit(h):
        # Each tuple records (degree, trace step, number of equal entries).
        # Galois conjugates of X have equal trace coordinates in its orbit sum.
        if total % 2 and h % 2:
            return False
        pair_size = 2 if total % 2 else 1
        coordinates = []
        for y, m in counts.items():
            if y == 1 or gcd(x, y) == 1:
                continue
            if p:
                fy = total // (y*y)
                if fy % (p*p):
                    continue  # A nonzero entry belongs to the principal p-block.
                if f == p*p and (x // gcd(x, p)) % y:
                    continue  # BDP Proposition 5.7 for nonzero rational traces.
            step = y // gcd(x, y)
            # The crossing field is contained in both character fields. A
            # y-field has degree at most its degree multiplicity.
            overlap = gcd(h, lcm(*range(1, m + 1)))
            step *= h // overlap
            if p and x % p == 0:
                # The integer (x/y) times the trace is zero modulo p.
                step = lcm(step, p*y // gcd(x, p*y))
            if y == x:
                coordinates.append((y, step, h))
                m -= h
            if m % pair_size:
                return False
            coordinates.extend([(y, step, pair_size)] * (m // pair_size))
        budget = h * (f-h)
        # Put expensive entries first; zero is their only admissible value
        # once their cost exceeds the remaining squared norm.
        coordinates.sort(key=lambda a: (a[1]*a[1]*a[2], a[0]), reverse=True)
        tail_gcd = [0]*(len(coordinates)+1)
        tail_norm = [0]*(len(coordinates)+1)
        for i in reversed(range(len(coordinates))):
            y,s,t = coordinates[i]
            tail_gcd[i] = gcd(y*s*t, tail_gcd[i+1])
            tail_norm[i] = t*y*y + tail_norm[i+1]

        @lru_cache(maxsize=25000)
        def solve(i, target, energy, needed):
            if i == len(coordinates):
                return target == 0 and energy == 0 and needed <= 0
            if target % tail_gcd[i] or max(target*target, needed*needed) > energy*tail_norm[i]:
                return False
            y,s,t = coordinates[i]
            maxval = isqrt(energy//(t*s*s))
            for val in range(-maxval,maxval+1):
                d = s*val
                if solve(i+1,target-y*t*d,energy-t*d*d,max(0,needed-y*t*abs(d))):
                    return True
            return False

        root = isqrt(total)
        ceiling = root if root*root == total else root+1
        needed = ceiling - 1 if h == 1 and ordinary else 0
        return solve(0, -h, budget, needed)

    if specified_h is not None:
        return for_orbit(specified_h)
    return any(for_orbit(h) for h in range(1,min(f-1,counts[x])+1))


def max_trace_possible(degrees):
    return nonunit_trace_row(degrees, max(degrees))
