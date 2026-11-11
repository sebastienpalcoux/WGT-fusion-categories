"""Necessary short-vector conditions on Galois orbits of character rows."""
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm

@lru_cache(maxsize=None)
def prime_power_base(f):
    """The prime base of a prime power, or None."""
    if f < 2:
        return None
    p = next((d for d in range(2, isqrt(f)+1) if f % d == 0), f)
    n = f
    while n % p == 0:
        n //= p
    return p if n == 1 else None


def trace_obstruction(degrees, cutoff=36):
    """Return an impossible codegree, or None. Skipping f > cutoff is safe."""
    mult = Counter(degrees)
    D = sum(x*x for x in degrees)
    for x, m in sorted(mult.items(), reverse=True):
        f = D // (x*x)
        if f > cutoff or x == 1:
            break
        p = prime_power_base(f)
        eligible = [y for y in mult if y > 1 and gcd(x, y) > 1
                    and (not p or D // (y*y) % (p*p) == 0)
                    and (not p or f != p*p or x // gcd(x, p) % y == 0)]
        if gcd(*eligible) != 1:
            return f

        def feasible(h):
            coordinates = []
            energy = h*(f-h)
            for y in eligible:
                copies = mult[y]
                step = y // gcd(x, y) * (h // gcd(h, lcm(*range(1, copies+1))))
                if p and x % p == 0:
                    step = lcm(step, p*y // gcd(x, p*y))
                if y == x:
                    coordinates.append((y, step, h))
                    copies -= h
                coordinates.extend([(y, step, 1)] * copies)
            coordinates = [(y,s,k) for y,s,k in coordinates if k*s*s <= energy]
            coordinates.sort(key=lambda a: (a[2]*a[1]**2, a[0]), reverse=True)
            length = len(coordinates)
            divisors, norms = [0]*(length+1), [0]*(length+1)
            for i in reversed(range(length)):
                y, step, k = coordinates[i]
                divisors[i] = gcd(divisors[i+1], k*y*step)
                norms[i] = norms[i+1] + k*y*y

            @lru_cache(maxsize=25000)
            def solve(i, target, left, needed):
                if i == length:
                    return target == left == needed == 0
                if target % divisors[i] or max(target*target, needed*needed) > left*norms[i]:
                    return False
                y, step, k = coordinates[i]
                cap = isqrt(left // (k*step*step))
                for multiple in range(-cap, cap+1):
                    z = step * multiple
                    if solve(i+1, target-k*y*z, left-k*z*z, max(0, needed-k*y*abs(z))):
                        return True
                return False
            return solve(0, -h, energy, isqrt(D-1) if h == 1 else 0)

        sizes = [h for h in range(1, min(m, f-1)+1) if feasible(h)]
        sums = {0}
        for n in range(1, m+1):
            if any(n-h in sums for h in sizes):
                sums.add(n)
        if m not in sums:
            return f
    return None
