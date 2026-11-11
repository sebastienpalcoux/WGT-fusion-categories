"""Necessary short-vector conditions on Galois orbits of character rows."""
from collections import Counter
from functools import lru_cache
from itertools import groupby
from math import gcd, isqrt, lcm


@lru_cache(maxsize=64)
def prime_power_base(f):
    """The prime base of a prime power, or None."""
    if f < 2:
        return None
    p = next((d for d in range(2, isqrt(f)+1) if f % d == 0), f)
    n = f
    while n % p == 0:
        n //= p
    return p if n == 1 else None


@lru_cache(maxsize=32)
def orbit_length_multiple(copies):
    return lcm(*range(1, copies+1))


@lru_cache(maxsize=256)
def equal_group(count, budget):
    """For a given sum and sum of squares retain the largest sum of absolute values."""
    options = {(0, 0): 0}
    for _ in range(min(count, budget)):
        extended = options.copy()  # one more coordinate may be zero
        for (total, squares), absolute in options.items():
            for a in range(1, isqrt(budget-squares)+1):
                q = squares+a*a
                value = absolute+a
                for w in (total-a, total+a):
                    if value > extended.get((w, q), -1):
                        extended[w, q] = value
        options = extended
    by_sum = {}
    for (w, q), absolute in options.items():
        by_sum.setdefault(w, []).append((q, absolute))
    return {w: sorted(choices) for w, choices in by_sum.items()}


@lru_cache(maxsize=256)
def group_energies(count, budget):
    return tuple({q for choices in equal_group(count, budget).values() for q, _ in choices})


@lru_cache(maxsize=4096)
def orbit_feasible(coordinates, h, energy, needed):
    """An orbit contributes squared norm h*(f-h) away from the unit."""
    groups = []
    for (moment, cost, norm), copies in groupby(coordinates):
        count = sum(1 for _ in copies)
        if count > 1 and energy//cost <= 48:
            options = equal_group(count, energy//cost)
            groups.append((moment, cost, count*norm, count, options))
        else:
            groups.extend([(moment, cost, norm, 1, None)]*count)
    length = len(groups)
    divisors = [0]*(length+1)
    norms = [0]*(length+1)
    energy_bits = [1]*(length+1)
    moduli = [1]*length
    inverses = [0]*length
    mask = (1 << (energy+1))-1
    for i in reversed(range(length)):
        moment, cost, norm, count, options = groups[i]
        divisors[i] = gcd(divisors[i+1], moment)
        norms[i] = norms[i+1]+norm
        if divisors[i+1]:
            modulus = divisors[i+1]//divisors[i]
            moduli[i] = modulus
            if modulus > 1:
                inverses[i] = pow(moment//divisors[i], -1, modulus)
        possible = energy_bits[i+1]
        if options:
            costs = group_energies(count, energy//cost)
        else:
            costs = (j*j for j in range(1, isqrt(energy//cost)+1))
        for q in costs:
            possible |= energy_bits[i+1] << (cost*q)
        energy_bits[i] = possible & mask

    @lru_cache(maxsize=8192)
    def solve(i, target, left, need):
        if i == length:
            return target == left == need == 0
        if not (energy_bits[i] >> left) & 1 or target % divisors[i]:
            return False
        bound = left*norms[i]
        if target*target > bound or need*need > bound:
            return False
        moment, cost, _, count, options = groups[i]
        if i+1 == length:
            j = target//moment
            if options:
                return any(cost*q == left and moment*absolute >= need
                           for q, absolute in options.get(j, ()))
            return cost*j*j == left and abs(target) >= need
        cap = isqrt(count*(left//cost))
        modulus = moduli[i]
        # Only these sums of coordinates leave a compatible scalar product.
        first = -cap
        if modulus > 1:
            first += (target//divisors[i]*inverses[i]+cap) % modulus
        for j in range(first, cap+1, modulus):
            contribution = moment*j
            if options:
                for q, absolute in options.get(j, ()):
                    if cost*q > left:
                        break
                    remaining = max(0, need-moment*absolute) if need else 0
                    if solve(i+1, target-contribution, left-cost*q, remaining):
                        return True
            else:
                remaining = max(0, need-abs(contribution)) if need else 0
                if solve(i+1, target-contribution, left-cost*j*j, remaining):
                    return True
        return False

    answer = solve(0, -h, energy, needed)
    solve.cache_clear()
    return answer


def trace_obstruction(degrees, cutoff=36, *, dimension=None, multiplicities=None):
    """Return an impossible codegree, or None. Skipping f > cutoff is safe."""
    mult = multiplicities if multiplicities is not None else Counter(degrees)
    D = dimension if dimension is not None else sum(x*x*n for x,n in mult.items())
    for x, m in sorted(mult.items(), reverse=True):
        f = D//(x*x)
        if f > cutoff or x == 1:
            break
        p = prime_power_base(f)
        eligible = [y for y in mult if y > 1 and gcd(x, y) > 1
                    and (not p or D//(y*y) % (p*p) == 0)
                    and (not p or f != p*p or x//gcd(x, p) % y == 0)]
        if gcd(*eligible) != 1:
            return f

        sums = 1  # Bit n means n labels can be partitioned into feasible orbits.
        for h in range(1, min(m, f-1)+1):
            energy = h*(f-h)
            coords = []
            for y in eligible:
                copies = mult[y]
                step = y//gcd(x, y)*(h//gcd(h, orbit_length_multiple(copies)))
                if p and x % p == 0:
                    step = lcm(step, p*y//gcd(x, p*y))
                if step*step > energy:
                    continue
                if y == x:
                    if h*step*step <= energy:
                        coords.append((h*y*step, h*step*step, h*y*y))
                    copies -= h
                coords.extend([(y*step, step*step, y*y)] * copies)
            coords.sort(key=lambda c: (c[1], c[0]), reverse=True)
            if orbit_feasible(tuple(coords), h, energy, isqrt(D-1) if h == 1 else 0):
                # Allow this orbit size to be used any number of times.
                shift = h
                while shift <= m:
                    sums |= sums << shift
                    shift *= 2
                if sums & (1 << m):
                    break
        else:
            return f
    return None
