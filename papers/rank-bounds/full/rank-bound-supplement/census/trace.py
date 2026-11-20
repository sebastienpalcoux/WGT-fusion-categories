"""Necessary short-vector conditions on Galois orbits of character rows."""
from collections import Counter
from functools import lru_cache
from itertools import groupby
from math import gcd, isqrt, lcm, prod
from typing import NamedTuple

from engine.primes import primes_up_to


class CyclotomicData(NamedTuple):
    degree: int
    trace_denominator: int
    prime_power_base: int


@lru_cache(maxsize=128)
def cyclotomic_data(n):
    """Degree and common denominator of normalized traces on Z[zeta_n]."""
    primes = [p for p in primes_up_to(n) if n % p == 0]
    denominator = prod(p-1 for p in primes)
    return CyclotomicData((n//prod(primes))*denominator, denominator,
                          primes[0] if len(primes) == 1 else 0)


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
        field = cyclotomic_data(f)
        p = field.prime_power_base
        eligible = [y for y in mult if y > 1 and gcd(x, y) > 1
                    and (not p or D//(y*y) % (p*p) == 0)
                    and (not p or f != p*p or x//gcd(x, p) % y == 0)]
        if gcd(*eligible) != 1:
            return f

        sums = 1  # Bit n means n labels can be partitioned into feasible orbits.

        for h in range(1, min(m, f-1)+1):
            # Conductor divides formal codegree, so K_X lies in Q(zeta_f).
            if field.degree % h:
                continue
            energy = h*(f-h)
            coords = []
            minimum_energy = 0
            for y in eligible:
                copies = mult[y]
                # Q(zeta_f) intersection Q(zeta_g) is Q(zeta_gcd(f,g)).
                # The normalized trace of an integral entry in Q(zeta_g)
                # has denominator dividing phi(rad g), its entry-field degree,
                # and hence the lcm of the possible orbit lengths at degree y.
                joint = cyclotomic_data(gcd(f, D//(y*y))).trace_denominator
                step = (y//gcd(x, y))*(h//gcd(h, joint, orbit_length_multiple(copies)))
                if p and x % p == 0:
                    step = lcm(step, p*y//gcd(x, p*y))
                if p and (h*y) % p:
                    # Every row in the orbit is p-principal. Its summed entry
                    # must equal h*y modulo p, on every copy of this degree.
                    if step % p == 0:
                        minimum_energy = energy + 1
                        break
                    residue = (h * y * pow(step, -1, p)) % p
                    least = min(residue, p-residue)*step
                    minimum_energy += copies*least*least
                    if minimum_energy > energy:
                        break
                if step*step > energy:
                    continue
                if y == x:
                    if h*step*step <= energy:
                        coords.append((h*y*step, h*step*step, h*y*y))
                    copies -= h
                coords.extend([(y*step, step*step, y*y)] * copies)
            if minimum_energy > energy:
                continue
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
