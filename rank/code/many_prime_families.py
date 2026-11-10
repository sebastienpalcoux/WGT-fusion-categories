"""Divisor-closed four- and five-repeated-prime rules from BDP, §§9–15.

The argument consists of (distinct prime, exponent) pairs. Exponent-one
decorations are immaterial. None means this predicate has no certificate.
The optional covered callback accepts an integer of the form p^H (abcd)^2
and must certify every divisor including arbitrary coprime squarefree
decoration. It is only called on strict predecessors.
"""
from itertools import combinations
from math import isqrt, prod

# BDP §14.2, output of both prime_projection.py and independent_projection.py.
S3 = frozenset((3, 5, 7, 11, 13, 17, 19, 31, 37))
S4 = frozenset((3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 109, 113, 127,
    131, 151, 163, 181, 191, 193, 197, 211, 251, 283))


def residue_survives(p, q, u, v, h, k):
    """Necessary simultaneous vacuum, size, and endpoint conditions."""
    x, y = p ** h, q ** k
    modulus = u * v
    residue = (x * y) % modulus
    if residue in (1, modulus - 1):  # BDP 9.6; h,k >= 2.
        return False
    if u * (p - 1) > p * min(x, y) or min(q, v) * (p - 1) > p * u:
        return False  # BDP 9.7.
    odd = p != 2
    if odd and v * (p - 1) > p * min(u, max(x // p, y // q)):
        return False  # BDP 10.12.
    if odd and k == 2:
        if p == 3 and 4 * u > 9 * (x // p) + 2:
            return False
        if p >= 5 and 2 * (p - 1)**2 * u > (2 * (p-1)**2 + p*p) * (x // p):
            return False  # BDP 10.16.
    j0 = pow(residue * v, -1, u)
    j1 = pow(residue * u, -1, v)
    for sign in (1, -1):
        ju, jv = sign * j0 % u, sign * j1 % v
        if ju * v + jv * u >= modulus:
            continue  # BDP 6.3 (vacuum budget); C>1 and equality impossible.
        if any((t-1) % j and (t+1) % j for t, j in ((u, ju), (v, jv))):
            continue  # BDP 6.5.
        if odd:
            if (u-1) % ju or (v-1) % jv:
                continue
            hu, hv = (u-1) // ju, (v-1) // jv
            if min(hu, hv) < 5 or not (hu % 2 and hv % 2):
                continue
            if 3*ju > u+1 or 3*jv > v+1:
                continue
            if ju+3 > hv and jv+3 > hu:
                continue  # BDP §14.1.
        bu, bv = ju * (u+1-ju), jv * (v+1-jv)
        if v*(p-1) > p*min(x, y)*bu or q*(p-1) > p*min(v*bu, u*bv):
            continue  # BDP 9.7, remaining two prime-size bounds.
        upper = max(2, (u+1) // ju - 1)
        if q > v * upper + 1 or v > min(x, y) * upper + 1:
            continue  # BDP 9.14 and interchange of the two higher primes.
        return True
    return False


def odd_residue_periods(p, q, u, v):
    """BDP §14.1 exact allowed k modulo d in G/<p> for odd positions.

    This optional calculation is independent of both exponent caps.
    Arithmetic modulo u*v realizes the two finite multiplicative groups.
    """
    def endpoint_indices(t):
        divisors = set()
        for a in range(1, isqrt(t-1) + 1):
            if (t-1) % a == 0:
                divisors.update((a, (t-1)//a))
        return [(d, (t-1)//d) for d in divisors if d >= 5 and d % 2
                and 3*((t-1)//d) <= t+1]

    modulus = u*v
    target = set()
    for hu, ju in endpoint_indices(u):
        for hv, jv in endpoint_indices(v):
            if ju+3 > hv and jv+3 > hu:
                continue
            # x*v=hu (mod u), x*u=hv (mod v) for the positive anomaly.
            a = hu * pow(v, -1, u) % u
            b = hv * pow(u, -1, v) % v
            x = a + u*((b-a)*pow(u, -1, v) % v)
            target.update((x, -x % modulus))
    subgroup = {1}
    x = p % modulus
    while x != 1:
        subgroup.add(x)
        x = x*p % modulus
    exponents = []
    inverse = pow(q, -1, modulus)
    d, inverse_coset = 0, 1
    while True:
        if any(x*inverse_coset % modulus in subgroup for x in target):
            exponents.append(d)
        inverse_coset = inverse_coset*inverse % modulus
        d += 1
        if inverse_coset in subgroup:
            return d, tuple(exponents)


def dyadic_character_box(q, u, v):
    """Two unbounded higher exponents, via the character of either singleton."""
    for t, w in ((u, v), (v, u)):
        if t % 8 != 1:
            continue
        if pow(q, (t-1)//2, t) != 1 or pow(w, (t-1)//2, t) != t-1:
            continue
        remaining = (t+1)//2
        divisor = 3
        while divisor*divisor <= remaining:
            if remaining % divisor == 0:
                if divisor % 4 != 1:
                    break
                while remaining % divisor == 0:
                    remaining //= divisor
            divisor += 2
        else:
            if remaining % 4 == 1:
                return True
    return False


def base_square_box(p, a, q, b, u, v):
    """p<q, u<v. Test all divisor-minimal roots of the two-square box."""
    if u < p or u < 5 or a < 6 or b < 4:
        return True  # BDP 9.3, 10.17 and 13.1.
    if p == 2 and dyadic_character_box(q, u, v):
        return True
    for h in range(3, a // 2 + 1):
        if p == 2 and ((h == 3 and q not in S3) or (h == 4 and q not in S4)):
            continue
        if p == 3 and h <= 4:
            continue
        for k in range(2, b // 2 + 1):
            if h + k < 6 or (p == 2 and h+k == 6 and h != k):
                continue
            if h == 4 and k == 2 and p in (3, 5, 7, 11, 13, 17, 19):
                continue
            if residue_survives(p, q, u, v, h, k):
                return False
    return True


def many_prime_family(repeated, covered=None):
    """Return a sufficient WGT family or None. Prime inputs are assumed."""
    items = sorted((p, e) for p, e in repeated if e >= 2)
    if len(items) == 4:
        if sum(e > 3 for _, e in items) <= 1:
            return "BDP 11.3: one unrestricted exponent, three cubes"
        for low in combinations(items, 2):
            (u, eu), (v, ev) = low
            if ev > 3 or eu > 3:
                continue
            (p, a), (q, b) = [entry for entry in items if entry not in low]
            cube_lift = min(eu, ev) == 2 or any(
                t % 2 and (t-1) % (p*q) for t in (u, v))
            if cube_lift and base_square_box(p, a, q, b, u, v):
                return "Two-square divisor descent and BDP 11.3 cube lift"
    if len(items) == 5:
        if max(e for _, e in items) <= 3:
            return "BDP 15.5: five cubes"
        for p, h in items:
            companions = [(t, e) for t, e in items if t != p]
            if any(e > 3 for _, e in companions):
                continue
            minimum = companions[0][0]
            lift = any(e == 2 for _, e in companions) or (
                minimum < p and any((t-1) % p for t, _ in companions[1:]))
            if not lift:
                continue
            if h <= 4:
                return "BDP 15.1, 15.6: quartic base and five-prime cube lift"
            if covered and any(e == 3 for _, e in companions):
                base = p**h * prod(t*t for t, _ in companions)
                if covered(base):
                    return "BDP 15.6: conditional five-prime cube lift"
    return None
