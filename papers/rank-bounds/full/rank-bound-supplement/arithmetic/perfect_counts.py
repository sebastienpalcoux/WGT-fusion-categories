"""Finite arithmetic in the classification of perfect categories through rank 44."""
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import isqrt


def prime_powers(limit, modulus):
    """Prime powers p^d <= limit congruent to 1 mod modulus, with d >= 3."""
    result = []
    for p in range(2, isqrt(limit) + 1):
        if p**3 > limit:
            break
        if any(p % d == 0 for d in range(2, isqrt(p) + 1)):
            continue
        value, degree = p**3, 3
        while value <= limit:
            if value % modulus == 1:
                result.append((p, degree))
            value, degree = value*p, degree+1
    return result


def prime_power(n):
    if n < 2:
        return False
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return True


def powers_of(base, bound):
    values, value = [], base
    while value <= bound:
        values.append(value)
        value *= base
    return values


# Elementary modules when a prime-order element acts freely off zero.
assert prime_powers(2017, 7) == [(2,3),(2,6),(2,9),(3,6),(11,3)]
assert prime_powers(4321, 5) == [(2,4),(2,8),(2,12),(3,4),(7,4),(11,3)]
assert prime_powers(1321, 5) == [(2,4),(2,8),(3,4)]

# Orbit size, w(L), involution fixed points, five-element fixed points.
# Every nonfree orbit costs at least four, including its contribution to t.
kinds = [(30,3,2,0),(20,8,0,0),(15,3,3,0),(12,24,0,2),
         (10,7,2,0),(6,15,2,1),(5,9,1,0),(1,17,1,1)]
patterns = {}
for count in range((44-19)//4 + 1):
    for pattern in combinations_with_replacement(kinds, count):
        n2 = 1 + sum(orbit[2] for orbit in pattern)
        if not any(orbit[3] for orbit in pattern):
            continue
        free = 44 - 19 - count - sum(orbit[1] for orbit in pattern) - 3*(n2 % 2)
        if free >= 0:
            patterns[tuple(sorted(orbit[0] for orbit in pattern))] = free
assert patterns == {(1,):7,(1,30):3,(1,15):0,(6,):6,(6,15):5,
                    (6,30):2,(6,15,30):1}

# Characteristic elementary subgroups are unions of neutral orbits.
singleton_sizes = {1+a+60*b for a in range(2) for b in range(8)}
assert sorted(n for n in singleton_sizes if prime_power(n)) == [2,61,121,181,241,361,421]
size30 = {1+a+30*b+60*c for a in range(2) for b in range(2) for c in range(4)}
assert sorted(n for n in size30 if n > 1 and n & (n-1) == 0) == [2,32]
assert sorted(1+a+15*b for a,b in product(range(2),repeat=2)
              if prime_power(1+a+15*b)) == [2,16,17]

# Fixed-point-free five-element, 256-element binary module: n2=16 and n3=n5=1.
triples = [(a,b,c) for a,b,c in product(range(5),range(9),range(18))
           if 60*a+30*b+15*c == 255 and 2*b+3*c == 15 and a+2*b+c <= 10]
assert triples == [(3,0,5)]
assert [v for v in powers_of(3,5761) if v % 5 == 1] == [81]
assert [v for v in powers_of(2,361) if v % 5 == 1] == [16,256]

# Possible prefixes of a nondecreasing three-term reciprocal-square sum of 1/9.
prefixes = []
for a in range(4,6):
    remaining = Fraction(1,9) - Fraction(1,a*a)
    b = a
    while Fraction(2,b*b) >= remaining:
        tail = remaining - Fraction(1,b*b)
        if tail > 0:
            prefixes.append((a,b))
            assert tail.numerator != 1 or isqrt(tail.denominator)**2 != tail.denominator
        b += 1
assert prefixes == [(4,5),(4,6),(5,5)]

# Second Chern class from the weights, then outer orbits on the three
# primary factors of H^4(A6,Z). The independent sign actions are proved
# in the manuscript. The A5 count applies to its doubles and deleted gaugings.
def second_chern(weights, modulus):
    return sum(w*v for i,w in enumerate(weights) for v in weights[i+1:]) % modulus


assert second_chern([0,1,2,3,0,2],4) == 3
assert second_chern([0,1,2,0,0,0],3) == 2
assert second_chern([0,1,2,0,1,2],3) == 1
parameters = product(range(4),range(3),range(5))
a6_orbits = {min((a,(s*b)%3,(t*c)%5) for s,t in product((-1,1),repeat=2))
             for a,b,c in parameters}
assert len(a6_orbits) == 24
assert len({(a,b,min(c,-c%5)) for a,b,c in product(range(2),range(3),range(5))}) == 18
assert len({min(t,49*t%120) for t in range(1,120,2)}) == 36

# Degrees in the A5 gauging of the deleted module, by conjugacy class.
electric = [1,3,3,4,5]
neutral = electric + [5]*3 + [15] + [10]*2 + [20]
degrees = [neutral, [30]*8+[60]*2, [40]*12, [48]*5, [48]*5]
assert [len(sector) for sector in degrees] == [12,10,12,5,5]
assert sum(map(len,degrees)) == 44
assert sum(d*d for sector in degrees for d in sector) == 57600
assert sum(sector.count(1) for sector in degrees) == 1

print('Verified rank-44 prime powers, seven orbit patterns, elementary subgroup sizes, '
      'reciprocal-square prefixes, cohomology orbit counts and deleted-gauging degrees.')
