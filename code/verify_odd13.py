"""Check rank-13 reciprocal tuples against the proved center exponent caps."""
from fractions import Fraction
from math import lcm


def tuples(pairs, previous, remainder=Fraction(1), prior=()):
    if not pairs:
        if remainder.numerator == 1:
            d = remainder.denominator
            if d >= previous and d % 2 and d % lcm(*prior) == 0:
                yield prior, d
        return
    low = max(previous, 1 + 2 // remainder)
    high = (2 * pairs + 1) // remainder
    for x in range(low | 1, high + 1, 2):
        yield from tuples(pairs - 1, x, remainder-Fraction(2, x), (*prior, x))


def exponents(n):
    e = []
    p = 3
    while p*p <= n:
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        if a:
            e.append(2*a)
        p += 2
    if n > 1:
        e.append(2)
    return e


def covered(n):
    ex = exponents(n)
    if len(ex) <= 2:
        return True
    small = sorted(ex)
    if len(ex) == 3:
        a,b,c = ex
        return a <= 5 or b <= 4 or c <= 3 or all(x <= y for x,y in zip(small,(4,5,7)))
    if len(ex) == 4:
        d,c,b,a = small
        return b<=3 or (d==2 and c<=3 and a//2+b//2<=5)
    if len(ex) == 5:
        return small[-1]<=3 or (small[-1]==4 and small[-2]<=3 and small[0]==2)
    return False


cases = list(tuples(6, 7))
above = [(xs, n) for xs, n in cases if n >= 121275 and xs.count(7) <= 2]
open_cases = [(xs,n) for xs,n in above if not covered(n)]
print(len(cases),len(above),open_cases)
assert (len(cases),len(above)) == (520,17)
assert sorted(open_cases) == sorted([((7,7,11,11,31,4785),1038345),((7,7,9,11,95,575),7571025)])
