"""Direct odd-dimensional rank bound from paired ordinary formal codegrees."""
from functools import lru_cache
from math import gcd, isqrt, lcm

@lru_cache(100000)
def factor(n):
    factors = []
    p = 3  # every denominator and final dimension in this script is odd
    while p*p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            factors.append((p, exponent))
        p += 2
    if n > 1:
        factors.append((n, 1))
    return factors

def divisors(q):
    result = [1]
    for prime, exponent in factor(q):
        old = result.copy()
        power = 1
        for _ in range(exponent):
            power *= prime
            result.extend(d * power for d in old)
    return result

def covered(dimension):
    if dimension < 6087156075:
        return True
    exponents = [a for p, a in factor(dimension) if a >= 2]
    if len(exponents) <= 2:
        return True
    sizes = sorted(exponents)
    if len(sizes) == 3:
        a, b, c = exponents
        return (a <= 5 or b <= 4 or c <= 3 or
                all(e <= cap for e, cap in zip(sizes, (4, 5, 7))))
    if len(sizes) == 4:
        d, c, b, a = sizes
        return b <= 3 or (d == 2 and c <= 3 and a//2 + b//2 <= 5)
    if len(sizes) == 5:
        return sizes[-1] <= 3 or (sizes[-1] <= 4 and sizes[-2] <= 3 and sizes[0] == 2)
    return False

def verify(rank, least=3):
    tuples = 0
    maximal_dimension = 0
    survivors = []
    def search(pairs, p, q, previous, L):
        nonlocal tuples, maximal_dimension
        if pairs == 1:
            for k in divisors(q):
                if (2*k + 1) % p:
                    continue
                D = q * ((2*k + 1) // p)
                if D % L or D//k < previous:
                    continue
                tuples += 1
                maximal_dimension = max(D, maximal_dimension)
                if not covered(D):
                    survivors.append(D)
            return
        first = max(previous, 2*q//p + 1) | 1
        last = (2*pairs + 1)*q//p
        for x in range(first, last + 1, 2):
            num, den = p*x - 2*q, q*x
            g = gcd(num, den)
            search(pairs - 1, num//g, den//g, x, lcm(L, x))
    search((rank - 1)//2, 1, 1, least, 1)
    assert not survivors, (rank, survivors)
    print(rank, least, tuples, maximal_dimension, len(survivors))

def square_prefix(pairs, p=1, q=1, previous=2, L=1):
    """Solve 1 = 42 sum(1/u^2) + 1/n^2, with u | n."""
    if not pairs:
        n = isqrt(q)
        if p == 1 and n*n == q and n >= previous and n % L == 0:
            yield n
        return
    low = max(previous, isqrt(42*q//p) + 1)
    high = isqrt((42*pairs + 1)*q//p)
    for u in range(low, high + 1):
        num, den = p*u*u - 42*q, q*u*u
        g = gcd(num, den)
        yield from square_prefix(pairs-1, num//g, den//g, u, lcm(L, u))

def modular_tail(weight, pairs, p, q, previous=1, L=1):
    """Return odd n with p/q = weight*sum(1/u^2) + 1/n^2 and all u | n."""
    if pairs == 0:
        n = isqrt(q)
        if p == 1 and n*n == q and n % 2 and n >= previous and n % L == 0:
            yield n
        return
    first = max(previous, isqrt(weight*q//p) + 1) | 1
    last = isqrt((weight*pairs + 1)*q//p)
    for u in range(first, last + 1, 2):
        num, den = p*u*u - weight*q, q*u*u
        g = gcd(num, den)
        yield from modular_tail(weight, pairs-1, num//g, den//g, u, lcm(L, u))

def free_group_case(order, pairs):
    dimensions = []
    # t = 2|H| sum(1/u^2) + 1/n^2; t is odd and square-free.
    for t in range(1, 2*order*pairs + 2, 2):
        if gcd(t, order) != 1 or any(a > 1 for _, a in factor(t)):
            continue
        for n in modular_tail(2*order, pairs, t, 1):
            if gcd(n, order) == 1:
                dimensions.append(order*t*n*n)
    print("Rank 15, group", order, "dimensions", sorted(dimensions))
    assert all(covered(D) for D in dimensions)

if __name__ == '__main__':
    print('rank minimum tuples maximum_dimension uncovered')
    for r in (3, 5, 7, 9, 11):
        verify(r)
    verify(13, 5)
    solutions = list(square_prefix(3))
    print('Rank 13, paired prefix 3,7,43, square equation:', solutions)
    assert not solutions
    verify(15, 7)
    for order, pairs in [(21, 5), (39, 4), (55, 4)]:
        free_group_case(order, pairs)
    single_stage = [55*n*n for n in modular_tail(10, 4, 1, 1) if n % 5]
    print('Rank 15, modular first quotient of dimension 11 times a square:', single_stage)
    assert all(covered(D) for D in single_stage)
