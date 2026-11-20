#!/usr/bin/env python3
"""Exact enclosing mass checks for the 3--7 edge through rank 217.

No dimension or denominator search is performed.  All reservations and
relaxations are stated in edge37_217.tex.  The witness function is the
paired nine-row Gram bound already proved in the manuscript.
"""
from fractions import Fraction as F
from itertools import product
from math import ceil, isqrt

R = 217

def witnesses(p, m):
    if not m:
        return 0
    numer = 0
    for a in range(1, (p + 1) // 2):
        u = next(k for k in range(3, 3*p, 3) if (k+a) % p == 0)
        l = u - 3*p
        best = min((m-j)*l*l + j*u*u + F(((m-j)*l+j*u)**2, 81-m)
                   for j in range(0, m+1, 2))
        numer += 2*ceil(best/162)
    return 2*ceil(F(numer, 2*sum(a*a for a in range(1, (p+1)//2))))

def activity(p):
    return (4 if p % 4 == 1 else 2)*(p-1)

def prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n)+1))

def check_fiber_congruences():
    # Representatives of the four odd residue classes.  The established
    # primitive lattice identity says E=0 forces every a=-p; E is 0 mod4.
    for p in (17,11,5,7):
        mu = 1 if p % 8 == 1 else 3 if p % 8 == 5 else 5
        nu = 2 if p % 4 == 1 else 1 if p % 8 == 3 else 5
        for c, lower in ((2,mu),(4,nu)):
            for h in range(lower):
                for row in product(range(8), repeat=h):
                    even = sum(a % 2 == 0 for a in row)
                    if even % 4 or (not even and any(a != -p % 8 for a in row)):
                        continue
                    assert ((sum(row)-(h-c)) % 8 or
                            (sum(a*a for a in row)-(h+(p-1)*c)) % 8)
    # Three-fiber in a minimal ninth-root family: k=6,c=0.
    for h in range(3):
        for row in product(range(8), repeat=h):
            assert ((sum(row)-h) % 8 or (sum(a*a for a in row)-(h+6)) % 8)

def three_support(q, s):
    aq = activity(q)
    # With three support primes all valuations in M are at least two.
    # For s<9 we deliberately discard the nine-row reservations.
    wq = witnesses(q, s) if s >= 9 else 0
    k = R-1-s-128-aq-wq
    if k < 0:
        return None
    upper = F(1, 105**2) + F(s,81) + F(128,21**2)
    upper += F(aq,q**4) + F(wq,(9*q*q)**2) + F(k,225)
    if s >= 9:
        delta3 = min(F(1,21**2)-F(1,189**2), F(1,15**2)-F(1,27**2))
        delta7 = F(1,21**2)-F(1,441**2)
        upper -= 8*delta3 + witnesses(7,s)*delta7
    return upper

def four_support(s, m=22):
    k = R-1-s-128-m-16-6-4
    if k < 0:
        return None
    three_floor = 99 if s else 33
    return (F(1,1155**2) + F(s,81) + F(128,441) + F(m,121)
            + F(16,625) + 2*(F(1,three_floor**2)+F(1,77**2)+F(1,275**2))
            + F(4,625) + F(k,225))

def main():
    check_fiber_congruences()
    # Several-edge cases.  The active-at-5 and active-at-7 sectors
    # are disjoint; four central-only labels are bounded at 15.
    multiedge = F(1,105**2)+F(72,225)+F(108,441)+F(4,225)+F(R-185,81)
    assert multiedge < 1

    # At least two isolated support primes outside {3,7}.
    excess = F(1,105**2)+F(128,441)+F(R-129,81)-1
    c5 = 16*(F(1,81)-F(1,625))
    c11 = 22*(F(1,81)-F(1,121))
    c13 = 48*(F(1,81)-F(1,169))
    assert c5+c11+c13 > excess
    assert c5+c13 > excess
    assert c11+c13 > excess
    assert c5+20*(F(1,81)-F(1,11**4)) > excess

    # Every possible isolated prime fits the activity reservation;
    # all larger primes already exceed the rank budget.
    results = []
    for q in range(5, R+2, 2):
        if q == 7 or not prime(q):
            continue
        for s in range(0,77,2):
            u = three_support(q,s)
            if u is not None:
                assert u < 1, (q,s,u)
                results.append((u,q,s))
    assert results

    # The sole remaining four-prime configuration is {3,5,7,11},
    # with denominator-eleven class marked.  Forty 5-active labels
    # would already violate the ordinary mass bound.
    dense5 = F(1,1155**2)+F(128,441)+F(22,121)+F(40,625)+F(R-191,81)
    assert dense5 < 1
    four = []
    for m in range(22,R,2):
        for s in range(0,77,2):
            u = four_support(s,m)
            if u is not None:
                assert u < 1, (s,m,u)
                four.append((u,s,m))
    assert four
    print("rank cap:", R)
    print("several-edge mass maximum:", multiedge)
    print("three-prime checks:",len(results),"maximum (mass,q,s):",max(results))
    print("dense-five mass maximum:",dense5)
    print("four-prime checks:",len(four),"maximum (mass,s,m):",max(four))
    print("all enclosing ordinary-mass bounds are strictly below one")

if __name__ == "__main__":
    main()
