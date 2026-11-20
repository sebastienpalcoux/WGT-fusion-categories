"""Primitive p-twist reservation from the Legendre component of the orbit sum."""
from math import ceil
from fractions import Fraction as Q

def two_low_possible(rank, nines, marked, shape):
    """Four support primes 3,5,p,q, with both p and q of valuation one."""
    p,q = sorted(marked)
    z,x,y,zc,xc,yc = shape
    # A witness omitting p,q serves both reservations; the others are disjoint.
    floors = (9*p,5*p,9*q,5*q,p*q)
    booked = 1 + nines + sum(m for _,m,_ in marked.values()) + 10 + z+x+y
    if booked>rank:
        return False
    mass = Q(1,(9*5*p*q)**2)+2*sum(Q(1,b*b) for b in floors)
    mass += Q(nines,81) + sum(Q(m,t*t) for t,(_,m,_) in marked.items())
    high = zc+xc+yc
    return mass + Q(rank-booked+z+x+y-high,225) + Q(high,2025)>=1

