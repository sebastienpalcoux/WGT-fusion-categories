"""Primitive p-twist reservation from the Legendre component of the orbit sum."""
from math import ceil
from fractions import Fraction as Q

def primitive_quota(p, active):
    if not active:
        return 0
    # Each square class has a multiple of p-1 labels (Galois and duality).
    imbalance = p-1 if active//(p-1) % 2 else 0
    return 2*ceil(Q(p**4 * active**2, 2*(p*p*active-imbalance**2)))

def marked_possible(rank, support, marked, edges, nines):
    # Shared-column Gram and weighted Weyl bounds: odd:eleven-nine.
    if (marked.get(11)==(1,22,20) and all(11 not in edge for edge in edges) and nines>=64):
        return False
    # Under the certified absence of three-prime twists, the edge cost
    # combined with orbit mass excludes this class below rank225.
    if rank<225 and (3,5) in edges and marked.get(5)==(2,4,4):
        return False
    for p, (v, m, active) in marked.items():
        if v!=2 or not active:
            continue
        booked = 1+nines+sum(t[1] for t in marked.values())
        for a,b in edges:
            # Every joint phase fibre has at least two dual pairs.
            booked += 4*(a-1)*(b-1)
        for t in support:
            if t==p or any(t in e for e in edges):
                continue
            old = nines if t==3 else marked.get(t,(0,0,0))[2]
            quota = (4 if t%4==1 else 2)*(t-1)
            booked += max(0, quota-old)
        if booked+primitive_quota(p,active)>rank:
            return False
    return True

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

def other_edge_possible(rank,support,marked,edges,nines,witnesses):
    """Sparse three-prime components accompanying many denominator-nine rows."""
    if edges==((3,11),):
        return 1+nines+16*(11-1)<=rank
    if edges!=((5,7),):
        raise ValueError("unhandled nonternary edge")
    reserve=172  # odd:edge-reservation; formerly only the first pair per fiber
    if len(support)==3:
        return 1+nines+reserve+witnesses(5,nines)<=rank
    booked=1+nines+reserve+sum(m for _,m,_ in marked.values())
    max_mass=Q(1,105**2)+Q(nines,81)+Q(reserve,35**2)
    return booked<=rank and max_mass+sum(Q(m,p*p) for p,(_,m,_) in marked.items())+Q(rank-booked,225)>=1
