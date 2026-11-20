"""Necessary rank and mass tests for small components of odd square data."""
from fractions import Fraction as Q
from itertools import product
from math import ceil
from edges import shapes

def isolated_minimum(p):
    return (4 if p%4==1 else 2)*(p-1)

def compatible(rank,P,positive,G,m,witnesses,output=None):
    if rank <= 217 and any(a==3 and t>=7 for a,t in G):
        return False  # odd:large-ternary-edge, a scoped ordinary-mass theorem.
    # Each edge supplies 2(p-1)(t-1) joint dual pairs and, in
    # each direction, one further pair per occupied primitive phase.
    if any(a!=3 or t not in (5,7) for a,t in G):
        occupied=1+m+sum(v[1] for v in positive.values())
        def edge_cost(a,t):
            if a<5:
                return 4*(a-1)*(t-1)+2*(a+t-2)
            mu=lambda p: 1 if p%8==1 else 3 if p%8==5 else 5
            return min(8*(a-1)*(t-1),
                       4*(a-1)*(t-1)+2*(t-1)*mu(a)+2*(a-1)*mu(t))
        # For several edges use the spherical-summand/Hall reservation.
        # The stronger single-edge phase bound is not added across vertices.
        reserve=(edge_cost(*G[0]) if len(G)==1
                 else sum(4*(a-1)*(t-1)+2*(a+t-2) for a,t in G))
        outside=[p for p in P if all(p not in e for e in G)]
        reserve+=sum(max(0,isolated_minimum(p)-(m if p==3 else positive.get(p,(0,0,0))[2])) for p in outside)
        return occupied+reserve<=rank
    adjacent=[t for _,t in G]
    isolated=[p for p in P if p not in [3]+adjacent]
    low=[p for p,v in positive.items() if v[0]==1]
    base=[(1,105),(m,9)] + [(s,p) for p,(v,s,a) in positive.items()]
    # New nonunit positions, in addition to the booked prime classes.
    base.extend((max(0,isolated_minimum(p)-positive.get(p,(0,0,0))[2]),p*p) for p in isolated)
    base_len=len(base)
    if low:
        p=low[0]
        if len(P)==4:
            base.extend((2,p*q**(2 if (q==3 and m) or (q in isolated and q not in positive) else 1)) for q in P if q!=p)
        else:
            base.append((6,3*p))
    def fits(pairs):
        slots=rank-sum(n for n,b in pairs)
        return slots>=0 and sum(Q(n,b*b) for n,b in pairs)+Q(slots,225)>=1
    isolated_absent=[p for p in isolated if p not in positive]
    edge_min=4+sum(16*(t-1) for t in adjacent)
    if len(P)!=3 and (not low or len(P)==4):
        neutral=0
        for p in isolated_absent:
            if rank-sum(n for n,b in base)+isolated_minimum(p)-edge_min<p*(p-1):
                neutral=max(neutral,4 if not low else 2)
        if neutral:base.append((neutral,min(isolated_absent)**3))
        if len(P)==4 and len(isolated_absent)==2 and not positive:
            if all(rank-1-m-edge_min-sum(isolated_minimum(q) for q in isolated_absent if q!=p)<p*(p-1) for p in isolated_absent):
                # At least eight neutral positions, an overlap count on the principal blocks.
                base.append((8-neutral,15))

    if len(P)==3 and m:
        for p in isolated_absent:
            base.append((witnesses(p,m), (9 if m>=9 else 3)*p*p))
    if len(P)==3 and len(adjacent)==1 and m>=9:
        # All eight three-free labels are distinct from hits omitting other primes.
        # If 9||M, the three-prime support makes gcd one in the principal
        # three-block impossible as all seen degrees divide M/9.
        extra=[(8,27),(witnesses(adjacent[0],m),9*adjacent[0]**2),(16*(adjacent[0]-1),3*adjacent[0])]
        if not fits(base+extra):return False
    if not fits(base):return False
    # Single-coordinate budgets hold at every leaf even if three has two neighbors.
    leaf=[]
    for t in adjacent:
        n=16*(t-1)
        if t in positive and positive[t][2]==t-1:
            # (t-1)/t^2 of singly active mass outside the sole primitive
            # square class at denominator t; tested together with mixed activity.
            n=max(n,4*(t-1)+((t-1)*225//(t*t)))
            # A pure higher t-power can be smaller than 3t only at t=3.
        leaf.append((n,3*t))
    if not fits(base+leaf+[(4,15)]):return False
    if len(P)==3 and m:
        gram=[(witnesses(t,m),(9 if m>=9 else 3)*t*t) for t in adjacent]
        if not fits(base+leaf+gram):return False
    if len(adjacent)!=1:
        avail=rank-sum(n for n,b in base)-4
        if avail<0:return False
        sets=[]
        for t in adjacent:
            w=t-1
            u=[]
            for n in range(8,avail//w+1):
                high=3 if n%3 else (1 if n%9 else 0)
                u.append((n*w,Q((n-high)*w,(3*t)**2)+Q(high*w,(9*t)**2)))
            sets.extend((u,u))
        booked=sum(Q(n,b*b) for n,b in base)+Q(4,225)
        return any(sum(n for n,v in row)<=avail and booked+sum(v for n,v in row)+Q(avail-sum(n for n,v in row),225)>=1 for row in product(*sets))

    t=adjacent[0];d=3*t
    if not fits(base+[(4+16*(t-1),d)]):return False
    projection = None
    if t in positive and positive[t][0] == 2:
        _, total, active = positive[t]
        projection = (total-active+Q((t+1)*active,(t-1)**2))/t**2
        booked_mass = Q(1,105**2) + sum(Q(n,b*b) for n,b in base[base_len:])
        need2 = max(0, 2*ceil((projection-booked_mass)*t**4/2))
        cap = rank-sum(n for n,b in base)-32
        if need2>cap:return False
        if active==t-1 and Q(need2,t**4)+Q(cap-need2,225)<Q(t-1,t*t):
            return False
    slots=rank-sum(n for n,b in base)

    upper=sum(Q(n,b*b) for n,b in base)
    found=False
    high_room=(upper+Q(slots,225)-1)/(Q(1,225)-Q(1,(3*d)**2))
    max_high=high_room.numerator//high_room.denominator
    for z,x,y,zc,xc,yc in shapes(t,slots,max_high):
        expensive=zc+xc+yc
        mass=upper+Q(z+x+y-expensive,d*d)+Q(expensive,(3*d)**2)+Q(slots-z-x-y,225)
        if mass<1:continue
        if projection is not None:
            other_slots=slots-z-x
            if other_slots<need2:continue
            # Mixed solely-active labels can meet the projection only at 3t².
            if booked_mass+Q(y,(3*t*t)**2)+Q(other_slots-y,t**4)<projection:continue
            both=min(yc,need2)
            own=Q(need2-both,t**4)+Q(other_slots-need2-(yc-both),225)+Q(yc,(3*d)**2)
            if active==t-1 and own<Q(t-1,t*t):continue
        if len(P)==3 and m:
            # The three-free positions and the witnesses omitting t have
            # three-primary twists, disjoint from every t-active reservation.
            if m>=9 and rank-(sum(n for n,b in base)+z+y)<8+witnesses(t,m):
                continue
            required=witnesses(t,m)
            witness_floor=(9 if m>=9 else 3)*t*t
            losses=[]
            for j in range(min(required,x)+1):
                if slots-z-x-y>=required-j:
                    i=min(xc,j)
                    losses.append(Q(i,(3*d)**2)+Q(j-i,d*d)+Q(required-j,225)-Q(required,witness_floor**2))
            if not losses or mass-min(losses)<1:continue
        if output is None:return True
        output.add((z,x,y,zc,xc,yc))
        found=True
    return found

