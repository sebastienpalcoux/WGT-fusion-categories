#!/usr/bin/env python3
"""Exact multiplicity bounds for odd nonsquare dimension, below rank 209."""
from bisect import bisect_left
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil, prod

RANK = 203  # the last admissible rank in the four-prime ternary case


def upper(slots, vacuum, demands):
    used = sum(n for n, _ in demands)
    if used > slots:
        return -1
    return vacuum + Q(slots-used, 225) + sum(Q(n,max(b,15)**2) for n,b in demands)


@lru_cache(None)
def signatures(m):
    """Twice neutral minus primitive, and positive minus negative primitive."""
    values = set()
    for higher in range(0,m+1,6):
        for neutral in range(0,m-higher+1,2):
            prim = m-higher-neutral
            for plus in range(0,prim+1,2):
                values.add((2*neutral-prim,2*plus-prim))
    return sorted(values)


@lru_cache(None)
def signed_blocks(p,m):
    """Signed sums at a pure allowed prime after averaging its phases."""
    if m==6*(p-1):
        return (Q(-4,p*p),Q(4,p*p))
    if p in (7,11) and m==6*p:
        return (Q(-8,p*p),0,Q(8,p*p))
    if m<6*(p-1):
        orbit = (p-1)&-(p-1)
        return tuple(Q(2*i,p*p) for i in range(-m//3,m//3+1,2*orbit))
    # Each ternary block contributes either ±2, ±2/(p-1) or zero, divided by p².
    limit = (m//3)*(p-1)
    return tuple(Q(2*i,(p-1)*p*p) for i in range(-limit,limit+1))


def saturated_fibers(energy,slots,vacuum,activity,reserve):
    """The minimal ramified-fiber lemma; all inequalities are outward exact."""
    assert sorted(activity)==[(24,15),(36,21)]
    assert slots==16 and sum(n for n,_ in reserve)==14
    neutral=upper(slots,vacuum,reserve)
    booked=sum(Q(n,b*b) for n,b in activity)
    residual=energy-booked
    assert residual==Q(617,363825)
    # Orbit sizes four and six force all booked labels onto their floors.
    losses=(4*(Q(1,225)-Q(1,625)),6*(Q(1,441)-Q(1,1225)))
    assert all(booked-loss+neutral<energy for loss in losses)
    assert -2*residual>Q(-2,81) and residual<Q(2,81)
    signed=(abs(sum(v)) for v in product(
        (-Q(4,81),Q(0),Q(4,81)),
        (-Q(8,49),Q(0),Q(8,49)),
        (-Q(8,121),Q(0),Q(8,121)),
        (-Q(4,225),Q(4,225)),(-Q(4,441),Q(4,441))))
    gap=min(signed)
    assert gap==Q(96356,12006225)
    assert gap>2*residual+Q(2,3**9*25*7*11)
    return True


def check(ps):
    ds = tuple(sorted({3,9}|set(ps)))
    base = 3*prod(ps)
    vacuum = Q(1,base**2)
    req = {p:6*(p-1) for p in ps}
    stats = [0,0,0,0]  # coarse mass, sector mass, height branches, surviving branches

    @lru_cache(None)
    def possible_full(p,m):
        if m<req[p] or (p==7 and m==42):
            return False
        return m//3!=p*p-p and any(
            j*(p-1)<=m//3 and p*abs(p*j-m//3)<=j*(p-j)+Q((p+1)*j*p,base)
            for j in range(2,p,2))

    def accept(ms,k,energy):
        if energy>upper(k,vacuum,[]):
            return
        stats[0]+=1
        m=dict(zip(ds,ms))
        activity=[]
        # Every allowed prime occurring as a denominator has valuation one.
        low={p for p in ps if m[p]}
        for p in ps:
            if m[p]<req[p]:
                count=max(req[p],ceil(Q(9*(p-1)*m[p],2*p)))
                activity.append((count,3*p))
            elif not possible_full(p,m[p]):
                activity.append((req[p],3*p))
        if energy>upper(k,vacuum,activity):
            return
        stats[1]+=1
        # Simultaneous disjoint-field projection (5 has valuation one).
        if m[3]>=8 and m.get(5,0)>=36:
            return
        if m[3]>=6 and m.get(5,0)>=48:
            return
        if m[3]==18 and m.get(7,0)>=42:
            return
        assert m[3] or m[9]
        for p in ps:
            if m[p]==req[p]:
                assert 1+m[p]+ceil(Q(9*(p-1)*m[p],2*p))+sum(req[s] for s in ps if s!=p)>RANK
        choices=[(activity,False)]
        if m.get(11)==66:
            choices=[(activity,True),(activity+[(60,33)],False)]
        for bookings,eleven_full in choices:
            if energy>upper(k,vacuum,bookings):
                continue
            if eleven_full and ps==(5,7,11):
                if m[5]>=36 or (m[5]==24 and m[3]>=4) or m[7]>=84 or m[3]==18:
                    continue
            finish(ms,k,energy,m,low,bookings,eleven_full)

    def finish(ms,k,energy,m,low,activity,eleven_full):
        slots=k-sum(n for n,_ in activity)
        signs={Q(0)}
        for p in ps:
            signs={a+b for a in signs for b in signed_blocks(p,m[p])}
        signs=sorted(signs)
        weighted=sum(Q(n,(b//3-1)*b*b) for n,b in activity)

        for g in (3,9,27):
            # Strict constant-pair-trace Bessel bound at minimal 3-valuation.
            if g==3 and (m[9] or m[3]>=18):
                continue
            balanced = g==27 or (g==9 and m[9]<=2)
            if balanced and (m[3]%6 or (m[3]==18 and m.get(7)==36)):
                continue
            # The range36 needs ternary balance or no quadratic rows;
            # the pre-existing range42 is unconditional.
            if (eleven_full and ps==(5,7,11) and m[3]==12 and m.get(7,0)>=36
                    and (m[7]>=42 or balanced or (g==9 and m[9]>=4))):
                if g<27:
                    continue
                if m.get(5,0) or m.get(7,0)>42:
                    continue
                g=3**9
            no_quadratic=g==9 and m[3]>3 and m[9]>=4
            if no_quadratic and (m[3]%6 or (m[3]==18 and m.get(7)==36)):
                continue
            if g==9 and energy<sum(Q(n,b*b) for n,b in activity if b//3 in low):
                continue
            lows=low|({3} if g==3 else set())
            reserve=[]
            for a,b in combinations((3,)+ps,2):
                if a in lows or b in lows:
                    floor=(g if a==3 else a)*b
                    if g==3**9 and (a==5 or b==5):
                        floor*=5
                    reserve.append((2 if a==3 or g>3 else 6,floor))
            free=sum(n for n,b in reserve if b%3==0)+sum(m[b] for b in ds if b%g==0)
            if g==9 and m[3]:
                all_free=(9,15,21) if m[3]>3 else (9,15,21,27)
                count=next((f for f in all_free if f>=1+free),None)
                if count is None:
                    continue
                reserve.append((count-1-free,45))
            else:
                reserve.append((max(0,8-free),g))
            if g>3 and slots<4+sum(n for n,_ in reserve):
                # Each witness at a denominator prime to three gives six equal labels.
                pairs={p*q for p,q in combinations(ps,2)}
                reserve=[(n,3*b if b in pairs else b) for n,b in reserve]
            if energy>upper(k,vacuum,activity+reserve):
                continue
            stats[2]+=1
            neutral=upper(slots,vacuum,reserve)
            if (g==3**9 and m[3]==12 and m[5]==0 and m[7]==42
                    and m[9]==6 and m[11]==66 and eleven_full):
                if saturated_fibers(energy,slots,vacuum,activity,reserve):
                    continue
            for d,u in ([(0,0)] if balanced or no_quadratic else signatures(m[3])):
                for e,v in signatures(m[9]):
                    z=Q(9*d+e,81)
                    if not -2*neutral<=z<=neutral:
                        continue
                    error=Q(6,g*base)+min(energy,(2*(neutral+weighted)+z)/3)
                    t=-Q(9*u+v,81)
                    i=bisect_left(signs,t-error)
                    if i<len(signs) and signs[i]<=t+error:
                        stats[3]+=1
                        print("unresolved profile",ps,m,g)
                        return

    def visit(i,ms,k,e):
        if e<=0:
            return
        if i==len(ds):
            accept(ms,k,e)
            return
        b=ds[i]
        if e>vacuum+Q(k,min(b,15)**2):
            return
        step=2 if b%3==0 else 3*((b-1)&-(b-1))
        for m in range(0,min(k,3*b*b-1)+1,step):
            if (b==3 and (m>25 or m==24)) or (b==5 and m>48) or (b in ps and m==3*b*(b-1)):
                continue
            visit(i+1,ms+(m,),k-m,e-Q(m,b*b))

    for p in ps:
        assert 30*(p-1)+sum(req[s] for s in ps if s!=p)>RANK
    visit(0,(),RANK-1,Q(3))
    return tuple(stats)


def global_bounds():
    """Small rational checks in the nonsquare support reduction."""
    assert 135+120*(1-Q(9,49))>208
    assert 189+112*(1-Q(9,25))>208
    assert 171+152*(1-Q(9,25))>208
    assert 12*Q(1600,89)>208 and 6*Q(7056,193)>208
    assert Q(62,25)+Q(143,81)<5
    # Four omission witnesses at each of the two other primes.
    assert (Q(104,25)+Q(60,441)+Q(20,81)+Q(12,625)
            +Q(4,2025)+Q(4,60025)+Q(1,11025**2)<5)
    # Uniform residual activity at3 and7 treats both quinary supports.
    assert (Q(20,9)+Q(20,81)+Q(60,441)+Q(100,121)
            +Q(4,25)+Q(1,1155**2)<5)
    # In the five-prime ternary branch,5,7,13 all need residual activity.
    assert min(1+12*(p-1)+sum(6*(s-1) for s in (5,7,11,13) if s!=p)
               for p in (5,7,13))==217>RANK
    assert (Q(10,9)+Q(24,225)+Q(36,441)+Q(60,121)
            +Q(72,1521)+Q(1,15015**2)<3)
    # The two proper-pair Gram cases, and the field-trace projection.
    x,y,z=Q(4,9),Q(2,7),Q(2,11)
    assert (1-x-y-z)*(16-4*x-64*z)-(-2*x+8*z)**2==Q(-64,693)
    assert Q(37,99)*Q(1,99)-Q(35,99)**2<0
    assert Q(4*6,27)+Q(48,300)==Q(236,225)>1
    # Enlarged full-eleven range: exact residual Gram determinant.
    gram=[[Q((1,2,2)[i]) if i==j else Q(0) for j in range(3)] for i in range(3)]
    for weight,vector in ((Q(4,9),(1,1,0)),(Q(12,49),(1,0,2)),
                          (Q(2,11),(1,-2,-1))):
        for i in range(3):
            for j in range(3):
                gram[i][j]-=weight*vector[i]*vector[j]
    a=gram
    determinant=(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
                 -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
                 +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    assert determinant==Q(-16,539)
    # The three free-trace budgets in the extended ternary balance proof.
    from math import isqrt
    assert [max(5*a+f-a*a for a in range(isqrt(f)+1)) for f in (25,16,7)]==[31,22,13]


if __name__=='__main__':
    global_bounds()
    receipts={
        (5,7,11):(54153,1537,50,0),
        (5,7,13):(28975,315,4,0),
        (5,7,17):(9935,7,0,0),
        (5,7,19):(54713,3,0,0),
        (5,7,23):(54559,0,0,0),
        (5,11,13):(31124,23,0,0),
        (5,11,17):(10847,0,0,0),
        (5,11,19):(57125,0,0,0),
        (5,13,17):(5680,0,0,0),
        (7,11,13):(69829,3,0,0),
        (7,11,17):(25364,0,0,0),
    }
    supports=[ps for ps in combinations((5,7,11,13,17,19,23),3)
              if 3+sum(6*(p-1) for p in ps)<=RANK]
    assert set(supports)==set(receipts)
    for ps in supports:
        found=check(ps)
        print(f"Support {(3,)+ps}: {found}",flush=True)
        if found!=receipts[ps]:
            raise ArithmeticError((ps,found,receipts[ps]))
