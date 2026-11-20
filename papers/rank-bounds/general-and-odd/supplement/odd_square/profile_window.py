"""Individual 15 and 45 counts and mixed mass in a fixed joint family at 15."""
from fractions import Fraction as F
from functools import lru_cache
from dense_phase import groups
from highprofile import merge


@lru_cache(None)
def tail(energy, first, n, p, floor):
    """Ordinary tail mass from its first two moments and coefficient lattice."""
    if not n:
        return ((F(), F()),) if not energy and not first else ()
    if energy <= 0:
        return ()
    lower=first*first/energy
    upper=min(F(n,floor*floor),energy)
    # (b-u)(b-u-p) >= 0 for all b in -1+pZ, u in -1+pZ.
    # Every one of the displayed inequalities is valid separately.
    for j in (-6,-4,-3,-2,-1,0,1,2,3,4,6):
        u=p*j-1;d=u*(u+p)
        threshold=((2*u+p)*first-energy)/d
        if d>0:
            lower=max(lower,threshold)
        else:
            upper=min(upper,threshold)
    return ((lower,upper),) if lower<=upper and upper>0 else ()


@lru_cache(None)
def groups_even(n,k,p,E,d):
    """First moment, squared moment, count of even coefficients at one level."""
    from math import isqrt
    mx=isqrt(E*k*k)
    states=[set() for _ in range(n+1)];states[0].add((0,0,0))
    for b in range(-mx,mx+1):
        if (b+1)%p:continue
        even=int((b+1)%2==0)
        tile=d if even else 1
        for j in range(tile,n+1):
            states[j].update((v+tile*b,e+tile*b*b,count+tile*even)
                             for v,e,count in states[j-tile] if e+tile*b*b<=E*k*k)
    return tuple(states[n])


@lru_cache(None)
def allowed_even(p,J,C,h):
    d=(p-1)&-(p-1)
    return tuple(e for e in range((-J)%8,h+1,8) if e%d==0
                 and (J!=(p-1)*C or 2*d*J >= (p-1)*(2*d+e*(d-1))))


@lru_cache(None)
def dyadic_moments(n,even,d):
    """Modulo 32 first and squared moments of free odd-square weights.

    An even coefficient occurs in d equal-weight, equal-coefficient
    copies; all other coefficients occur singly. A cofactor may have
    any odd square inverse modulo 32.
    """
    if even<0 or even>n or even%d:return frozenset()
    if not n:return frozenset({(0,0)})
    if even<n:
        tile=1;bs=range(0,32,2) # odd a gives even b=a-1
        prev=dyadic_moments(n-1,even,d)
    else:
        tile=d;bs=range(1,32,2)
        prev=dyadic_moments(n-d,even-d,d)
    shifts={((tile*b*w)%32,(tile*b*b*w)%32) for b in bs for w in (1,9,17,25)}
    return frozenset(((m+t)%32,(e+v)%32) for m,e in prev for t,v in shifts)


@lru_cache(None)
def tail_by_parity(energy,first,n,p,floor,emin,emax):
    base=tail(energy,first,n,p,floor)
    if not n or not base:return base
    lower,upper=base[0]
    max_even=F(emax,floor*floor)
    max_odd=F(n-emin,floor*floor)
    square=(p-1)**2
    upper=min(upper,(energy+(square-1)*max_even)/square)
    for parity,opposite in ((0,max_odd),(1,max_even)):
        for j in range(-3,4):
            u=p*(2*j+parity)-1
            coeff=u*(u+2*p)
            v=((2*u+2*p)*first-energy-p*p*opposite)/coeff
            if coeff>0:lower=max(lower,v)
            else:upper=min(upper,v)
    return ((lower,upper),) if lower<=upper and upper>0 else ()


@lru_cache(None)
def one(p,J,C,h,v,w):
    return one_deep(p,J,C,h,(v,w))


@lru_cache(None)
def shallow_window(z,side,total,at15=None,at45=None):
    """Union of possible total ordinary masses, an empty tuple signifies exclusion.

    z=32/64 primitive or 48 order-45, joint denominators all fifteen.
    side='three' or 'five' describes which prime is active, not the Weyl prime.
    at15 and at45 are exact full label counts, omitted for unconstrained.
    """
    if z not in (32,48,64) or side not in ('three','five'):
        raise ValueError((z,side))
    if side=='five':
        p=3;m=2;nphase=2
        J,C={32:(4,2),48:(6,0),64:(8,4)}[z]
    else:
        p=5;m=3 if z==48 else 1;nphase=1 if z==48 else 2
        J,C=(16,4) if z==64 else (8,2)
    unit=2*m
    if any(a is not None and (a<0 or a%unit) for a in (total,at15,at45)):return ()
    K=total//unit
    v0=None if at15 is None else at15//unit
    w0=None if at45 is None else at45//unit
    intervals=[]
    def choices(value,limit):
        return range(limit+1) if value is None else (value,) if 0<=value<=limit else ()
    if nphase==1:
        for v in choices(v0,K):
            for w in choices(w0,K-v):
                intervals.extend(one(p,J,C,K,v,w))
    else:
        for h in range(1,K//2+1):
            rem=K-h
            for v in range(min(h,K if v0 is None else v0)+1):
                for w in range(min(h-v,K if w0 is None else w0)+1):
                    aa=one(p,J,C,h,v,w)
                    if not aa:continue
                    for v2 in choices(None if v0 is None else v0-v,rem):
                        for w2 in choices(None if w0 is None else w0-w,rem-v2):
                            bb=one(p,J,C,rem,v2,w2)
                            intervals.extend((a+c,b+d) for a,b in aa for c,d in bb)
    factor=F(unit,225)
    return tuple((factor*l,factor*u) for l,u in merge(intervals))

@lru_cache(None)
def one_deep(p,J,C,h,counts):
    """Sharper window with fixed counts at all observed levels in this phase."""
    known=(1,3,5,9)[:len(counts)]
    floor={2:5,3:9,4:15}[len(counts)]
    r=h-sum(counts)
    if r<0 or h<=0:return ()
    d=(p-1)&-(p-1)
    E=J+2*C
    base=sum((F(n,k*k) for k,n in zip(known,counts)),F())
    allowed=allowed_even(p,J,C,h)
    if not allowed:return ()
    if J==(p-1)*C and allowed==(0,):
        V=F(C,p+1)
        if V.denominator%2==0 or h%8!=(V.numerator*pow(V.denominator,-1,8))%8:
            return ()
        return ((V,V),) if base<=V<=base+F(r,floor*floor) and (not r or base<V) else ()
    # At each step use the total remaining weight, including the exact
    # multiplicities at the still unexamined levels.
    scale=2025
    budget=scale*E
    target=-scale*C
    states={(0,0,0)}
    for ix,(k,n) in enumerate(zip(known,counts)):
        if not n:continue
        unit=scale//(k*k)
        unpaid=sum(v*(scale//(b*b)) for v,b in zip(counts[ix+1:],known[ix+1:]))
        future=(unpaid*floor**2+r*scale)
        terms=groups_even(n,k,p,E,d)
        if not terms:return ()
        if len(states)*len(terms)>35000 and len(counts)>2:
            aa=one(p,J,C,h,counts[0],counts[1])
            return merge((max(a,base),min(b,base+F(r,floor**2))) for a,b in aa
                         if max(a,base)<=min(b,base+F(r,floor**2)))
        left=r+sum(counts[ix+1:])
        feasible={e for e in range(max(allowed)+1) if any(e<=a<=e+left for a in allowed)}
        states={(m+t*unit,e+v*unit,c+a) for m,e,c in states for t,v,a in terms
                if c+a in feasible and e+v*unit<=budget and
                (target-m-t*unit)**2*floor**2 <= (budget-e-v*unit)*future}
        if not states:return ()
    intervals=[]
    for m,e,c in states:
        mr=((target-m)*pow(scale,-1,32))%32
        er=((budget-e)*pow(scale,-1,32))%32
        ec=[a-c for a in allowed if c<=a<=c+r and (mr,er) in dyadic_moments(r,a-c,d)]
        if not ec:continue
        for a,b in tail_by_parity(F(budget-e,scale),F(target-m,scale),r,p,floor,min(ec),max(ec)):
            intervals.append((base+a,base+b))
    return merge(intervals)


@lru_cache(None)
def phase_window(z,side,total,at15=None,at45=None,at75=None,at135=None):
    """Count-conditioned mass intervals, including two optional exact higher levels."""
    if at75 is None:
        return shallow_window(z,side,total,at15,at45)
    if at15 is None or at45 is None:
        return shallow_window(z,side,total,at15,at45)
    if side not in ('three','five') or z not in (32,48,64):
        raise ValueError((z,side))
    if side=='five':
        p=3;m=2;nphase=2;J,C={32:(4,2),48:(6,0),64:(8,4)}[z]
    else:
        p=5;m=3 if z==48 else 1;nphase=1 if z==48 else 2;J,C=(16,4) if z==64 else (8,2)
    unit=2*m
    counts=(at15,at45,at75) if at135 is None else (at15,at45,at75,at135)
    if any(n<0 or n%unit for n in (total,)+counts) or sum(counts)>total:return ()
    K=total//unit
    counts=tuple(n//unit for n in counts)
    intervals=[]
    if nphase==1:
        intervals=one_deep(p,J,C,K,counts)
    else:
        from itertools import product
        # Larger allocations get the shallow necessary condition.
        size=1
        for n in counts:size*=n+1
        if size*(K-sum(counts)+1)>3500:
            return shallow_window(z,side,total,at15,at45)
        for left in product(*(range(n+1) for n in counts)):
            right=tuple(n-a for n,a in zip(counts,left))
            for h in range(max(sum(left),1),min(K//2,K-sum(right))+1):
                aa=one_deep(p,J,C,h,left)
                if not aa:continue
                bb=one_deep(p,J,C,K-h,right)
                intervals.extend((a+c,b+d) for a,b in aa for c,d in bb)
    factor=F(unit,225)
    return tuple((l*factor,u*factor) for l,u in merge(intervals))
