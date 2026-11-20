from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from math import prod,lcm,isqrt
from bisect import bisect_left,bisect_right
from nine_gram import possible as nine_possible, ceiling as nine_ceiling
from highprofile import merge
from profile_window import phase_window

@lru_cache(None)
def overlaps(z, x):
    """Maximum feasible 3-free and 5-free overlaps in an exact rigid population."""
    if z not in (32,48) or x not in ((12,28) if z==32 else (18,66)):return None
    def fib(h):
        # Reciprocal squares of 3,5-smooth cofactors in one third-root fiber.
        def visit(left, residual, low, n3, n5, both):
            if not left:return {(n3,n5,both)} if not residual else set()
            if residual<=0:return set()
            cap=isqrt(left*residual.denominator//residual.numerator);values=[];a=1
            while a<=cap:
                b=a
                while b<=cap:
                    if b>=low:values.append(b)
                    b*=5
                a*=3
            out=set()
            for b in values:
                for j in range(1,min(left,int(residual*b*b))+1):
                    out.update(visit(left-j,residual-F(j,b*b),b+1,n3+(j if b%9==0 else 0),n5+(j if b%15==0 else 0),both+(j if b%45==0 else 0)))
            return out
        return visit(h,F(1,3),3,0,0,0)
    if z==48:
        return {(6*a,6*b,6*(a+b-c)) for a,b,c in fib(x//6)}
    return {(2*(a+A),2*(b+B),2*(a+A+b+B-c-C)) for a,b,c in fib(3) for A,B,C in fib(x//2-3)}
from highprofile import sums
from full_blocks import count as free_count, neutral_pair_floors

@lru_cache(None)
def first_pure(spare, three_free, scale):
    """Charges, minimum counts, for the first pure levels at 25 and 27."""
    entries=[]
    for five in range(spare//4+1):
        for three in range((spare-4*five+three_free)//2+1):
            outside3=max(0,2*three-three_free)
            extra=spare-4*five-outside3
            modulus=5*five*(scale//625)+3*three*(scale//729)
            low=4*five*(scale//625)+2*three*(scale//729)
            up=4*five*(scale//625)+outside3*(scale//729)
            future=spare-4*five+three_free-2*three
            entries.append((modulus,low,up,extra,future))
    entries.sort()
    return tuple(a[0] for a in entries),entries


@lru_cache(None)
def box(z,side,total,high):
    if high:
        return merge((a,b) for n in range(total-high+1) for a,b in phase_window(z,side,total,n))
    return phase_window(z,side,total)


def sparse_possible(row,shape,witnesses,seven=None):
    P,classes,edges,s=row
    R=209
    z,x,y,zh,xh,yh=shape
    if zh or z not in (32,48,64,96,128) or s<9 or tuple(map(tuple,edges))!=((3,5),):return True
    marked={p:tuple(v) for p,v in zip(P[1:],classes) if v[0]}
    if any(v!=(1,2*p,2*(p-1)) for p,v in marked.items()):return True
    Q=[p for p in P if p not in (3,5) and p not in marked]
    if len(P) not in (3,4,5) or not Q or (len(P)==3 and len(Q)!=1):return True
    if len(P)==5 and 7 not in Q:return True
    fixed=1+s+z+x+y+sum(v[1] for v in marked.values())
    if R<fixed+sum(2*(p-1) for p in Q):return False
    hs={3:9,5:5,**{p:p**3 for p in Q},**{p:p for p in marked}}
    pair_floor={p:p*p for p in Q}
    if seven is not None:
        if Q!=[7]:raise ValueError("seven signature requires exactly one unmarked isolated prime")
        total,small,middle=seven
        if total not in (12,24,36) or min(small,middle)<0 or small%12 or middle%12 or small+middle>total:return False
        largest=2401 if total>small+middle else 343 if middle else 49
        hs[7]=largest**2//7
        pair_floor[7]=49 if small else 343 if middle else 2401
    root=prod(hs.values())
    if len(P)==3:
        q,=Q;root=27*25*hs[q]
        wq=witnesses(q,s);w5=witnesses(5,s);fq=max(wq,free_count(q,s)) if q<=13 else wq
        oms=[]
        for o3 in range(0,7,2):
            for o5 in range(0,w5+1,2):
                if o3+o5>x:continue
                r3=6-o3; r5=w5-o5
                base=(F(r3,27**2)+F(2,(27*q**2)**2)+F(r5,225**2)
                      +F(wq,(9*hs[q])**2)+F(fq-wq,hs[q]**2))
                n=8-o3+r5+fq
                charge=F(3*r3,2*27**2)
                # A neutral five-free pair either occurs among the selected
                # Gram witnesses or uses two further disjoint labels.
                if r5>=2:
                    oms.append((o3,o5,n,base-F(2,225**2)+F(2,(225*q**2)**2),charge))
                oms.append((o3,o5,n+2,base+F(2,(25*q**2)**2),charge))
    else:
        ns=neutral_pair_floors(P,marked,hs,s)
        oms=[(0,0,2*len(ns),sum(F(2,b*b) for b in ns),F())]
    ovs=overlaps(z,x) if not zh else None
    if ovs is not None:
        oms=[(i,j,n,m,g) for i,j,n,m,g in oms if any(i<=a and j<=b and i+j<=c for a,b,c in ovs)]
    M=F(s,81)+F(z,225)+sum(F(v[1],p*p) for p,v in marked.items())
    gp=sum(F(v[2],p*(p-1)) for p,v in marked.items())
    lam=F(64,45)
    # Every remaining Q-activity is sparse, handled separately otherwise.
    # For fixed total reserve a degree-p-divisible pure inactive pair when
    # the proved upper quota is small enough to make it disjoint.
    pair={p:(2 if R-fixed-sum(2*(t-1) for t in Q)<2*p*p*(1-F(1,p-1)-F(1,6*p)) else 0) for p in Q}
    # One common multiple makes every stored endpoint exact, including the moments.
    # All divisions in the continuous interval intersection round outward.
    scale=2**8*3**16*5**16*prod(p**14*(p-1) for p in P if p>5)
    def units(v):
        q,r=divmod(v.numerator*scale,v.denominator)
        assert r==0
        return q
    vac=scale//root
    M=units(M)
    nine=nine_ceiling(R,P,s,shape,marked,witnesses)
    cs=[c for c in range(s//2+1) if nine_possible(s,c,c%2,nine)]
    if not cs:return False
    for a,b,c in product(range(0,z//16+1,2),repeat=3):
        if b>c or 16*a+24*(b+c)!=z:continue
        gj=gp+F(14*a+24*(b+c),225)
        X=sums(5,tuple((g,4*t,t,order) for g,t,order in [(1,a,3),(1,a,3),(3,b,9),(3,c,9)]),x,xh)
        Y=sums(3,((2,2*a+3*(b+c),a,5),)*2,y,yh)
        if not X or not Y:continue
        def down(q):return (q*scale).__floor__()
        def up(q):return (q*scale).__ceil__()
        if z in (32,48,64):
            # These supports have only one three-primary order.
            x_boxes=tuple((down(l),up(h),down(l) if z==48 else 0,up(h) if z==48 else 0)
                          for l,h in box(z,'three',x,xh))
            y_boxes=tuple((down(l),up(h)) for l,h in box(z,'five',y,yh))
        else:
            x_boxes=tuple(tuple(units(q) for q in interval) for interval in X)
            y_boxes=merge((units(yl),units(yu)) for yl,yu,_,_ in Y)
        if not x_boxes or not y_boxes:continue
        ends=tuple(hi for lo,hi in y_boxes)
        gj=units(gj)
        for o3,o5,n,ne,chg in oms:
            if seven is not None:
                # This profile and all current omission groups are disjoint
                # from the pure seven shell.
                quota=R-fixed-n-seven[0]
                pair[7]=2 if quota<2*7**2*(1-F(1,6)-F(1,42)) else 0
            cap=R-fixed-n-sum(pair.values())
            if cap<sum(2*(p-1) for p in Q):continue
            mx=units(F(x-o3-o5-max(0,xh-o3-o5),225)+F(max(0,xh-o3-o5),45**2)+F(o3,135**2)+F(o5,225**2))
            ntr=units(ne+F(1,root**2)+sum(F(pair[p],pair_floor[p]**2) for p in Q))
            chg=units(chg)
            for steps in product(*(range(2*(p-1),min(cap,(p-1)**2)+1,2*(p-1)) for p in Q)):
                if seven is not None and steps!=(seven[0],):continue
                spare=cap-sum(steps)
                if spare<0:continue
                for starts in product(*(range(0,h+1,2*(p-1)) for p,h in zip(Q,steps))):
                    if seven is not None and starts!=(seven[1],):continue
                    mq=sum(t*(scale//p**4) for p,t in zip(Q,starts))
                    gq=sum(t*p*(scale//((p-1)*p**4)) for p,t in zip(Q,starts))
                    mh=sum((h-t)*(scale//p**6) for p,h,t in zip(Q,steps,starts))
                    gh=sum((h-t)*p*(scale//((p-1)*p**6)) for p,h,t in zip(Q,steps,starts))
                    exact_shell=0
                    if seven is not None:
                        total,small,middle=seven
                        mq+=middle*(scale//343**2)
                        gq+=7*middle*(scale//(6*343**2))
                        remaining=total-small-middle
                        mh=remaining*(scale//2401**2)
                        gh=7*remaining*(scale//(6*2401**2))
                        if total==12 and small+middle==total:
                            exact_shell=pair[7]*(scale//pair_floor[7]**2)
                    mass=M+mq
                    r3=6-o3 if len(P)==3 else 0
                    charges,pure=first_pure(spare,r3,scale)
                    for xl,xu,x9l,x9u in x_boxes:
                        xu=min(xu,mx)
                        if xl>xu:continue
                        mhi=scale-mass-xl-exact_shell
                        mlo=scale-mass-xu-ntr-mh-spare*(scale//225)
                        if mlo>ends[-1] or mhi<y_boxes[0][0]:continue
                        for c9 in cs:
                            g=gj+gq+c9*(scale//27)
                            # Positive signs express residual charge needed after X.
                            rlow=xl if z==48 else (3*xl-x9u)//2
                            rhigh=xu if z==48 else -((x9l-3*xu)//2)
                            glo=scale-vac-g-rhigh-gh-chg
                            ghi=scale+vac-g-rlow
                            lo=max(mlo,4*(glo-spare*(scale//486))//5,(45*mlo+64*glo)//125)
                            hi=min(mhi,-((-4*ghi)//5))
                            pos=bisect_left(ends,lo)
                            if lo>hi:continue
                            for l,h in y_boxes[pos:]:
                                if l>hi:break
                                ylo,yhi=max(l,lo),min(h,hi)
                                # Complete Gauss residual from pure new levels.
                                lower=glo+chg-(5*yhi+3)//4
                                upper=ghi-5*ylo//4
                                q0=bisect_left(charges,lower-(spare+r3)*(scale//4374))
                                q1=bisect_right(charges,upper)
                                for q,pl,pu,k,tail_count in pure[q0:q1]:
                                    lm=scale-mass-xu-ntr-mh-pu-k*(scale//225)
                                    hm=mhi-pl
                                    lg=glo+chg-q
                                    hg=ghi-q
                                    # Future pure labels at 81 or 125 have charge <= 1/4374.
                                    yl=max(ylo,lm,4*(lg-tail_count*(scale//4374))//5)
                                    yu=min(yhi,hm,-((-4*hg)//5))
                                    if yl<=yu:return True
    return False
