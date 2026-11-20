from fractions import Fraction as F
from math import lcm, prod
from dense_phase import yintervals
from full_blocks import count as free_count, neutral_pair_floors

def possible209(P,marked,s,shape,p,witnesses,explain=False):
    r=209;z,x,y,zh,xh,yh=shape
    if p!=7 or 5 not in P or len(P) not in (3,4) or zh or z not in(32,48):
        return True
    nfiber=2 if z==32 else 3
    if (x//2-3*nfiber)%8 or x!=6*nfiber:
        # Higher feasible populations require individual denominator choices.
        return (x//2-3*nfiber)%8==0
    # x=6*n_fiber => 3 pairs per fiber, all at denominator 45.
    J=F(z,225); A=J/24; Cj=F(7,8) if z==32 else F(1)
    Ca=F(3,2) if z==32 else F(1)
    marked_m=sum((F(m,p*p) for p,(v,m,a) in marked.items()),F())
    marked_g=sum((F(a,p*(p-1)) for p,(v,m,a) in marked.items()),F())
    fixed=1+s+sum(m for _,m,a in marked.values())+z+x+y
    outside=[t for t in P if t not in (3,5,7) and t not in marked]
    if fixed+42+sum(2*(t-1) for t in outside)>r:return False
    qcap=r-fixed-sum(2*(t-1) for t in outside)
    if qcap>=84:return True
    hs={3:9,5:5,7:343,**{p:p for p in marked},**{t:t**3 for t in outside}}
    # A_cap(t)<2t(t-1) for t>=11 throughout this invocation.
    root=prod(hs.values())
    if len(P)==3:
        q5=witnesses(5,s);q7=witnesses(7,s);fq=max(q7,free_count(7,s))
        # All x positions have b=45, unlike the 5-free and 3-free witnesses.
        groups=[(q7,9*343),(fq-q7,343),(q5,225),(8,27)]
        root=27*25*343
    else:
        groups=[(2,b) for b in neutral_pair_floors(P,marked,hs,s)]
    base_slots=fixed+sum(m for m,b in groups)+sum(2*(t-1) for t in outside)
    # For the lower mass use only contributions fixed exactly.
    M0=F(s,81)+J+A+marked_m
    extra_ordinary=sum((F(m,b*b) for m,b in groups),F())+F(1,root*root)
    # Outside active bookings are disjoint; their ordinary and Gauss rates are bounded above.
    extra_ordinary+=sum(F(2*(t-1),t**4) for t in outside)
    extra_charge=sum(F(2*t,t**4) for t in outside)+(F(12,27**2) if len(P)==3 else 0)
    lam=F(64,45); best=[]
    for u in range(qcap//12+1):
        for v in range(2):
            n7=12*u+42*v
            if not 42<=n7<=r-base_slots:continue
            spare=r-base_slots-n7
            for uh in range(u+1):
                for vh in range(v+1):
                    m7=F(12*uh+42*vh,49**2)
                    g7=F(14*uh+42*vh,49**2)
                    hm7=F(n7-12*uh-42*vh,343**2)
                    hg7=F(14*(u-uh)+42*(v-vh),343**2)
                    for c in range(s//2+1):
                        G0=F(c,27)+Cj*J+Ca*A+marked_g+g7
                        for blo,bhi in yintervals(y,yh,z):
                            lo=max(blo,1-M0-m7-extra_ordinary-hm7-F(spare,225))
                            hi=min(bhi,1-M0-m7)
                            lo=max(lo,(1-F(1,root)-G0-extra_charge-hg7-F(spare,486))*F(4,5))
                            hi=min(hi,(1+F(1,root)-G0)*F(4,5))
                            lo=max(lo,((1+lam-F(lam,root)-M0-m7-G0*lam-extra_ordinary-hm7
                                          -lam*(extra_charge+hg7)-F(spare,225))/(1+lam*F(5,4))))
                            if lo<=hi:
                                if not explain:return True
                                best.append((n7,uh,vh,c,spare,lo,hi))
    return best if explain else False


def dense_possible(rank,support,marked,nines,shape,prime,witnesses):
    """Supplementary 209 test after the original dense_budget has retained."""
    if (rank!=209 or prime!=7 or nines<9 or support[0:2] not in ((3,5),[3,5])
            or any(v!=(1,2*p,2*(p-1)) for p,v in marked.items())
            or not ({3,5,7}<=set(support)) or len(support) not in (3,4)):
        return True
    z,x,y,zh,xh,yh=shape
    if zh:return True
    if z==64:
        if len(support)!=3 or y!=8 or yh or nines<42:return True
        if rank-1-nines-z-x-y>=84:return True
        from dense64 import try64
        return bool(try64(nines,x,xh,witnesses))
    if z not in (32,48) or (len(support)==3 and nines<42):
        return True
    remaining=possible209(support,marked,nines,shape,prime,witnesses,explain=True)
    if isinstance(remaining,bool):return remaining
    for n7,u,v,c,slots,lo,hi in remaining:
        if (set(support)!={3,5,7,11} or 11 not in marked or slots>14
                or (n7,u,v) not in ((42,0,1),(42,0,0),(48,3,0))):
            return True
        from dense_phase import y_five_pair_refinement
        intervals=(y_five_pair_refinement() if z==32 and y==40 and lo<hi
                   else ((lo,hi),))
        for left,right in intervals:
            a,b=max(lo,left),min(hi,right)
            if a>b:continue
            J=F(z,225);A=J/24
            common=1-F(c,27)-(F(7,8) if z==32 else 1)*J-(F(3,2) if z==32 else 1)*A-F(2,11)
            common-=F(14*u+42*v,49**2)
            extra=(42 if (n7,u,v)==(42,0,0) else 14 if n7==48 else 0)
            seven=( [(0,0,343)] if not extra else [
                (F(extra,343**2),F(extra,343**2),343**2//7),
                (F(extra,2401**2),F(extra,2401**2),2401**2//7),
                (0,F(extra,16807**2),16807**2//7)])
            for seven_lo,seven_hi,seven_part in seven:
                if short_pure_charge(common-F(5,4)*b-seven_hi,
                                     common-F(5,4)*a-seven_lo,slots,9*5*11*seven_part):
                    return True
    return False


def short_pure_charge(lo,hi,slots,root):
    # Actual packets of pure p-power denominator B contribute pc/B² and
    # use at least (p-1)c labels. At most 14 slots in applications.
    levels=sorted((b,p) for p,bs in ((3,(27,81,243,729,2187)),(5,(25,125,625,3125))) for b in bs)
    def visit(i,k,total,L):
        if total>hi+F(1,L):return False
        if i==len(levels):
            return total+F(3*k,2*6561**2)>=lo-F(1,L)
        b,p=levels[i]
        if total+F(3*k,2*b*b)<lo-F(1,L):return False
        for c in range(k//(p-1)+1):
            if visit(i+1,k-c*(p-1),total+F(p*c,b*b),L if not c else lcm(L,b)):
                return True
        return False
    return visit(0,slots,F(),root)
