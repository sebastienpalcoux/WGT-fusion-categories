"""Bessel bounds for a complete class of denominator-nine rows."""
from fractions import Fraction as Q
from full_blocks import count, neutral_pair_floors
from math import prod

def ceiling(rank,support,nines,shape,marked,witnesses):
    z,x,y,*_=shape
    isolated=[p for p in support if p not in (3,5) and p not in marked]
    activity=sum(2*(p-1) for p in isolated)
    root=(27*25 if len(support)==3 else 45)*prod(p**3 for p in isolated)*prod(marked)
    outside=2*len(isolated)+1+nines+sum(v[1] for v in marked.values())+z+y+activity
    extra=Q(1,root**2)
    if len(support)>=4:
        parts={3:9,5:5,**{p:p for p in marked},**{p:p**3 for p in isolated}}
        floors=neutral_pair_floors(support,marked,parts,nines,nine=True)
        extra+=2*sum(Q(1,b*b) for b in floors)
        outside+=2*len(floors)
    else:
        q,=isolated
        w5,wq=(witnesses(t,nines) for t in (5,q))
        full_q=max(wq,count(q,nines))
        outside+=full_q
        extra+=Q(full_q,(9*q**3)**2)
        # The other two witness sets may overlap the x mixed positions.
        values=[]
        for i in range(0,min(x,w5)+1,2):
            for j in range(0,min(x-i,8)+1,2):
                spare=rank-outside-w5-8-x+i+j
                if spare>=0:
                    values.append(Q(spare,27**2)+Q(w5,225**2)+Q(j,135**2)
                                  +Q(8-j,27**2)+Q(x-i-j,45**2))
        return extra+max(values,default=0)
    return extra+Q(x,45**2)+Q(rank-outside-x,27**2)

def possible(nines,c,i,extra):
    # i records the unsigned primitive square-class imbalance.
    estimate=81*extra
    # Use the real part as well: its projected squared norm is exact.
    C=nines-3*c
    H=nines+estimate
    u=min(estimate,max(-estimate/2,H/4-C))
    if 2*(C+u)**2>H*(estimate+u):
        return False
    if 3*i*i>estimate**2:
        return False
    length2=(nines-3*c)**2+3*i*i
    margin=length2-nines*estimate
    return margin<=0 or margin**2<=4*estimate*estimate*length2
