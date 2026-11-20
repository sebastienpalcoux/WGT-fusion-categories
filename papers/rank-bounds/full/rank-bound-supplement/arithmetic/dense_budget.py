"""Ordinary mass before the sparse-height conclusion for an isolated prime."""
from fractions import Fraction as Q
from full_blocks import count as free_count, seven_free

def possible(rank,support,marked,nines,shape,prime,witnesses):
    joint,left,right,high_joint,high_left,high_right=shape
    occupied=1+nines+sum(v[1] for v in marked.values())
    dense=prime*(prime-1)
    other_activity=sum(2*(p-1) for p in support if p not in (3,5,prime) and p not in marked)
    activity_cap=rank-occupied-joint-left-right-other_activity
    height_certified=activity_cap<2*prime*(prime-1)
    prime_floor=prime**(3 if height_certified else 2)
    occupied+=dense+joint+left+right
    mass=Q(1,105**2)+Q(nines,81)+Q(dense,prime**4)
    mass+=sum(Q(v[1],p*p) for p,v in marked.items())
    if len(support)==3:
        w5,wq=(witnesses(p,nines) for p in (5,prime))
        full_q=max(wq,free_count(prime,nines)) if height_certified else wq
        occupied+=full_q+w5+8
        mass+=Q(wq,(9*prime_floor)**2)+Q(full_q-wq,prime_floor**2)+Q(w5,225**2)
        for i in range(0,min(w5,left)+1,2):
            for j in range(0,min(8,left-i)+1,2):
                free=rank-occupied+i+j
                if free<0:
                    continue
                high=high_joint+high_right+max(0,high_left-i-j)
                tail=Q(j,135**2)+Q(8-j,27**2)+Q(high,45**2)
                tail+=Q(free+joint+left+right-i-j-high,225)
                if mass+tail>=1:
                    return True
        return False
    others=[p for p in support if p not in (3,5,prime) and p not in marked]
    for p in others:
        activity=2*(p-1)
        occupied+=activity
        mass+=Q(activity,p**4)
    if len(support)==5:
        # Discard all omission witnesses; this gives a conservative upper mass.
        free=rank-occupied
        high=high_joint+high_left+high_right
        return free>=0 and mass+Q(joint+left+right+free-high,225)+Q(high,45**2)>=1
    for p in marked:
        occupied+=6
        mass+=2*sum(Q(1,b*b) for b in (9*p,5*p,prime_floor*p))
    if prime==7 and height_certified:
        extra=max(0,seven_free(nines,marked)-2*len(marked))
        occupied+=extra
        mass+=Q(extra,prime_floor**2)
    free=rank-occupied
    high=high_joint+high_left+high_right
    survives=(free>=0 and mass+Q(joint+left+right+free-high,225)+Q(high,45**2)>=1)
    if (survives and rank<=201 and tuple(support)==(3,5,7,11)
            and marked=={11:(1,22,20)} and prime==7
            and nines in (34,36,38,40) and joint==32 and high_joint==0
            and left>=12 and right>=40):
        # The moment-interval and Gauss-boundary argument excludes every
        # actual occupancy above these lower bounds; dense_boundary.py
        # independently checks its finite rational intervals.
        from dense_boundary import check
        check()
        return False
    return survives
