"""Ordinary mass and pair residues at a terminal signed Gauss root."""
from full_blocks import neutral_pair_floors
from functools import lru_cache


@lru_cache(None)
def roots(M, support):
    remainder=M
    parts={}
    divisors=[1]
    for p in support:
        power=1
        new=list(divisors)
        while remainder%p==0:
            remainder//=p
            power*=p
            new.extend(power*d for d in divisors)
        if power==1:
            return None
        parts[p]=power
        divisors=new
    if remainder!=1:
        return None
    allowed=tuple(sorted((b for b in divisors if b>=15 and
                         sum(M//b%p==0 for p in support)>=2),reverse=True))
    return parts,allowed


def possible(M,support,marked,L,nines,k,neutral,o3,o5,mass,scale,witnesses,pure_floors,exact_pure):
    data=roots(M,support)
    if data is None:
        return False
    parts,allowed=data
    if any(parts[p]!=p for p in marked):
        return False
    for p in support:
        if p in (3,5) or p in marked:
            continue
        if L%(parts[p])!=0:
            return False
    # The charged classes and the vacuum have been removed. All remaining
    # labels occur in dual pairs, have composite degree, and have denominator >=15.
    num,den=divmod(M*M*mass,scale)
    if den:
        return False
    deficit=M*M-1-num
    if deficit<0 or deficit%2:
        return False
    target=deficit//2
    isolated=[p for p in support if p not in (3,5) and p not in marked]
    remaining_count=k+neutral-2*sum(exact_pure)
    neutral-=2*len(isolated)
    if len(support)==3:
        q=support[2]
        hits=witnesses(q,nines)
        floors=([9*parts[q]]*(hits//2)+[parts[q]]*((neutral-hits)//2)
                +[max(15,parts[3])]*(o3//2)+[9*parts[5]]*(o5//2))
    else:
        floors=neutral_pair_floors(support,marked,parts,nines)
    floors+=[floor for floor,exact in zip(pure_floors,exact_pure) if not exact]
    if len(floors)>remaining_count//2:
        return False
    # Solve the last reciprocal-square interval modulo eight and at the
    # largest few allowed degrees. Unresolved small terms are relaxed to
    # a closed interval, so that absence is a valid exclusion.
    for n in range(target%8,remaining_count//2+1,8):
        if n<len(floors):
            continue
        limits=sorted(floors+[15]*(n-len(floors)),reverse=True)
        # A larger denominator costs less ordinary mass. Visiting the
        # largest degrees first gives a sharp upper bound after a skip.
        degrees=sorted({(M//b)**2 for b in allowed if (M//b)**2<=target},reverse=True)
        ceilings=[(M//floor)**2 for floor in reversed(limits)]
        def interval(i,left,slots,step):
            if not slots:
                return left==0
            if i>=len(degrees) or left<slots or left>sum(min(degrees[i],v) for v in ceilings[n-slots:]):
                return False
            if step==0:
                return True
            value=degrees[i]
            cap=min(slots,left//value)
            cap=min(cap,sum(v>=value for v in ceilings[n-slots:]))
            return any(interval(i+1,left-count*value,slots-count,step-1) for count in range(cap,-1,-1))
        if interval(0,target,n,6):
            return True
    return False
