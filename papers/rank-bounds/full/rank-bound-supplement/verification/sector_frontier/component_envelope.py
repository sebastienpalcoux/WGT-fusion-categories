"""Joint suffix reservations in disconnected twist components."""
from functools import lru_cache
from itertools import combinations
from math import isqrt,prod
from neutral_cap import cap

def squares(n):
    if isqrt(n)**2==n:return 1
    if any(isqrt(n-a*a)**2==n-a*a for a in range(1,isqrt(n)+1)):return 2
    while n%4==0:n//=4
    return 4 if n%8==7 else 3

@lru_cache(maxsize=None)
def compplans(S,F,q,C,Bq,low,least_multiple):
    P=sorted(set().union(*(part for part,k in F)))
    comps=[]
    for p in P:
        if any(p in c for c in comps):continue
        c={p}
        while True:
            new=c.union(*(part for part,k in F if c&part))
            if new==c:break
            c=new
        comps.append(c)
    # A row has (quota, eligible divisibility, orbit generators, future floor,
    # pure prime). A group of rows has disjoint booked witnesses. Different
    # groups are alternative estimates on the same component's population.
    plans=[]
    for comp in comps:
        alternatives=[]
        facet=[]
        isolated=next(iter(comp)) if len(comp)==1 and q==1 and 2 not in comp else 0
        def policy(v,k,h,pure=0):
            if pure:
                floor=pure
                while floor<low:floor*=pure
            else:floor=least_multiple(S,low,k)
            return (v,k,tuple(h),floor,pure)
        for face,k in F:
            if not(comp&face):continue
            for n in range(1,len(face)+1):
                for K in combinations(face,n):
                    size=prod(t-1 for t in K)
                    h=prod((p-1)//2 for p in K if p!=2)
                    facet.append(policy(C*size,k,[h],isolated))
            for p,t in ((p,t) for p in face if p!=2 for t in face if p!=t):
                h=4 if t==2 else 3 if t%4==3 else 2
                alternatives.append([policy(h*(p-1),p*t,[(p-1)//2])])
                if p>=5 and t>2:
                    v=(p-1+squares(p+1))*(t-1)+p-1
                    alternatives.append([policy(v,p*t,[(p-1)//2,(t-1)//2])])
        alternatives.append(facet)
        alternatives.append([policy(Bq*(p-1),p,[(p-1)//2],isolated) for p in comp if p!=2])
        if isolated:alternatives.append([policy(isolated,isolated,[1],isolated)])
        plans.append(alternatives)
    levels=sorted({row[3] for comp in plans for group in comp for row in group})
    return tuple((d,tuple(tuple(tuple(row for row in group if row[3]>=d) for group in comp) for comp in plans)) for d in levels)

def envelope(S,F,q,hist,low,C,Bq,base,Mfloor,least_multiple):
    h=len(base)
    for d,comps in compplans(S,F,q,C,Bq,low,least_multiple):
        missing=0
        for comp in comps:
            worst=0
            for group in comp:
                if not group:continue
                target=sum(row[0] for row in group)
                used=0
                for b,n in hist.items():
                    sizes=set()
                    for count,k,parts,floor,pure in group:
                        if b%k or pure and b!=pure**_v(b,pure):continue
                        sizes.update(parts)
                    used+=cap(n,tuple(sorted(sizes)))
                worst=max(worst,target-used)
            missing+=worst
        if missing>=h:return None
        row=sorted([Mfloor]+[d]*missing+[low]*(h-1-missing))
        base=list(map(max,base,row))
    return base

@lru_cache(maxsize=None)
def _v(n,p):
    e=0
    while n%p==0:e+=1;n//=p
    return e
