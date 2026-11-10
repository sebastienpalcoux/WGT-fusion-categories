"""An independent short certificate for the three-prime square class ten."""
from fractions import Fraction
from math import isqrt, lcm, prod


def exists(u,rank,cut=50,dyadic=8):
    P=(2,5,u)
    a=(dyadic,25,u**(3 if u==3 else 2))
    small=set()
    for b in range(1,cut+1):
        z=b
        for p in P:
            while z%p==0: z//=p
        if z==1: small.add(b)
    values=sorted(small)
    def rec(i,C,S,count,last):
        if S>=10: return False
        L=lcm(*C)
        height=[]
        for p,v in zip(P,a):
            while L%(p*v)==0: v*=p
            height.append(v)
        # All p-free witnesses form disjoint sets, with the three cardinalities below.
        required=[]
        for v,n in zip(height,(3,3,2)):
            missing=max(0,n-sum(m for d,m in C.items() if d%v==0))
            required.extend([v]*missing)
        rem=rank-1-count
        if len(required)>rem: return False
        required=sorted([1]*(rem-len(required))+required)
        missing=max(0,u-1-sum(m for d,m in C.items() if d%u==0))
        if missing>rem: return False
        for j in range(rem-missing,rem): required[j]=max(required[j],u)
        upper=S+Fraction(1,prod(height)**2)
        for d in required:
            # smallest P-smooth integer above both d and the processed prefix
            t=max(d,last+1)
            while True:
                z=t
                for p in P:
                    while z%p==0: z//=p
                if z==1: break
                t+=1
            upper+=Fraction(1,t*t)
        if upper<10:return False
        # Formal defect one at 5 consists of sets of five equals and two-plus-two.
        for c in range(1,last+1,2):
            if c%5==0:continue
            b=c; debt=0
            while b<=last:
                n=C.get(b,0)-debt
                # 2*k+5*s=n, debt at the next b is 2*k. Choose least k.
                options=range(0,n//2+1)
                k=next((k for k in options if (n-2*k)%5==0),None)
                if k is None:return False
                debt=2*k; b*=2
        if i==len(values):return True
        d=values[i]
        maxi=min(rem,7 if d==1 else rem)
        for m in range(maxi,-1,-1):
            C2=C.copy()
            if m:C2[d]=m
            if rec(i+1,C2,S+Fraction(m,d*d),count+m,d):return True
        return False
    return rec(0,{},Fraction(0),0,0)

if __name__=='__main__':
    for kind, cutoff, dyadic in [('ordinary',35,8),('reduced',33,4)]:
        P=[u for u in range(3,cutoff+1,2) if u!=5 and all(u%k for k in range(3,isqrt(u)+1,2))]
        admitted=[u for u in P if exists(u,cutoff,dyadic=dyadic)]
        print(kind, 'rank at most',cutoff, 'admitted primes',admitted)
        if admitted: raise SystemExit(1)
