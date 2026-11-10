"""Independent signed reciprocal enumeration; close the fourth term by lcm."""
from fractions import Fraction
from itertools import combinations, product
from math import gcd,isqrt,lcm

def divisors(n):
    for i in range(1,isqrt(n)+1):
        if n%i==0:
            yield i
            if i*i!=n:
                yield n//i

def verify():
    out=set()
    def search(ks,r):
        ell=4-len(ks)
        if r==0:return
        if ell==1:
            L=lcm(*ks)
            C=int(L*r)
            for g in divisors(L):
                for e,delta in product((-1,1),repeat=2):
                    if (e*L+delta*g)%C:continue
                    k=(e*L+delta*g)//C
                    if k<ks[-1] or k%3==0 or gcd(L,k)!=g:continue
                    K=(*ks,k);n=lcm(*K);ds=[n//j for j in K]
                    if len(set(K))<3 or min(ds)<6 or gcd(*ds)!=1:continue
                    if all(gcd(x,y)>1 for x,y in combinations(ds,2)):
                        assert abs(r-Fraction(e,k))==Fraction(1,n)
                        out.add(n)
            return
        hi=int((ell+Fraction(1,6))/abs(r))
        for k in range(ks[-1] if ks else 2,hi+1):
            if k%3:
                for e in [-1,1]:search((*ks,k),r-Fraction(e,k))
    search((),Fraction(1))
    return sorted(out)

if __name__=='__main__':print(verify())
