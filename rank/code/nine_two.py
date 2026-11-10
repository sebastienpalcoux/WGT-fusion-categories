from fractions import Fraction
from math import gcd,isqrt
import itertools

def solutions():
    def rec(prefix, residual):
        left=4-len(prefix)
        if not left:
            if abs(residual.numerator)!=1:return
            n=residual.denominator
            K=[k for k,s in prefix]
            if n%3==0 or any(n%k or n//k<6 for k in K) or len(set(K))<3:return
            degrees=[n//k for k in K]
            if gcd(*degrees)!=1 or any(gcd(u,v)==1 for u,v in itertools.combinations(degrees,2)):return
            yield n,prefix
            return
        if not residual:return
        low=prefix[-1][0] if prefix else 2
        cap=(Fraction(left)+Fraction(1,6))/abs(residual)
        for k in range(low,int(cap)+1):
            if k%3==0:continue
            for sign in [1,-1]:
                yield from rec(prefix+[(k,sign)], residual-Fraction(sign,k))
    return rec([],Fraction(1))
if __name__=='__main__':
 from collections import Counter
 s=list(solutions()); print(len(s),Counter(n for n,kk in s));print(s[:20]);
