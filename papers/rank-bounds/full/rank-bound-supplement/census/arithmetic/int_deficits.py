"""Exact integer implementation of phase_deficits.sparse_budget.

Only arithmetic representation and counters differ from the reference recursion.
All bounds are rounded down only when compared to an integer denominator.
A deadline raises TimeoutError and never certifies an unfinished profile.
"""
from math import gcd, isqrt, lcm
from itertools import product
from fractions import Fraction as F
from time import monotonic

class Stats:
    def __init__(self, seconds=None):
        self.start=monotonic(); self.deadline=None if seconds is None else self.start+seconds
        self.nodes=0; self.yields=0; self.max_depth=0; self.starts=0
    def tick(self, depth):
        self.nodes+=1
        if depth>self.max_depth:self.max_depth=depth
        if self.nodes%4096==0 and self.deadline is not None and monotonic()>self.deadline:
            raise TimeoutError('structural profile deadline')
    def info(self):
        return dict(nodes=self.nodes, yields=self.yields,max_depth=self.max_depth,
                    starts=self.starts,seconds=monotonic()-self.start)

def prime_free(n,p):
    while n%p==0:n//=p
    return n

def sparse_budget(residual,multiple,odd_primes,budget,stats,max_eights=None,protected=()):
    primes=(2,)+tuple(odd_primes)
    allmask=(1<<len(odd_primes))-1
    missing_cost=tuple(sum(p-1 for i,p in enumerate(odd_primes) if mask>>i&1)
                       for mask in range(1<<len(odd_primes)))
    missing_adjust=tuple(sum(p-2 for i,p in enumerate(odd_primes) if mask>>i&1)
                         for mask in range(1<<len(odd_primes)))
    bits={p:1<<i for i,p in enumerate(odd_primes)}
    def visit(a,b,required,least,left,missing,eights,depth):
        stats.tick(depth)
        if a==0:
            assert protected
            return
        aa=abs(a)
        if not missing and aa==1 and b%required==0:
            stats.yields+=1
            yield a,b
        max_terms=left-missing_adjust[missing]
        if max_terms<=0 or aa*least*required>b*(max_terms*required+least):return
        if a<0 and aa*required>=b:return
        for p in primes:
            cost=max(1,p-1)
            rest_missing=missing&~bits.get(p,0)
            if left-cost<missing_cost[rest_missing]:continue
            coprime=prime_free(required,p)
            ceiling=max((2*max_terms*b)//aa,(4*b*b)//((coprime*aa)**2))
            if aa*required>b:
                ceiling=min(ceiling,(max_terms*b*required)//(aa*required-b))
            d=32 if p==2 else p**3
            factor=p if p==2 else p*p
            while d<least:d*=factor
            while d<=ceiling:
                underlying=isqrt(d)
                if underlying*underlying!=d:underlying=isqrt(p*d)
                k=eights+int(p==2 and underlying==8)
                if max_eights is None or k<=max_eights:
                    next_least=d+1 if p==2 and d>32 and max_eights is None else d
                    na=a*d-b;nb=b*d;g=gcd(na,nb)
                    yield from visit(na//g,nb//g,lcm(required,d if p in protected else underlying),
                                     next_least,left-cost,rest_missing,k,depth+1)
                d*=factor
    yield from visit(residual.numerator,residual.denominator,multiple,1,budget,allmask,0,0)

def pure_check(rank,rows,support,modes,neutral,stats,progress=None):
    prime=[p for p,z in zip(support,modes) if z==2]
    odd=tuple(p for p in support if p not in prime)
    M0=8
    for p,z in zip(support,modes):M0*=p if z==2 else p**(3 if p>=5 else 2)
    out=set()
    for js in product(*((1,) if p==5 else range(1,p) for p in prime)):
        ms=tuple(p if j==1 else (p-1)*j for p,j in zip(prime,js))
        budget=rank-1-rows-sum(ms)-neutral
        if budget<sum(p-1 for p in odd):continue
        rem=1-F(rows,16)-sum(F(n,p*p) for n,p in zip(ms,prime))
        protected={p for p in odd if budget-sum(t-1 for t in odd if t!=p)<p*(p-1)}
        for a in range(-rows,rows+1):
            val=1-F(rows-a,16)-sum(F(j,p) for j,p in zip(js,prime))
            if val+F(1,M0)<=0 or val>2*rem+F(1,M0):continue
            stats.starts+=1
            if progress:progress(js,a,stats.info())
            for sign,M in sparse_budget(val,M0,odd,budget,stats,protected=protected):
                if M%M0:continue
                rest=M
                for p in (2,)+tuple(support):rest=prime_free(rest,p)
                if rest==1:out.add(M)
    return tuple(sorted(out))
