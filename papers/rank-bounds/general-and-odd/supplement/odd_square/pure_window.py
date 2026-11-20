"""Charge intervals for the unassigned pure-prime activity at a search node."""
from functools import lru_cache
from itertools import product
from math import lcm

@lru_cache(maxsize=20000)
def power_brackets(levels, free, need, cap, isolated):
    # Common scale for the first two powers. Each prime at which a label
    # can occur is a divisor of the dimension. No upper bound is imposed.
    primes=(3,5)+tuple(isolated)
    unit=4*lcm(*(p*b for p,b in zip(primes,levels)))**2
    rate={p:p*(unit//b**2) for p,b in zip(primes,levels)}
    next_rates={p:(t//(p*p)+p-2)//(p-1) for p,t in rate.items()}
    slots=free+sum(need)
    # First-level numbers may be zero. Unspent slots beyond it contribute
    # anywhere between zero and the stated next-power upper bound.
    intervals=[]
    options=[range(0,c+1,2*(p-1)) for p,c in zip(isolated,cap)]
    for first in product(*options):
        leftover=tuple(max(0,n-t) for n,t in zip(need,first))
        space=slots-sum(first)-sum(leftover)
        if space<0:continue
        qcharge=sum(rate[p]*(t//(p-1)) for p,t in zip(isolated,first))
        booked=sum(next_rates[p]*n for p,n in zip(isolated,leftover))
        best=max([next_rates[3],next_rates[5]]+[next_rates[p] for p,c,t,n in zip(isolated,cap,first,leftover) if c>t+n])
        for i in range(space//2+1):
            for j in range((space-2*i)//4+1):
                low=qcharge+i*rate[3]+j*rate[5]
                intervals.append((low,low+booked+(space-2*i-4*j)*best))
    merged=[]
    for a,b in sorted(intervals):
        if merged and merged[-1][1]>=a:
            merged[-1][1]=max(b,merged[-1][1])
        else:
            merged.append([a,b])
    return unit,tuple(b for a,b in merged),tuple(a for a,b in merged)
