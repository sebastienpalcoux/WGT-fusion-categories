"""Enumerate ordinary or reduced perfect simple dimension types using exact integers."""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import groupby
from math import gcd, isqrt, lcm
import json
import sys
import time

def primes_through(bound):
    return [n for n in range(2, bound + 1)
            if all(n % j for j in range(2, isqrt(n) + 1))]

def is_smooth(n, allowed):
    for p in allowed:
        while n % p == 0:
            n //= p
    return n == 1

def smooth_interval(lo, hi, allowed):
    if hi - lo < 96:
        return (b for b in range(lo, hi+1) if is_smooth(b, allowed))
    # All products in the interval, each once, using successive largest primes.
    def visit(index, n):
        if index == 0:
            b = n
            if b < lo:
                b <<= max(0, ((lo - 1) // b).bit_length())
            while b <= hi:
                yield b
                b *= 2
            return
        p = allowed[index]
        while n <= hi:
            yield from visit(index - 1, n)
            n *= p
    return visit(len(allowed) - 1, 1)

def enumerate_types(rank_cap, mode, min_codegree, with_orbits=True):
    # Odd prime twist orbits in the positive half-table are preserved up to fermion.
    primes = primes_through(2 * rank_cap + 1)
    orbit_sieve = with_orbits
    bands = [sum(1 << i for i,p in enumerate(primes) if (p-1)//2 <= n)
             for n in range(rank_cap+1)]
    @lru_cache(maxsize=131072)
    def factor_mask(n):
        return sum(1 << i for i,p in enumerate(primes) if n % p == 0)
    stats = Counter()
    per_rank = {}
    accepted = []
    def count(rank,stage):
        per_rank.setdefault(rank, Counter())[stage] += 1
    @lru_cache(maxsize=16384)
    def divisors(n):
        # Ordinary branches are smooth. Reduced search uses unrestricted integers.
        result = [1]
        trial = 2
        while trial*trial <= n:
            power = 1
            old = list(result)
            while n % trial == 0:
                n //= trial
                power *= trial
                result.extend(d*power for d in old)
            trial += 1 if trial == 2 else 2
        if n > 1:
            result.extend(d*n for d in list(result))
        return tuple(result)
    def eligible_orbits(q, prefix, slots):
        # For a known prime its equal-denominator active run can be in the prefix,
        # extend the last run, or occupy entirely new positions.
        if not with_orbits:
            return True
        runs = [(a, len(list(group))) for a,group in groupby(prefix)]
        for p in primes:
            h = (p-1)//2
            if h <= slots or (q % p and all(a % p for a,_ in runs)):
                continue
            if not any(q*a % p == 0 and n + (slots if i+1==len(runs) else 0) >= h
                       for i,(a,n) in enumerate(runs)):
                return False
        return True
    def candidate(q, M, denominators):
        rank=len(denominators)+1
        degrees=tuple(M//a for a in denominators)
        count(rank,"lcm_vacuum")
        # A number with two prime factors has a nontrivial unitary divisor.
        if any(not any(1<p<d and gcd(p,d//p)==1 for p in divisors(d)) for d in set(degrees)):
            return
        count(rank,"mixed_degrees")
        if any(gcd(*[y for i,y in enumerate(degrees) if i!=j])!=1 for j in range(len(degrees))):
            return
        count(rank,"double_avoidance")
        if not eligible_orbits(q,denominators,0):
            return
        count(rank,"orbit_multiplicities")
        D=q*M*M
        if mode=="reduced" and D % 2:
            return
        degree_set=set(degrees)
        for x in degree_set:
            possible=[y for y in degree_set if gcd(x,y)>1 and lcm(x,y)**2+x*x<=D]
            if gcd(*possible,0)!=1:
                return
        count(rank,"row_gcd")
        accepted.append([q,M,list(degrees)])
    def visit(q, qmask, L, E, floor, prefix, seen, witnessed, run, active):
        # E/L² is the remaining mass. All completed denominators divide L.
        slots=rank_cap-len(prefix)-1
        stats["visited"]+=1
        forced = seen & ~witnessed & ~bands[slots] if orbit_sieve else 0
        if forced & ~(active & bands[slots+run]):
            stats["orbit_pruned"]+=1
            return
        if slots<=0:
            return
        low=max(floor,isqrt(L*L//E)+1)
        # Deleting any nonunit leaves lcm M: no prime has a unique p-free
        # nonunit label. In particular, before the last denominator M=L.
        d = isqrt(E-1)
        if d >= 6 and d*d+1==E and L%d==0:
            b=L//d
            if b>=floor:
                candidate(q,L,prefix+(b,))
        if slots<2:
            return
        high=isqrt((36*slots+1)*L*L//(36*E))
        if low>high:
            return
        if slots==2 and high-low > len(divisors(L))**2//4:
            # Both tail degrees divide L. Each prime in k=M/L is omitted
            # by both, since all prefix degrees are multiples of k.
            tail_degrees=[d for d in divisors(L) if d>=6 and
                          any(1<t<d and gcd(t,d//t)==1 for t in divisors(d))]
            for i,one in enumerate(tail_degrees):
                for other in tail_degrees[i:]:
                    d,e=max(one,other),min(one,other)
                    k2,rem=divmod(1+d*d+e*e,E)
                    if rem: continue
                    k=isqrt(k2)
                    if k*k!=k2 or gcd(k,d*e)!=1: continue
                    M=L*k
                    b,c=M//d,M//e
                    if b<floor or lcm(L,b,c)!=M: continue
                    if orbit_sieve and not is_smooth(k,primes):continue
                    candidate(q,M,prefix+(b,c))
            return
        if forced:
            choices=(floor,) if low<=floor<=high else ()
        elif orbit_sieve:
            allowed=[p for p in primes if p<=2*slots+1 or q*L%p==0]
            choices=smooth_interval(low, high, allowed)
        else:
            choices=range(low,high+1)
        for b in choices:
            g=gcd(L,b)
            k=b//g
            run_b = run+1 if b==floor and prefix else 1
            bmask = factor_mask(b) if orbit_sieve else 0
            active_b = bmask|qmask
            visit(q,qmask,L*k,E*k*k-(L//g)**2,b,prefix+(b,),
                  seen|active_b,witnessed | (active_b & bands[run_b]),run_b,active_b)
    for q in range(1,rank_cap):
        if any(q%(p*p)==0 for p in primes):
            continue
        floor=isqrt((min_codegree-1)//q)+1
        t=time.process_time()
        qmask=factor_mask(q) if orbit_sieve else 0
        visit(q,qmask,1,q,floor,(),qmask,0,0,0)
        print(f"q={q} CPU={time.process_time()-t:.3f}, nodes={stats['visited']}",file=sys.stderr,flush=True)
    accepted.sort()
    return {"mode":mode,"max_rank":rank_cap,"min_codegree":min_codegree,
            "nodes":dict(stats),"rank_counts":dict(sorted(per_rank.items())),
            "row_feasible":accepted}
if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("rank",type=int)
    parser.add_argument("output")
    parser.add_argument("--mode",choices=["ordinary","reduced"],default="reduced")
    parser.add_argument("--min-codegree",type=int,default=3)
    parser.add_argument("--no-orbits",action="store_true")
    args=parser.parse_args()
    data=enumerate_types(args.rank,args.mode,args.min_codegree,not args.no_orbits)
    with open(args.output,"w") as stream:
        json.dump(data,stream,separators=(",",":"))
        stream.write("\n")
