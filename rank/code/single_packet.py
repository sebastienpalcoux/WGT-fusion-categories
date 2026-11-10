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
    # The orbit bounds are used only in the ordinary modular case.
    primes = primes_through(2 * rank_cap + 1)
    orbit_sieve = with_orbits and mode == 'ordinary'
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
        if not with_orbits or mode != 'ordinary':
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
        if not eligible_orbits(q,denominators,0):
            return
        count(rank,"orbit_multiplicities")
        D=q*M*M
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
        # One last nonunit degree d=L/g, where g=gcd(b,L).
        # With b=k*g, the residual equation reads E*k²=d²+1.
        last_hi=isqrt(37*L*L//(36*E))
        # A prime whose required orbit no longer fits in all remaining slots
        # forces the current run of equal denominators to continue.
        if forced and low>floor:
            return
        if low<=last_hi:
            if forced:
                last_hi = min(last_hi, floor)
            if last_hi-low+1 < len(divisors(L)):
                for b in range(low,last_hi+1):
                    g=gcd(b,L)
                    M=L*(b//g)
                    if L//g>=6 and E*(b//g)**2==(L//g)**2+1:
                        if mode=='ordinary' and with_orbits and not is_smooth(b,primes):
                            continue
                        candidate(q,M,prefix+(b,))
            else:
                for d in divisors(L):
                    if d<6: continue
                    sq, rem=divmod(d*d+1,E)
                    if rem: continue
                    k=isqrt(sq)
                    if k*k!=sq or gcd(k,d)!=1: continue
                    b=k*(L//d)
                    if b>=floor:
                        if mode=='ordinary' and with_orbits and not is_smooth(b,primes):
                            continue
                        candidate(q,L*k,prefix+(b,))
        if slots<2:
            return
        high=isqrt((36*slots+1)*L*L//(36*E))
        if low>high:
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
    parser.add_argument("--mode",choices=["ordinary","reduced"],default="ordinary")
    parser.add_argument("--min-codegree",type=int,default=5)
    parser.add_argument("--no-orbits",action="store_true")
    args=parser.parse_args()
    data=enumerate_types(args.rank,args.mode,args.min_codegree,not args.no_orbits)
    with open(args.output,"w") as stream:
        json.dump(data,stream,separators=(",",":"))
        stream.write("\n")
