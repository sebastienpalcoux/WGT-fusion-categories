"""Fresh reciprocal-square enumeration, with exact closure of the last two terms."""
import argparse
import json
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm
from pathlib import Path
from families import dimension_family, factors, family


@lru_cache(maxsize=32768)
def mixed_divisors(n):
    items = [(1,0)]
    for p,e in factors(n):
        extra = []
        power = 1
        for _ in range(e):
            power *= p
            extra.extend((d*power,w+1) for d,w in items)
        items.extend(extra)
    return tuple(sorted(d for d,w in items if w >= 2))


def smooth(lo, hi, ps):
    if lo > hi:
        return
    if hi - lo < 32:
        for n in range(lo, hi+1):
            if all(p in ps for p,_ in factors(n)):
                yield n
        return
    def visit(i, n):
        if i == 0:
            if n < lo:
                n <<= ((lo-1)//n).bit_length()
            while n <= hi:
                yield n
                n *= 2
            return
        while n <= hi:
            yield from visit(i-1,n)
            n *= ps[i]
    yield from visit(len(ps)-1,1)


def enumerate_rank(R, min_f=5, use_families=True, strong_packets=False, q_only=None):
    """Ordinary simple perfect modular types of rank at most R."""
    primes = tuple(p for p in range(2,2*R+2) if factors(p)==((p,1),))
    bands = [sum(1<<i for i,p in enumerate(primes) if (p-1)//2 <= m) for m in range(R+1)]
    @lru_cache(maxsize=100000)
    def mask(n):
        return sum(1<<i for i,p in enumerate(primes) if n%p==0)
    @lru_cache(maxsize=32768)
    def allowed(q,L,h):
        limit = h+1 if strong_packets else 2*h+1
        return tuple(p for p in primes if p==3 or p<=limit or (q*L)%p==0)
    stats = Counter()
    accepted = []

    def emit(q,M,den):
        D=q*M*M
        if use_families and dimension_family(D):
            return
        ds=tuple(sorted(M//b for b in den))
        if any(M%b or len(factors(M//b))<2 for b in set(den)):
            return
        if any(sum(d%p!=0 for d in ds)<2 for p,_ in factors(M)):
            return
        mult=Counter(den)
        for p,_ in factors(D):
            h=(p-1)//2
            if p>=5 and not any((q*b)%p==0 and c>=h for b,c in mult.items()):
                return
        for x in set(ds):
            hit=[y for y in set(ds) if gcd(x,y)>1 and x*x+lcm(x,y)**2<=D]
            if gcd(*hit,0)!=1:
                return
        accepted.append((1,)+ds)

    def last_two(q,L,E,lower,den):
        # The two unidentified degrees are divisors of L, prime to k=M/L.
        # A new prime in k can only be 3. Two missing labels cannot
        # introduce 2 modulo 4 or 5 by its same-degree orbit modulo 5.
        if E%4==0:
            return
        base=dict(factors(L))
        possible=sorted(({3} | set(base) | {p for p,_ in factors(q)})-{2})
        cap=isqrt((2*L*L+1)//E)
        caps=dict(factors(q*L*L))
        for p in possible:
            u=L//p**base.get(p,0)
            upper=isqrt((2*u*u+1)//E)
            v=0
            while upper>=p:
                upper//=p
                v+=1
            caps[p]=caps.get(p,0)+2*v
        if use_families and family(sorted(caps.items())):
            stats['whole_tail_family']+=1
            return
        def scales(i,k,u):
            if i==len(possible):
                yield k,u
                return
            p=possible[i]
            yield from scales(i+1,k,u)
            u//=p**base.get(p,0)
            upper=min(cap,isqrt((2*u*u+1)//E))
            k*=p
            while k<=upper:
                yield from scales(i+1,k,u)
                k*=p
        for k,u in scales(0,1,L):
            M=L*k
            if use_families and dimension_family(q*M*M):
                continue
            target=E*k*k-1
            if target < 72:
                continue
            lo=isqrt((target-1)//2)+1
            hi=min(isqrt(target-36),M//lower,u)
            if lo>hi:
                continue
            for d in mixed_divisors(u):
                if d<lo: continue
                if d>hi: break
                e=isqrt(target-d*d)
                if e>=6 and e*e==target-d*d and u%e==0:
                    emit(q,M,den+(M//d,M//e))

    def visit(q,L,E,lower,den,seen,done,run,active):
        h=R-1-len(den)
        stats['prefixes']+=1
        forced=seen & ~done & ~bands[h]
        if forced & ~(active & bands[h+run]):
            return
        d=isqrt(E-1)
        if d>=6 and d*d+1==E and L%d==0 and L//d>=lower:
            emit(q,L,den+(L//d,))
        if h==2:
            last_two(q,L,E,lower,den)
            return
        lo=max(lower,isqrt(L*L//E)+1)
        hi=isqrt((36*h+1)*L*L//(36*E))
        if forced:
            choices=(lower,) if lo<=lower<=hi else ()
        else:
            choices=smooth(lo,hi,allowed(q,L,h))
        for b in choices:
            g=gcd(b,L)
            k=b//g
            length=run+1 if den and b==lower else 1
            bearing=mask(q)|mask(b)
            visit(q,L*k,E*k*k-(L//g)**2,b,den+(b,),seen|bearing,
                  done|(bearing&bands[length]),length,bearing)
    for q in range(1,R):
        if q_only is not None and q!=q_only: continue
        if any(e>1 for _,e in factors(q)): continue
        lower=isqrt((min_f-1)//q)+1
        if lower==1 and factors(q)==((q,1),): lower=2
        visit(q,1,q,lower,(),mask(q),0,0,0)
    rows=sorted(set(accepted),key=lambda d:(len(d),d))
    assert len(rows)==len(accepted)
    return {'rank_cap':R,'minimum_codegree':min_f,'strong_packets':strong_packets,
            'families':use_families,'q_only':q_only,'counts':dict(stats),'types':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('rank',type=int)
    p.add_argument('output',type=Path)
    p.add_argument('--minimum',type=int,default=5)
    p.add_argument('--without-families',action='store_true')
    p.add_argument('--strong-packets',action='store_true')
    p.add_argument('--q',type=int,help='one square-free kernel (a shard only)')
    args=p.parse_args()
    if args.rank<4:p.error('rank cap must be at least 4')
    out=enumerate_rank(args.rank,args.minimum,not args.without_families,args.strong_packets,args.q)
    args.output.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['counts'],len(out['types']))
