"""Necessary rank tests for an ordinary fusion-simple perfect integral modular category."""
from bisect import bisect_left
from functools import lru_cache
from itertools import combinations
from math import isqrt, lcm, prod


def is_prime(p):
    return p >= 2 and all(p % i for i in range(2, isqrt(p)+1))


@lru_cache(None)
def witness_count(p, other_residues, odd, rank):
    """Minimum number of nonunit labels, or rank+1 if no allowed count fits."""
    modulus = p*p
    step = 2 if odd else 1
    generators = {res*res % modulus for res in other_residues}
    H = {1}
    queue = [1]
    while queue:
        x = queue.pop()
        for g in generators:
            y = x*g % modulus
            if y not in H:
                H.add(y)
                queue.append(y)
    full = (1 << modulus)-1
    allowed = tuple(step*h % modulus for h in H)
    reachable = 1
    for j in range(1, (rank-1)//step + 1):
        nxt = 0
        for a in allowed:
            nxt |= ((reachable << a) | (reachable >> (modulus-a))) & full
        if j >= 2 and nxt >> (modulus-1) & 1:
            return step*j
        if nxt == reachable and not (nxt >> (modulus-1) & 1):
            break
        reachable = nxt
    return rank+1


@lru_cache(None)
def divisors_only(primes, ceiling):
    ds = [1]
    for p in primes:
        added = []
        for d in ds:
            x = d*p
            while x <= ceiling:
                added.append(x)
                x *= p
        ds.extend(added)
    return sorted(ds)


def multiplicity_cap(f, largest, prime_support, q, local_cap):
    if f == 8:
        if largest and len(prime_support)<4: return 0
        return 2 if largest else 4
    if f == 9:
        if largest and len(prime_support)<4: return 0
        return 4
    if f < 8 or is_prime(f):
        return 0
    cap = f-1
    if local_cap and any(f%p for p in prime_support):
        for p in prime_support:
            if p%2 and f%p == 0 and f%(p*p):
                cap = min(cap, f-1-(p-1)//2)
    return cap


def optimistic_mass(primes, heights, q, rank, B, local_cap=True, prefix=(), after=0):
    """Upper bound the allowable squared-Egyptian mass for a whole valuation box."""
    from fractions import Fraction
    k = len(primes)
    odd = 2 not in primes
    slots = rank-1-len(prefix)
    if odd and rank%8 != q%8:
        return Fraction(0)
    # [primes, heights] describe a *lower* M-box with fixed squareclass,
    # so the prime-power witness bounds can also cover every larger exponent.
    if prefix:
        current=lcm(*prefix)
        updated=[]
        for p,a in zip(primes,heights):
            while current%(p**(a+1))==0:
                a+=1
            updated.append(a)
        heights=updated
    M0 = prod(p**a for p,a in zip(primes, heights))
    if q*B*B >= rank:
        return Fraction(0)
    if not multiplicity_cap(q*B*B, True, primes, q, local_cap):
        return Fraction(0)
    witnesses = []
    tails = {}
    for p, a in zip(primes, heights):
        other = tuple(sorted({i % (p*p) for i in primes if i != p}))
        w = witness_count(p, other, odd, rank)
        w=max(0,w-sum(b % p**a == 0 for b in prefix))
        if w > slots:
            return Fraction(0)
        witnesses.append((p**a,w))
        if p>2 and q%p:
            tails[p] = max(0,(2 if odd else 1)*(p-1)-sum(b%p==0 for b in prefix))
    required = [max(B,after+1)]*slots
    for threshold in sorted({s for s,w in witnesses} | set(tails)):
        w = sorted((w for s,w in witnesses if s>=threshold), reverse=True)
        # Each slot can miss at most k-2 primes. The strongest subset
        # uses all witnesses or just the largest individual count.
        cover = max(w[0], (sum(w)+k-3)//(k-2)) if w else 0
        cover = max(cover,max((need for t,need in tails.items() if t>=threshold), default=0))
        if cover > slots: return Fraction(0)
        for i in range(slots-cover, slots):
            required[i] = max(required[i], threshold)
    # Replace each threshold by the next *possible* smooth denominator,
    # respecting upper bounds on multiplicities of small denominators.
    capcache = {}
    mass = Fraction(1,M0*M0)+sum((Fraction(1,b*b) for b in prefix),Fraction(0))
    previous = 0
    counts = {}
    def next_smooth(lower):
        # At very large thresholds allow every integer (a safe relaxation).
        if lower > 500:
            return lower
        ceiling = 2*lower
        while True:
            ds=divisors_only(primes,ceiling)
            j=bisect_left(ds,lower)
            if j<len(ds):
                return ds[j]
            ceiling *= 2
    for L in required:
        b = next_smooth(max(previous,L))
        while True:
            if b not in capcache:
                c=(slots if q*b*b > 2*rank+2 else
                   multiplicity_cap(q*b*b,b==B,primes,q,local_cap))
                capcache[b]=c-c%2 if odd else c
            if counts.get(b,0)<capcache[b]: break
            b = next_smooth(b+1)
        if not after and not previous and b!=B:
            return Fraction(0)
        counts[b] = counts.get(b,0)+1
        previous = b
        mass += Fraction(1,b*b)
    return mass


def root_possible(primes, heights, q, rank, B, local_cap=True):
    return optimistic_mass(primes,heights,q,rank,B,local_cap) >= q


def three_prime_supports(q, rank, local_cap=True):
    """Enumerate every possible non-BDP three-prime support at exact rank."""
    pp=[p for p in range(2,2*rank+2) if is_prime(p)]
    for P in combinations(pp,3):
        if prod(p for p in P if q%p==0)!=q: continue
        # BDP covers e1<=5, e2<=4, e3<=3 in increasing prime order.
        heights=tuple((e-(q%p==0)+1)//2 for p,e in zip(P,(6,5,4)))
        if any(p>2 and q%p and p-1>rank-1 for p in P): continue
        Bs = [b for b in divisors_only(P,isqrt((rank-1)//q))
              if root_possible(P,heights,q,rank,b,local_cap)]
        if Bs: yield P,Bs


def squarefree(n):
    return all(n%(p*p) for p in range(2,isqrt(n)+1))


def rank_report(rank):
    return {str(q):[[list(P),B] for P,B in roots]
            for q in range(1,rank) if squarefree(q)
            if (roots:=list(three_prime_supports(q,rank)))}


if __name__ == '__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ranks',nargs='+',type=int)
    args=parser.parse_args()
    print(json.dumps({str(r):rank_report(r) for r in args.ranks},indent=2))
