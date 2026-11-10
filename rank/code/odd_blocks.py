"""Odd-prime block and squared-mass feasibility for ordinary integral modular types."""
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt


def odd_valuation_primes(D):
    # A supplied type has d_max^2|D; factor only its much smaller codegree.
    p = 3
    while D % 2 == 0:
        D //= 2
    while p*p <= D:
        v = 0
        while D % p == 0:
            D //= p
            v += 1
        if v % 2:
            yield p
        p += 2
    if D > 1:
        yield D


def defect_one_tiles(degrees, p):
    """Whether degrees of p-defect one admit the two-orbit block partition."""
    if not degrees:
        return True
    d = sorted(degrees)
    counts = tuple(degrees[y] for y in d)
    index = {y:i for i,y in enumerate(d)}
    tiles = set()
    for i,y in enumerate(d):
        for a in range(1,p):
            if (p-1)%a: continue
            for b in {(p-1)//2,p-1}:
                if p==3 and a+b!=3: continue
                x2, rem = divmod((p-1)*a*y*y,b)
                x = isqrt(x2)
                if rem or x*x!=x2 or x not in index: continue
                tile=[0]*len(d)
                tile[i]+=a
                tile[index[x]]+=b
                if all(t<=m for t,m in zip(tile,counts)):
                    tiles.add(tuple(tile))
    at=[[tile for tile in tiles if tile[i]] for i in range(len(d))]
    @lru_cache(None)
    def search(state):
        if not any(state): return True
        i=min((i for i,m in enumerate(state) if m),key=lambda j:len(at[j]))
        for tile in at[i]:
            rest=tuple(m-t for m,t in zip(state,tile))
            if min(rest)>=0 and search(rest): return True
        return False
    return search(counts)


def signed_sum(groups, target):
    """Solve exact bounded signed-sum problem, using shifted suffix gcd and interval bounds."""
    groups=sorted((sorted(set(g)) for g in groups),key=lambda g:g[-1]-g[0],reverse=True)
    n=len(groups)
    low=[0]*(n+1); high=[0]*(n+1); mod=[0]*(n+1)
    for i in reversed(range(n)):
        g=groups[i]; low[i]=low[i+1]+g[0]; high[i]=high[i+1]+g[-1]; m=mod[i+1]
        for c in g: m=gcd(m,c-g[0])
        mod[i]=m
    @lru_cache(None)
    def search(i,rem):
        if rem<low[i] or rem>high[i]: return False
        if mod[i] and (rem-low[i])%mod[i]: return False
        if i==n: return rem==0
        return any(search(i+1,rem-c) for c in groups[i])
    return search(0,target)


def mass_feasible(high_defect, p, odd, principal):
    # Defect-one blocks may be removed entirely: each satisfies the slice separately.
    step=2 if odd else 1
    h=(p-1)//2
    opts=[]
    for y,m in high_defect.items():
        coefficients=set()
        for t in range(0,m+1,step*h):
            unassigned=m-t
            if principal:
                # Unselected labels and labels of higher twist contribute zero.
                for inactive in range(0,unassigned+1,step):
                    coefficients.add((t-(p-1)*inactive)*y*y)
            else:
                for u in range(0,unassigned+1,step*p*h):
                    coefficients.add((t-(p-1)*(unassigned-u))*y*y)
        opts.append(coefficients)
    return signed_sum(opts,p-1)


def odd_block_obstruction(ds):
    D=sum(y*y for y in ds)
    degrees=Counter(ds)
    del degrees[1]
    for p in odd_valuation_primes(D//max(ds)**2):
        if D%(p*p)==0:
            high={y:m for y,m in degrees.items() if (D//(y*y))%(p*p)==0}
        else:
            # Without squarefullness the vacuum itself can belong to a defect-one block.
            continue
        low={y:m for y,m in degrees.items() if y not in high}
        if not defect_one_tiles(low,p): return 'defect-one tiles at '+str(p)
        if not mass_feasible(high,p,D%2,False): return 'odd-slice mass at '+str(p)
        if not mass_feasible(high,p,D%2,True): return 'principal slice at '+str(p)
    return None


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path,nargs='+',help='JSONL types, or a JSON object with a survivors array')
    args=ap.parse_args()
    for path in args.input:
        s=path.read_text(); rows=json.loads(s)['survivors'] if s.lstrip().startswith('{') else [json.loads(t) for t in s.splitlines() if t]
        cnt=Counter(odd_block_obstruction(t) or 'passes' for t in rows)
        print(path.name,json.dumps(dict(cnt),sort_keys=True))
