"""Even-valuation primary square classes, including higher-level activity."""
from collections import Counter
from functools import lru_cache
from math import gcd


def max_mass(D, mult, p):
    """Upper bound on the active squared mass in a p-inertia stable block union/fiber."""
    step=((p-1)//2)*(2 if D%2 else 1)
    choices=[]
    for d,m in mult.items():
        if d==1 or (D//(d*d))%p or m<step:continue
        n=m//step
        high=(D//(d*d))%(p*p)==0
        # Signed imbalance k and largest active count compatible with k.
        options={}
        for k in range(-n,n+1):
            left=n-abs(k)
            # Idle slots can pair off, or form a group of p for high twist.
            missing=left%2 and (not high or left<p)
            options[k*d*d]=step*(n-int(missing))*d*d
        choices.append(options)
    def span(G):return max(G)-min(G)
    choices.sort(key=span,reverse=True)
    n=len(choices)
    upper=[0]*(n+1);div=[0]*(n+1)
    for i in reversed(range(n)):
        upper[i]=upper[i+1]+max(choices[i]); div[i]=gcd(div[i+1],*choices[i])
    @lru_cache(None)
    def visit(i,target):
        if abs(target)>upper[i] or (div[i] and target%div[i]):return -1
        if i==n:return 0 if target==0 else -1
        best=-1
        for k,active in choices[i].items():
            z=visit(i+1,target-k)
            if z>=0:best=max(best,z+active)
        return best
    return visit(0,0)


def fiber_mass(D, mult, p, odd_primes, state_limit=100000):
    """Conditional prime-fiber bound, or None if the general fallback is needed."""
    from itertools import product
    from math import prod
    step=((p-1)//2)*(2 if D%2 else 1)
    groups=[(d,m//step) for d,m in mult.items() if d!=1 and m>=step and (D//(d*d))%p==0]
    if any(n>=p and (D//(d*d))%(p*p)==0 for d,n in groups):return None
    if prod(n+1 for d,n in groups)>state_limit:return None
    ds=[d for d,n in groups]; caps=[n for d,n in groups]
    allowed={}
    for packet in product(*(range(n+1) for n in caps)):
        if not any(packet):continue
        w=sum(a*d*d for a,d in zip(packet,ds))
        g=gcd(*(d for a,d in zip(packet,ds) if a))
        if any((w//(g*g))%t for t in odd_primes):continue
        allowed.setdefault(w,[]).append(packet)
    for w in sorted(allowed,reverse=True):
        for i,left in enumerate(allowed[w]):
            for right in allowed[w][i:]:
                if all(a+b<=n for a,b,n in zip(left,right,caps)):return 2*step*w
    return 0


def verdict(ds,p,projections=None,odd_primes=()):
    D=sum(d*d for d in ds); counts=Counter(ds); x=min(d for d in counts if d>1)
    act=max_mass(D,counts,p)
    if act==0:return 'empty'
    if odd_primes:
        crossed=fiber_mass(D,counts,p,odd_primes)
        if crossed==0:return 'empty crossed fiber'
        if crossed is not None and 2*p*x*crossed<(p-1)*D:return 'crossed mass'
    if 2*p*x*act<(p-1)*D:return 'mass'
    if projections is not None:
        for pairs in projections:
            local=Counter(dict(pairs));local[1]=1
            split=max_mass(D,local,p)+max_mass(D,counts-local,p)
            if 2*p*x*split>=(p-1)*D:break
        else:return 'projected mass'
    return 'passes'


def prefix_balance(q, p, closed, remaining, next_min, odd=False):
    """Even-valuation necessary test on sealed runs b:count for a known p|M, p∤q.

    next_min strictly exceeds all closed denominators. One might have zero
    remaining slots (then next_min is just any larger integer).
    """
    from math import lcm
    from families import factors
    step=(p-1)//2*(2 if odd else 1)
    if q%p==0 or p==2 or remaining<0: raise ValueError('expected odd p of even valuation and nonnegative remaining count')
    if remaining>=2*step:return True
    denominators=[]
    for b,m in closed.items():
        if b%p:continue
        n=m//step
        if n>=2:return True
        if n:denominators.append(b)
    if not denominators:return False
    L=lcm(*denominators)
    sums={0}
    for b in denominators:
        w=(L//b)**2
        shifted={s+w for s in sums}
        if not sums.isdisjoint(shifted):return True
        sums.update(shifted)
    if remaining<step:return False
    ds=[1]
    for ell,v in factors(L):
        ds=[a*ell**k for a in ds for k in range(v+1)]
    for b in ds:
        if b<next_min or b%p:continue
        w=(L//b)**2
        if any(s+w in sums for s in sums):return True
    return False


if __name__ == '__main__':
    import argparse, json
    from pathlib import Path
    from families import factors
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('types',type=Path,help='A JSON list of types, a survivors object, or JSONL.')
    parser.add_argument('--project',action='store_true',help='Also test full principal block projections.')
    parser.add_argument('--fibers',action='store_true',help='Apply the other odd-valuation primes on each phase fiber.')
    args=parser.parse_args()
    raw=args.types.read_text()
    try:
        rows=json.loads(raw)
        if isinstance(rows,dict):rows=rows['survivors']
        elif rows and isinstance(rows[0],int):rows=[rows]
    except json.JSONDecodeError:
        rows=[json.loads(s) for s in raw.splitlines() if s.strip()]
    result=Counter()
    for ds in rows:
        D=sum(d*d for d in ds)
        outcome=[]
        for p,e in factors(D):
            if p==2 or e%2:continue
            blocks=None
            if args.project:
                from principal import projections
                blocks=projections(ds,p,e)
            odd_primes=tuple(t for t,a in factors(D) if t%2 and a%2) if args.fibers else ()
            outcome.append((p,verdict(ds,p,blocks,odd_primes)))
        result[tuple(outcome)]+=1
    for k,n in sorted(result.items()):print(n,k)
