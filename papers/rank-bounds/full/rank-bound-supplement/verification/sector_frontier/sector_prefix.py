"""Extend each rank-45 reservation by complete classes; stop at fixed depth."""
from collections import Counter,defaultdict
from functools import lru_cache
from fractions import Fraction
from math import gcd,lcm,prod,isqrt
import argparse,json
from neutral_cap import improve
from component_envelope import envelope
from sector_frontier import frontier,costs,_part,_free_min,primes_to

primes=primes_to(45)

@lru_cache(maxsize=None)
def smoothseq(S,low):
    # Only four next eligible values needed for a cut or denominator floor.
    result=[]
    n=low
    while len(result)<4:
        m=n
        for p in S:
            while m%p==0:m//=p
        if m==1:result.append(n)
        n+=1
    return tuple(result)

def next_multiple(S,low,k):
    return k*smoothseq(S,(low+k-1)//k)[0]

@lru_cache(maxsize=None)
def cut_floor(S,mask,low):
    left=[p for i,p in enumerate(S) if mask>>i&1]
    right=[p for p in S if p not in left]
    ks=[p*t for p in left for t in right]
    first=set()
    for k in ks:
        first.update(k*t for t in smoothseq(S,(low+k-1)//k)[:3])
    return tuple(ks),tuple(sorted(first)[:3])

def isolated_primes_with_class(S,F,q,hist):
    return () if q!=1 else tuple(p for a,p in F if len(a)==1 and p%2 and p in hist)


def bound(S, F, q, hist, low, floors, Mfloor):
    """Optimistic upper mass on all remaining positions, including the unit."""
    total=sum(hist.values());h=45-total
    if h<=0:return None
    ell=smoothseq(S,low)[0]
    baseline=[ell]*(h-1)+[max(ell,Mfloor)]
    C,Bq=costs(q,S)
    rows=[]
    def charge(n,x):
        return [x]*max(0,n)
    activity=[]
    for p in S:
        if p==2 or q%p==0:continue
        packet=(p-1)//2
        used=sum((n//packet)*packet for b,n in hist.items() if b%p==0)
        activity+=charge(Bq*(p-1)-used,next_multiple(S,low,p))
    rows.append(activity)
    mandatory=[];pure_entries=[]
    for part,k in F:
        isolated=q==1 and len(part)==1 and k%2
        if isolated:
            used=sum(n for b,n in hist.items() if _part(b,k)==b)
            floor=1
            while floor<low:floor*=k
            group=charge(k-used, floor)
            pure_entries+=group
        else:
            used=sum(n for b,n in hist.items() if b%k==0)
            group=charge(C*k-used,next_multiple(S,low,k))
        mandatory+=group
    rows.append(mandatory)
    rows.append([Mfloor]+pure_entries)
    # At valuation one the proper p-free classes are all neutral, and
    # are disjoint from pure labels at every isolated prime.
    for p in isolated_primes_with_class(S,F,q,hist):
        distinct=sum(1 for c in hist if c!=p and c%p==0)
        n=max(0,3-distinct)
        if not n:continue
        low_k=(low+p-1)//p
        vals=smoothseq(tuple(x for x in S if x!=p),max(2,low_k))[:n]
        witnesses=[p*x for x in vals]
        if p==5:
            # The four proper signs have k congruent to 1 or -1 mod 5.
            # Distinct and total requirements can be imposed separately.
            good=[]
            k=max(2,low_k)
            allowed=tuple(x for x in S if x!=p)
            while len(good)<n:
                k=smoothseq(allowed,k)[0]
                if k%5 in (1,4):good.append(p*k)
                k+=1
            witnesses=good
            missing=max(0,4-sum(v for c,v in hist.items() if c!=p and c%p==0))
            rows.append([Mfloor]+pure_entries+[good[0]]*missing)
        rows.append([Mfloor]+pure_entries+witnesses)
    for mask in range(1,2**len(S)-1,2):
        ks,least=cut_floor(S,mask,low)
        known=sum(any(b%k==0 for k in ks) for b in hist)
        rows.append([Mfloor]+pure_entries+list(least[:max(0,3-known)]))
    charges=[]
    for t,power in floors.items():
        need=max(0,_free_min(t)-sum(n for b,n in hist.items() if b%power==0))
        k=power*max(x for x in S if x!=t)
        need_large=not any(b%k==0 for b in hist)
        group=charge(max(0,need-int(need_large)),next_multiple(S,low,power))
        if need_large:group.append(next_multiple(S,low,k))
        rows.append([Mfloor]+group)
        charges+=group
    rows.append([Mfloor]+sorted(charges,reverse=True)[::len(S)-2])
    for row in rows:
        if len(row)>h:return None
        row=sorted([ell]*(h-len(row))+row)
        baseline=list(map(max,baseline,row))
    baseline=envelope(S,F,q,hist,ell,C,Bq,baseline,Mfloor,next_multiple)
    if baseline is None:return None
    baseline=improve(S,F,q,hist,floors[2],ell,Mfloor,baseline,C,Bq,next_multiple)
    if baseline is None:return None
    return sum((Fraction(n,d*d) for d,n in Counter(baseline).items()),Fraction())

def sphere_possible(S,q,hist,divisor):
    # Simultaneous unknown M/divisor modulo each principal prime, retaining
    # zero unless valuation one is forced by a present prime denominator.
    for p in S:
        if q==1 and p in hist and _part(divisor,p)>p:
            return False
        closed=[(b,n) for b,n in hist.items() if _part(q*b*b,p)==q*b*b]
        if not closed:continue
        residues=range(1,p) if q==1 and p in hist else range(p)
        good=False
        for residue in residues:
            for b,m in closed:
                groups=[m]
                if q==1 and p==b and m==p:
                    # There is one neutral member and two square-class orbits.
                    groups += [1,(p-1)//2]
                for a in groups:
                    cost=0
                    for c,n in hist.items():
                        u=b//gcd(b,c)
                        if u%p and c%p==0:
                            g=(divisor//lcm(b,c))*residue
                            z=(a*g)%p
                            cost+=n*u*u*min(z,p-z)**2
                    if cost>a*(q*b*b-a):break
                else:continue
                break
            else:
                good=True
                break
        if not good:return False
    return True

@lru_cache(maxsize=None)
def _block_shapes(t):
    # Entries (labels at c, denominator multiplier, labels at k*c).
    if t==2:return ((2,1,0),)
    vals=set()
    for a in range(1,t):
        if (t-1)%a:continue
        for b in {(t-1)//2,t-1}:
            numerator=(t-1)*a
            if numerator%b:continue
            k=isqrt(numerator//b)
            if k*k*b!=numerator:continue
            vals.add((a+b,1,0) if k==1 else (b,k,a))
    return tuple(sorted(vals))


def block_prefix(t,hist,low,slots,forbidden_equal=frozenset()):
    labels=tuple(sorted((b,n) for b,n in hist.items() if b%t))
    shapes=_block_shapes(t)
    def visit(i,demand):
        if i==len(labels):
            return all(c>=low for c,n in demand.items() if n) and sum(demand.values())<=slots
        b,m=labels[i]
        needs=demand.copy()
        left=m-needs.pop(b,0)
        if left<0:return False
        def tile(j,remaining):
            if j==len(shapes):
                return not remaining and visit(i+1,needs)
            consume,k,other=shapes[j]
            if k==1 and b in forbidden_equal:
                return tile(j+1,remaining)
            for copies in range(remaining//consume+1):
                old=needs.get(k*b,0)
                if copies and other:needs[k*b]=old+copies*other
                if tile(j+1,remaining-copies*consume):return True
                if copies and other:
                    if old:needs[k*b]=old
                    else:del needs[k*b]
            return False
        return tile(0,left)
    return visit(0,{})


def closure(S,F,q,hist,low,floors):
    isolated={w for part,w in F if q==1 and len(part)==1 and w%2}
    for t in S:
        if q%t:continue
        # At q=t a full equal t-block at an odd prime denominator p
        # supplies an effective row of codegree p^2. Its degree is p-free.
        forbidden=frozenset(p for p in hist if p in S and p!=2 and p!=t
                            and floors[p]>p) if q==t else frozenset()
        if not block_prefix(t,hist,low,44-sum(hist.values()),forbidden):
            return False
    for b,m in hist.items():
        f=q*b*b
        for p in S:
            if _part(f,p)==f and (m>f-isqrt(p)-1 or m==f-p and f>=p*p):
                return False
        if q!=1 or b not in S:continue
        p=b
        if p<5 or p not in isolated or floors[p]>p:return False
        if m%p:return False
        j=m//p
        if j>=p-1 or (p-1)%j or p==5 and j!=1 or p==7 and j>2:
            return False
        proper=[(c,n) for c,n in hist.items() if c!=p and c%p==0]
        cost=0
        for c,n in proper:
            a=pow(c//p,-1,p)
            lifts=(a,a-p)
            if n==1:lifts=tuple(v for v in lifts if v%j==0)
            if not lifts:return False
            cost+=n*min(v*v for v in lifts)
        if cost+max(0,3-len(proper))>j*(p-j):return False
        if p==7 and j==2:
            # The three proper degrees have multiplicity 2,1,1, with
            # (absolute) Weyl coefficients 1,2,2 respectively.
            if len(proper)>3 or sum(n==1 for c,n in proper)>2:return False
            if sum(n==2 for c,n in proper)>1:return False
            if any(n not in (1,2) or c//p%7 not in
                   ((1,6) if n==2 else (3,4)) for c,n in proper):return False
    for p in isolated:
        if p>=low or p in hist:continue
        extra=sum(k for a,k in F if k!=p)
        if 45-1-extra<p*(p-1):
            if any(c%p==0 and c%(p*p) for c in hist):return False
            floors[p]=max(floors[p],p**3)
    return True


def extend(plan,depth,counts,exception=None):
    q,B,m,S,F=plan
    # With four dimension primes in nonsquare dimension, at most one
    # unramified prime can have v_p(M)=1. The designated exception is
    # permitted to have larger height. It remains fixed along this prefix.
    initial={p:max(_part(B,p),p**(1+(len(S)==3 or
                   (exception is not None and q%p and p!=exception)))) for p in S}
    def visit(hist,nsteps,prefix_mass):
        counts['prefix_tests']+=1
        low=max(hist)+1
        floors={p:max(initial[p],*(_part(b,p) for b in hist)) for p in S}
        if not closure(S,F,q,hist,low,floors):return False
        divisor=prod(floors.values())
        if not sphere_possible(S,q,hist,divisor):return False
        # At least two support primes in every nonunit degree.
        Mfloor=max(divisor,prod(S[:2])*max(hist))
        upper=bound(S,F,q,hist,low,floors,Mfloor)
        if upper is None or prefix_mass+upper<q:return False
        if nsteps==0:return True
        N=45-sum(hist.values())-1
        if N<=0:return False  # This routine stops well before termination.
        residual=Fraction(q)-prefix_mass
        if residual<=0:return False
        # M >= 6*b for the next denominator b. All N nonunits are >=b.
        mx=isqrt(((36*N+1)*residual.denominator)//(36*residual.numerator))
        b=smoothseq(S,low)[0]
        while b<=mx:
            for n in range(1,N+1):
                mass=prefix_mass+Fraction(n,b*b)
                if mass>=q:break
                if n==N:
                    # A last closed class leaves only the unit.
                    rem=Fraction(q)-mass
                    if rem.numerator==1 and isqrt(rem.denominator)**2==rem.denominator:
                        return True
                    continue
                if visit({**hist,b:n},nsteps-1,mass):return True
            b=smoothseq(S,b+1)[0]
        return False
    return visit({B:m},depth,Fraction(m,B*B))

def run(depth):
    cover=frontier(include_plans=True);by=defaultdict(Counter); counts=Counter();left=defaultdict(int)
    for q,b,m,S,F in cover:
        exceptions=[p for p in S if q%p] if q!=1 and len(S)==4 else [None]
        if any(extend((q,b,m,S,F),depth,counts,p) for p in exceptions):
            left[q,b,m,S]+=1
    for (q,b,m,S),n in left.items():
        by[f'{q},{b}']['plans']+=n
        by[f'{q},{b}']['starts']+=1
    return dict(depth=depth,initial_plans=len(cover),prefix_tests=counts['prefix_tests'],
                starts=len(left),plans=sum(left.values()),by_q_b=dict(by),
                surviving_starts=[[q,b,m,list(S),n] for (q,b,m,S),n in left.items()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('depth',type=int,choices=(0,1,2))
    p.add_argument('output')
    args=p.parse_args();result=run(args.depth)
    with open(args.output,'w') as f:json.dump(result,f,indent=2);f.write('\n')
    print({k:v for k,v in result.items() if k!='surviving_starts'})
