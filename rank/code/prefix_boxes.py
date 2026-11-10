"""Rank-specific valuation boxes with a bounded, exact small-denominator prefix."""
from collections import Counter
from fractions import Fraction
from rank_family import divisors_only, multiplicity_cap, optimistic_mass, three_prime_supports


def tilable_prefix(prefix, last, q):
    """Necessary chain recurrence on confirmed denominator counts at p=3,5."""
    count=Counter(prefix)
    for p in (3,5):
        if q%p: continue
        h=(p-1)//2
        for c in range(1,last+1,2):
            if c%p==0: continue
            b=c; incoming=0
            while b<=last:
                remaining=count[b]-2*incoming
                if remaining<0: return False
                incoming=(pow(h,-1,p)*remaining)%p
                if h*incoming>remaining: return False
                b*=2
    return True


def feasible(P, q, R, B, cut):
    heights=tuple((e-(q%p==0)+1)//2 for p,e in zip(P,(6,5,4)))
    choices=[b for b in divisors_only(P,cut) if b>=B]
    def dfs(i,prefix,mass,last):
        if mass>=q or not tilable_prefix(prefix,last,q): return None
        if optimistic_mass(P,heights,q,R,B,prefix=prefix,after=last)<q: return None
        if i==len(choices): return prefix
        b=choices[i]
        cap=min(R-1-len(prefix),multiplicity_cap(q*b*b,b==B,P,q,True))
        step=1 if 2 in P else 2
        for m in range(cap-cap%step,(step-1) if b==B else -1,-step):
            found=dfs(i+1,prefix+(b,)*m,mass+Fraction(m,b*b),b)
            if found is not None: return found
        return None
    return dfs(0,(),Fraction(0),0)


if __name__=='__main__':
    import argparse,json
    from rank_family import squarefree
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('ranks',nargs='+',type=int); ap.add_argument('--q',type=int); ap.add_argument('--cut',type=int,default=8)
    args=ap.parse_args(); out={}
    for r in args.ranks:
        roots=[]
        for q in ([args.q] if args.q else [x for x in range(1,r) if squarefree(x)]):
            for P,Bs in three_prime_supports(q,r):
                for b in Bs:
                    v=feasible(P,q,r,b,args.cut)
                    if v is not None: roots.append({'q':q,'P':P,'B':b,'prefix':v})
        out[r]=roots
    print(json.dumps(out,indent=2))
