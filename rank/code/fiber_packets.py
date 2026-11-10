"""Packet and rank lower bounds from odd-valuation complementary twist fibers."""
from collections import Counter
from math import prod


def factors(n):
    p=2
    while p*p<=n:
        e=0
        while n%p==0:
            n//=p
            e+=1
        if e:yield p,e
        p=3 if p==2 else p+2
    if n>1:yield n,1


def fiber_size(q, omit=0, dyadic=False):
    """q is the square-free part of the dimension of the (ordinary or reduced) table."""
    odd=[p for p,_ in factors(q) if p%2 and p!=omit]
    if dyadic and q%2==0:return 2*prod((p+1)//2 for p in odd)
    if not odd:return 1
    first,*rest=odd
    return (first+3)//2 * prod((p+1)//2 for p in rest)


def need(q,p,dyadic=False):
    """Required equal-degree active packets at odd support p."""
    return (1 if q%p==0 else 2)*fiber_size(q,omit=p,dyadic=dyadic)


def rank_floor(q,p,odd_dimension=False,dyadic=False):
    """Rank lower bound at an odd support prime p, from all other odd-valued primes."""
    delta=2 if odd_dimension else 1
    phases=(p-1)//2 if q%p==0 else p-1
    return fiber_size(q,omit=p,dyadic=dyadic)*(1+delta*phases)


def inspect(degrees,dyadic=False):
    D=sum(d*d for d in degrees)
    profile=list(factors(D))
    q=prod(p for p,e in profile if e%2)
    counts=Counter(degrees)
    results={}
    for p,e in profile:
        if p==2:continue
        h=(p-1)//2*(2 if D%2 else 1)
        available=sum(m//h for d,m in counts.items() if d>1 and (D//(d*d))%p==0)
        results[str(p)]={'available':available,'required':need(q,p,dyadic),'rank_floor':rank_floor(q,p,D%2,dyadic)}
    return {'dimension':D,'square_class':q,'primes':results}


def obstruction(degrees,dyadic=False):
    result=inspect(degrees,dyadic)
    for p,record in result['primes'].items():
        if record['available']<record['required']:return int(p)
    return None


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path,help='JSON list of types or JSON object with a survivors list')
    ap.add_argument('--output',type=Path,help='write the details for every input type')
    ap.add_argument('--dyadic','--ordinary-dyadic',action='store_true',help='Use the dyadic torus bound at even square-free part.')
    args=ap.parse_args()
    rows=json.loads(args.input.read_text())
    if isinstance(rows,dict):rows=rows['survivors']
    records=[inspect(ds,args.dyadic) for ds in rows]
    totals=Counter(tuple(int(p) for p,r in row['primes'].items() if r['available']<r['required']) for row in records)
    print('Failed-prime sets and counts:',dict(totals))
    if args.output:args.output.write_text(json.dumps(records,indent=2)+'\n')
