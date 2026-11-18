#!/usr/bin/env python3
"""Standard-library certificate for the rank14 two-degree-two branch.

Independent enumeration of canonical multiplicities, conductor-only linear
field determinant choices, and the five final forbidden 24/200 patterns.
No optimization or source-row feasibility computation is required.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd,prod,isqrt
from pathlib import Path
import argparse,json
N=5400
HERE=Path(__file__).resolve().parent

def prime_support(n):
    return sum(n%p==0 for p in (2,3,5))

PALETTE=tuple(sorted({2**i*3**j*5**k for i,j,k in product(range(4),range(4),range(3))
                      if 10<=2**i*3**j*5**k<=900 and prime_support(2**i*3**j*5**k)>=2},reverse=True))

def profiles():
    # Count-vector recursion, separately arranged from the attachment's
    # entry-by-entry recursion. All degrees are nonincreasing on output.
    @lru_cache(None)
    def solve(i,slots,total):
        if slots==0:return ((),) if total==0 else ()
        if i==len(PALETTE) or total<10*slots or total>PALETTE[i]*slots:return ()
        a=PALETTE[i];answer=[]
        for count in range(min(slots,total//a)+1):
            rem=total-count*a
            for tail in solve(i+1,slots-count,rem):answer.append((a,)*count+tail)
        return tuple(answer)
    return tuple(sorted(solve(0,9,N-1),reverse=True))

def squarefree_part(n):
    n=abs(n);out=1;p=2
    while p*p<=n:
        parity=0
        while n%p==0:n//=p;parity^=1
        if parity:out*=p
        p+=1
    return out*n

def prime_factors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out

def quadratic_discriminants(c):
    # Every squarefree radicand has prime support in c. Enumerate this
    # support instead of scanning all integers in [-c,c].
    ps=prime_factors(c);ds=[]
    for mask in range(1<<len(ps)):
        rad=prod(p for i,p in enumerate(ps) if (mask>>i)&1)
        for sign in (-1,1):
            d=sign*rad
            if d==1:continue
            D=d if d%4==1 else 4*d
            if c%abs(D)==0:ds.append(D)
    return tuple(sorted(ds))

def phi(n):
    return sum(gcd(x,n)==1 for x in range(1,n+1))

def enumerate_records():
    ps=profiles();assert len(ps)==26
    records=[]
    for index,p in enumerate(ps):
        counts=Counter(p)
        candidates=sorted([a for a,c in counts.items() if c>=2],reverse=True)
        for i,a in enumerate(candidates):
            for b in candidates[i:]:
                doubles=(a,b);nonlinear=Counter(doubles)
                if any(2*j>counts[d] for d,j in nonlinear.items()):continue
                linear={d:c-2*nonlinear[d] for d,c in counts.items() if c>2*nonlinear[d]}
                assert sum(linear.values())==5
                baseline=squarefree_part(prod(d**c for d,c in linear.items()))
                repeated=[(d,c) for d,c in linear.items() if c>1]
                assert len(repeated)<=1
                choices=[1]
                if repeated:
                    d,c=repeated[0];f=N//d
                    assert c in (2,3)
                    if c==3:assert phi(2*f)%3!=0
                    choices += list(quadratic_discriminants(2*f))
                surviving=[D for D in choices if squarefree_part(baseline*abs(D))==1]
                pattern=False
                if surviving:
                    assert linear.get(24)==linear.get(200)==1
                    assert 24 not in nonlinear and 200 not in nonlinear
                    assert all(d in (24,200) or d%225==0 for d in counts)
                    pattern=True
                records.append({'index':index,'profile':p,'nonlinear_degrees':doubles,
                                'linear_counts':sorted(linear.items(),reverse=True),
                                'linear_baseline_squareclass':baseline,
                                'repeated_linear_block':repeated,
                                'conductor_only_discriminant_choices':choices,
                                'determinant_survivors':surviving,
                                'forbidden_24_200_pattern':pattern,
                                'excluded':not surviving or pattern})
    return records

def check_modular_implication():
    # Exact residue proof: the two congruences force 75|(d-h).
    residues=[]
    for d,h,g in product(range(9),repeat=3):
        if (d+24*h+200*g)%9==0 and (d*d+24*h*h+200*g*g)%9==0:
            assert (d-h)%3==0
            residues.append((d,h,g))
    for d,h,g in product(range(25),repeat=3):
        if (d+24*h+200*g)%25==0:assert (d-h)%25==0
    # The residual object has dimension 23, so a source constituent must have
    # 1<=d<=23; even this weak range is incompatible with strict|h|<d.
    assert not [(d,h) for d in range(1,24) for h in range(-d+1,d) if (d-h)%75==0]
    return len(residues)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'two_blocks_verification.json')
    args=parser.parse_args();records=enumerate_records()
    assert len(records)==32
    survivors=[(r['index'],tuple(r['nonlinear_degrees'])) for r in records if r['determinant_survivors']]
    assert survivors==[(12,(900,675)),(20,(900,900)),(20,(900,450)),(25,(900,675)),(25,(675,675))]
    assert all(r['excluded'] for r in records)
    residue_count=check_modular_implication()
    result={'status':'PASS','dimension':N,'complex_algebra':'C^6 + M2(C) + M2(C)',
            'canonical_length':10,'palette':PALETTE,'profiles':26,'placements':32,
            'after_conductor_only_determinant':5,'after_forbidden_pattern':0,
            'mod9_solution_classes_checked':residue_count,
            'survivors_before_pattern':survivors,'records':records,
            'scope':'Excludes this rank14 algebra shape under the Round32 center hypotheses; the other rank14 shape is untouched.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','palette')},indent=2))
if __name__=='__main__':main()
