#!/usr/bin/env python3
"""Independent exact certificate excluding canonical length ten at5400.

Enumerates all multiplicity partitions, all relevant conductor subfields via
quotient permutation discriminants, and applies the forbidden24/200 pattern
and final degree-three matrix skew determinant. Standard library only.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd,prod
from pathlib import Path
import argparse,json
from verify_two_blocks import (profiles,quadratic_discriminants,squarefree_part,
                              check_modular_implication)
HERE=Path(__file__).resolve().parent
N=5400

@lru_cache(None)
def partitions(k):
    # Multiplicity-count vectors, independent of the root's recursive
    # bounded-part partition enumeration.
    answers=[]
    for counts in product(*(range(k//m+1) for m in range(1,k+1))):
        if sum(m*c for m,c in enumerate(counts,1))==k:
            answers.append(tuple(sorted((m for m,c in enumerate(counts,1) for _ in range(c)),reverse=True)))
    return tuple(sorted(answers,reverse=True))

def kronecker(D,a):
    answer=1
    while a%2==0:
        a//=2
        if D%8 in (3,5):answer=-answer
    D%=a
    while D:
        while D%2==0:
            D//=2
            if a%8 in (3,5):answer=-answer
        D,a=a,D
        if D%4==a%4==3:answer=-answer
        D%=a
    assert a==1
    return answer

@lru_cache(None)
def generated_subgroups(c):
    units={a for a in range(c) if gcd(a,c)==1}
    seen={frozenset([1])};todo=list(seen)
    while todo:
        h=todo.pop()
        for a in units-h:
            k=set(h)|{a}
            while True:
                new=k|{x*y%c for x in k for y in k}
                if new==k:break
                k=new
            k=frozenset(k)
            if k not in seen:seen.add(k);todo.append(k)
    return tuple(sorted(tuple(sorted(h)) for h in seen))

def permutation_sign(p):
    seen=set();cycles=0
    for i in range(len(p)):
        if i in seen:continue
        cycles+=1;j=i
        while j not in seen:seen.add(j);j=p[j]
    return -1 if (len(p)-cycles)%2 else 1

@lru_cache(None)
def subfield_orbits(c,k):
    units=tuple(a for a in range(c) if gcd(a,c)==1)
    ds=quadratic_discriminants(c);options=[]
    for h in generated_subgroups(c):
        e=len(units)//len(h)
        if e>k:continue
        cosets=sorted({tuple(sorted(a*b%c for b in h)) for a in units})
        where={C:i for i,C in enumerate(cosets)}
        signs={a:permutation_sign([where[tuple(sorted(a*b%c for b in C))] for C in cosets]) for a in units}
        if all(s==1 for s in signs.values()):D=1
        else:
            matches=[D for D in ds if all(kronecker(D,a)==s for a,s in signs.items())]
            assert len(matches)==1,(c,h,signs,matches)
            D=matches[0]
        options.append((e,squarefree_part(abs(D)),D,h))
    return tuple(options)

@lru_cache(None)
def block_classes(a,k):
    if not k:return (1,)
    if k==1:return (1,)
    f=N//a;c=gcd(2*f,f*f)
    if k>=3:
        assert c in (9,12,16),(a,k,c)
        options={(e,s) for e,s,D,h in subfield_orbits(c,k)}
    else:
        options={(1,1)}|{(2,squarefree_part(abs(D))) for D in quadratic_discriminants(c)}
    states=[set() for _ in range(k+1)];states[0]={1}
    for used in range(k):
        for old in states[used]:
            for degree,disc in options:
                if used+degree<=k:states[used+degree].add(squarefree_part(old*disc))
    return tuple(sorted(states[k]))

def records():
    out=[]
    for index,p in enumerate(profiles()):
        counts=Counter(p);aa=sorted(counts,reverse=True);base=squarefree_part(prod(p))
        for parts in product(*(partitions(counts[a]) for a in aa)):
            for part in parts:
                cc=Counter(part)
                assert all(c==1 for m,c in cc.items() if m>1 and m%2)
            rank=1+sum(m*m for part in parts for m in part)
            comm=all(m==1 for part in parts for m in part)
            classes={1}
            for a,part in zip(aa,parts):
                classes={squarefree_part(old*new) for old in classes for new in block_classes(a,part.count(1))}
            det=base in classes
            # Original lower bounds suffice because ranks atL10 areeven.
            low=rank<(11 if comm else 13)
            pattern=counts[24]==counts[200]==1 and all(a in (24,200) or a%225==0 for a in aa)
            remaining=det and not low and not pattern
            out.append({'index':index,'profile':p,'multiplicities':list(zip(aa,parts)),
                        'rank':rank,'commutative':comm,'baseline':base,
                        'linear_field_squareclasses':sorted(classes),'determinant_pass':det,
                        'excluded_by_inherited_rank':low,'forbidden_24_200':pattern,
                        'remaining_before_skew':remaining})
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'length10_verification.json')
    args=parser.parse_args();rr=records()
    assert len(rr)==208
    after_det=[r for r in rr if r['determinant_pass']]
    after_rank=[r for r in after_det if not r['excluded_by_inherited_rank']]
    after_pattern=[r for r in after_rank if not r['forbidden_24_200']]
    assert (len(after_det),len(after_rank),len(after_pattern))==(37,21,1)
    final=after_pattern[0]
    assert final['index']==10 and final['rank']==16
    assert final['multiplicities']==[(900,(3,1,1)),(540,(1,)),(200,(1,)),(135,(1,)),(24,(1,))]
    # The final pair oflinear f6characters must have discriminantclass3.
    final_fields=[D for D in quadratic_discriminants(12) if squarefree_part(abs(D))==3]
    assert final_fields==[-3,12]
    branches=[]
    for D in final_fields:
        q=int(D<0);b=q+3
        linear_skew=Fraction(2*abs(D if D%4==1 else D//4),6) if D<0 else Fraction(1)
        matrix_skew=Fraction(1,3)**3
        total=squarefree_part((linear_skew*matrix_skew).numerator*(linear_skew*matrix_skew).denominator)
        target=squarefree_part(2**b)
        assert total!=target
        branches.append({'quadratic_discriminant':D,'source_dual_pairs':b,
                         'total_skew_squareclass':total,'required_skew_squareclass':target,
                         'excluded':True})
    check_modular_implication()
    result={'status':'PASS','profiles':26,'multiplicity_assignments':208,
            'after_determinant':37,'after_inherited_rank_bounds':21,
            'after_forbidden_24_200_pattern':1,'after_matrix_skew_determinant':0,
            'final_case':final,'final_skew_branches':branches,
            'field_block_classes':{f'{a}:{k}':block_classes(a,k) for a,k in sorted({(a,part.count(1)) for r in rr for a,part in r['multiplicities'] if part.count(1)>=2})},
            'records':rr,
            'conclusion':'Canonical length10 is impossible under the standing hypotheses; the inherited length bound therefore strengthens to L>=11.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','field_block_classes')},indent=2))
if __name__=='__main__':main()
