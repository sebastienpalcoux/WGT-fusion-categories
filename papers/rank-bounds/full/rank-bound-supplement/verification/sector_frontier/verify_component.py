"""Check component reservations against synthetic exact phase populations."""
from collections import Counter
from itertools import combinations
from math import prod, isqrt
from random import Random
from component_envelope import envelope, squares
from sector_prefix import next_multiple, smoothseq
from sector_frontier import costs

def verify():
    square_checks=0
    for n in range(2,47):
        values={0}
        for h in range(1,5):
            values={x+a*a for x in values for a in range(1,isqrt(n)+1) if x+a*a<=n}
            if n in values:break
        assert squares(n)==h
        square_checks+=1
    # Physical label multiplicities in exact phase sectors.  For each case
    # validate the hypotheses, then check every closed-denominator cutoff
    # after shifting full Galois bundles to larger eligible denominators.
    scenarios=[
        (1,(2,3,5),((2,3),(5,)),{(2,3):{():1,(2,):1,(3,):4,(2,3):6},
                                      (5,):{():1,(5,):8}}),
        (2,(2,3,5,7),((3,5),(7,)),{(3,5):{():1,(3,):3,(5,):8,(3,5):12},
                                        (7,):{():2,(7,):12}}),
        (3,(2,3,5),((2,5),),{(2,5):{():2,(2,):2,(5,):8,(2,5):8}}),
        (1,(2,3,5,7),((2,3),(3,5),(7,)),{(2,3):{():1,(2,):1,(3,):4,(2,3):6},
                                          (3,5):{():1,(3,):4,(5,):8,(3,5):12},
                                          (7,):{():1,(7,):6}}),
    ]
    rng=Random(45213); checks=0
    for q,S,facets,counts in scenarios:
        C,Bq=costs(q,S)
        F=tuple((frozenset(face),prod(face)) for face in facets)
        for face in facets:
            for k in range(len(face)+1):
                for sector in combinations(face,k):
                    assert counts[face][sector]>=C*prod(p-1 for p in sector)
            for p in face:
                if p==2:continue
                assert sum(m for s,m in counts[face].items() if p in s)>=Bq*(p-1)
                for t in face:
                    if t==p:continue
                    h=4 if t==2 else 2 if t%4==1 else 3
                    assert sum(m for s,m in counts[face].items() if p in s)>=h*(p-1)
                    if p>=5 and t>2:
                        assert sum(m for s,m in counts[face].items() if p in s or t in s)>=(p-1+squares(p+1))*(t-1)+p-1
        for trial in range(14):
            labels=[]
            for face,groups in counts.items():
                for sector,m in groups.items():
                    H=prod((p-1)//2 for p in sector if p!=2)
                    assert m%H==0
                    for _ in range(m//H):
                        factor=face[0] if q==1 and len(face)==1 else rng.choice(S)
                        b=prod(face)*factor**rng.randrange(3)
                        labels.extend([b]*H)
            M=prod(S)*max(labels)*max(labels)
            labels.sort()
            for low in [1]+sorted({b+1 for b in labels[:-1]}):
                hist=Counter(b for b in labels if b<low)
                rest=[b for b in labels if b>=low]+[M]
                ell=smoothseq(S,low)[0]; baseline=[ell]*(len(rest)-1)+[M]
                calculated=envelope(S,F,q,hist,ell,C,Bq,baseline,M,next_multiple)
                assert calculated is not None and all(x<=y for x,y in zip(calculated,rest))
                checks+=1
    return {'square_counts':square_checks,'complete_artificial_prefixes':checks}
if __name__=='__main__':print(verify())
