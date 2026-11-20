"""Complete short enumeration of mixed four-prime omission signatures at44.

This certifies only the mixed signatures, not the all-primary Gauss roots.
The new isolation and sparse-height proofs must be applied first.
"""
import sys,json,hashlib
from pathlib import Path
here=Path(__file__).resolve().parent
from large_four import preliminary,profiles,MIXED

initial=[];kept=[];removed=[]
for m,five,support,modes in preliminary(44):
    if MIXED not in modes:continue
    if next(profiles(44,m,five,support,modes),None) is None:continue
    item={'rows':m,'five':five,'support':support,'modes':modes}
    initial.append(item)
    if 5 in support and modes[support.index(5)]==MIXED:
        assert (m,five,support,modes) in (
            (10,False,(3,5,7),(0,1,2)),(11,False,(3,5,7),(0,1,2)))
        removed.append(item)
    else:
        assert modes[0]==MIXED and MIXED not in modes[1:]
        kept.append(item)
assert (len(initial),len(removed),len(kept))==(41,2,39)
from collections import Counter
out={'rank':44,'scope':'mixed four-prime omission signatures only',
     'initial':len(initial),'excluded_by_two_five_lemma':removed,
     'remaining':len(kept),'counts_by_four_multiplicity':dict(Counter(c['rows'] for c in kept)),
     'signatures':kept}
here.joinpath('high-four44-mixed.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k!='signatures'})
