"""Audit the 66 supplied residuals. Pass the original rank28_completed.zip."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
from zipfile import ZipFile
import json

from generator_activity import mass_obstruction,saturated_twelfth_floor


def check(path):
    with ZipFile(path) as z:
        raw = z.read('family_result.json')
        report = json.loads(raw)
        candidates = {tuple(json.loads(l)) for l in z.read('candidates.jsonl').splitlines()}
    survivors = list(map(tuple,report['survivors']))
    assert len(candidates) == report['profiles'] == 9458
    assert len(set(survivors)) == len(survivors) == 66 and set(survivors)<=candidates
    cases = []
    for ds in survivors:
        D = sum(x*x for x in ds)
        M = isqrt(D//3)
        g = min(x for x in ds if x>1)
        assert len(ds)==28 and D==3*M*M and ds.count(M//2)==10 and ds.count(M//4)>=2
        # Independent valuation-only computation: 5 divides the codegree
        # exactly when the degree omits at least one copy of 5 in M.
        top = 1
        while M%(5*top)==0:
            top *= 5
        allowed = [x for x in ds if x>1 and x%top]
        capacity = sum(x*x for x in allowed)
        counts = Counter(allowed)
        paired = sum((m//2)*2*d*d for d,m in counts.items())
        assert mass_obstruction(ds,False) == (5,capacity)
        assert 10*g*paired <= 10*g*capacity < 4*D
        B = M//max(allowed)
        rank_floor = 13+2*((3*B*B+5*g-1)//(5*g))
        assert rank_floor==saturated_twelfth_floor(ds,5) and rank_floor>28
        cases.append((g,B,rank_floor,Fraction(10*g*capacity,4*D),Fraction(10*g*paired,4*D)))
    print('family_result.json SHA256',sha256(raw).hexdigest())
    print('verified supplied residuals',len(cases))
    print('smallest certified rank floor',min(x[2] for x in cases))
    print('maximum raw capacity / required mass',max(x[3] for x in cases))
    print('maximum packet capacity / required mass',max(x[4] for x in cases))
    print('minimum eligible denominator histogram',dict(sorted(Counter(x[1] for x in cases).items())))
    print('smallest-degree histogram',dict(sorted(Counter(x[0] for x in cases).items())))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('zip',type=Path)
    check(parser.parse_args().zip)
