"""Exact rational row test for a top singleton in a ternary mixed block."""
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt
import json
from pathlib import Path


def test(degrees):
    n = max(degrees)
    D = sum(d*d for d in degrees)
    if n%2 or D != 12*n*n or degrees.count(1)!=1:
        return None
    ternary = Counter(d for d in degrees if D//(d*d)%9 != 0)
    if ternary != {n:10, n//2:2}:
        return None
    groups = [(0, 1, [(0,0,0), (0,3,2*n)])]
    for d,m in sorted(Counter(degrees).items()):
        if d in (1,n//2,n):
            continue
        step = d//gcd(d,n)
        states = {(0,0):0}
        for _ in range(m):
            updated = {}
            for (signed,cost),value in states.items():
                limit = isqrt(11-cost)//step
                for j in range(-limit,limit+1):
                    u = j*step
                    key = (signed+u,cost+u*u)
                    updated[key] = max(updated.get(key,0),value+abs(u))
            states = updated
        choices = [(d*sgn,cost,d*value) for (sgn,cost),value in states.items()]
        groups.append((d,step,choices))
    suffix_gcd = [0]*(len(groups)+1)
    for i in reversed(range(len(groups))):
        suffix_gcd[i] = gcd(suffix_gcd[i+1],groups[i][0]*groups[i][1])
    @lru_cache(None)
    def best(i, target, budget):
        if i==len(groups):
            return 0 if target==0 and budget==0 else -1
        if target%suffix_gcd[i]:
            return -1
        score = -1
        for signed,cost,absolute in groups[i][2]:
            if cost<=budget:
                tail = best(i+1,target-signed,budget-cost)
                if tail>=0:
                    score = max(score,tail+absolute)
        return score
    best_value = best(0,-1,11)
    if best_value<0:
        return 'orthogonality'
    if (1+best_value)**2<D:
        return 'modular inequality'
    return 'open'


if __name__=='__main__':
    import argparse,zipfile
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,help='completed zip, family-result JSON, or types JSONL')
    path=parser.parse_args().input
    if path.suffix=='.zip':
        with zipfile.ZipFile(path) as bundle:
            types=json.loads(bundle.read('family_result.json'))['survivors']
    elif path.suffix=='.json':
        types=json.loads(path.read_text())['survivors']
    else:
        types=[json.loads(line) for line in path.read_text().splitlines() if line]
    print(dict(Counter(test(t) or 'outside hypotheses' for t in types)))
