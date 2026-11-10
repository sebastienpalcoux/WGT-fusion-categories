"""List square-dimension jobs that survive the rank-envelope root test."""
import argparse
import json
import sys
from itertools import combinations
from math import isqrt
from pathlib import Path

from odd_frontier import odd_smooth_between, primes_through
for base in (Path(__file__).resolve().parent, Path(__file__).resolve().parents[1]/'code',
             Path(__file__).resolve().parents[1]/'package'/'code'):
    if (base/'rank_family.py').exists():
        sys.path.insert(0,str(base))
        break
from rank_family import root_possible


def jobs(cutoff):
    cap = cutoff - (cutoff-1)%8
    primes = [p for p in primes_through(2*cap+1) if 2*(p-1)<=cap-1]
    rows = []
    for k in range(3,len(primes)+1):
        for support in combinations(primes,k):
            heights=(3,3,2) if k==3 else (1,)*k
            heads=[b for b in odd_smooth_between(5,isqrt(cap-1),support)
                   if root_possible(support,heights,1,cap,b)]
            if heads:
                rows.append({'primes':support, 'first':sorted(heads)})
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cutoff',type=int)
    args=parser.parse_args()
    print(json.dumps(jobs(args.cutoff),indent=2))
