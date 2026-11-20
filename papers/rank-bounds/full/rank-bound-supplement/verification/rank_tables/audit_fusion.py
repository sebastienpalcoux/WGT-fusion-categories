"""Independent, small degree-partition and phase-cap audit of rejection witnesses.

Usage: python audit_fusion.py path/to/reference-census.zip. This checks rejections on saved rows; it does
not generate degree lists. No third-party packages or number-field routines.
"""
import argparse
from collections import Counter
from gzip import decompress
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path
from zipfile import ZipFile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'census'))
from fusion_support import bounded_sum, fusion_failure


def require(condition, message):
    if not condition:
        raise ValueError(message)


def independent_witness(row, proof):
    """Rebuild every allowed per-label multiplicity by scanning integers."""
    nonunit = row[1:]
    a, b, vacuum = (proof[k] for k in ('a', 'b', 'vacuum'))
    require(a in nonunit and b in nonunit, 'missing product input')
    require((vacuum == 1 and a == b) or
            (vacuum == 0 and (a != b or nonunit.count(a) >= 2)),
            'unrealized vacuum multiplicity')
    dimension = sum(x*x for x in row)
    prime = next((t for t in range(2, isqrt(dimension)+1) if dimension % t == 0), dimension)
    m = min(nonunit)
    target = a*b-vacuum
    require(target == proof['target'], 'incorrect target')
    supplies = []
    for c in sorted(set(nonunit)):
        if c > target:
            continue
        upper = 0
        for n in range(1, target//c + 1):
            feasible = True
            for x,y,z in ((a,b,c),(a,c,b),(b,c,a)):
                rest = x*y - n*z
                if (rest < (1 if x==y else m) or
                        (gcd(x,y) == 1 and n*z*prime > x*y)):
                    feasible = False
                    break
            if not feasible or (not vacuum and a*b-n*c<m):
                break
            upper = n
        if upper:
            supplies.append([c, upper*nonunit.count(c)])
    require(proof['capacities'] == supplies, 'incorrect coin bounds')
    reachable = {0}
    for c, upper in supplies:
        reachable = {v+k*c for v in reachable for k in range(min(upper, (target-v)//c)+1)}
    require(target not in reachable, 'rejected a feasible degree partition')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    # Compare both outcomes of the fast bounded sum against direct enumeration.
    arithmetic_checks = 0
    for limits in product(range(4), repeat=3):
        supplies = list(zip((3, 5, 7), limits))
        attainable = {3*i+5*j+7*k for i,j,k in product(*(range(n+1) for n in limits))}
        for target in range(51):
            require(bounded_sum(target, supplies) == (target in attainable), 'bounded-sum regression')
            arithmetic_checks += 1
    by_rank = {}
    with ZipFile(args.archive) as z:
        for rank in range(1,45):
            counts = Counter()
            for raw in decompress(z.read(f'ranks/{rank:02}/raw-census.jsonl.gz')).splitlines():
                r, _, _, _, row = json.loads(raw)
                require(r == rank and len(row) == rank, 'rank mismatch')
                proof = fusion_failure(row)
                counts['rows'] += 1
                if proof is not None:
                    independent_witness(row, proof)
                    counts['excluded'] += 1
            if counts:
                by_rank[rank] = dict(counts)
    print(json.dumps({'bounded_sum_tests': arithmetic_checks, 'reference_ranks': by_rank}, indent=2))


if __name__ == '__main__':
    main()
