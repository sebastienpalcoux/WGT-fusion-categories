"""Count forms and obstruction witnesses on ordinary/reduced JSONL degree types."""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
import json
from math import gcd, isqrt
from pathlib import Path

from support_certificates import b, missing_pair


@lru_cache(maxsize=65536)
def factor(n):
    found = []
    p = 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            found.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        found.append((n, 1))
    return tuple(found)


def one_degree(ds, primes):
    for d in sorted(set(ds) - {1}):
        l = max(p**e for p, e in factor(d))
        for r in primes:
            if d % r and r > b(d//l):
                return d, r
    return None


def audit(path, reduced):
    counts = Counter()
    combinations = Counter()
    for row in Path(path).read_text().splitlines():
        ds = tuple(json.loads(row))
        if ds.count(1) != 1 or any(x < 1 for x in ds):
            raise ValueError("A positive table type with a single unit is required.")
        delta = sum(x*x for x in ds)
        if any(delta % (x*x) for x in ds):
            raise ValueError("Squared divisibility fails on an input type.")
        primes = tuple(p for p, _ in factor(delta))
        a = one_degree(ds, primes)
        c = missing_pair(ds, delta, primes, reduced)
        counts['rows'] += 1
        counts['crossing_4_17'] += bool(a)
        counts['pair_16_1'] += bool(c)
        counts['either'] += bool(a or c)
        combinations[(bool(a), bool(c))] += 1
    return dict(counts), {f"{a},{c}": n for (a, c), n in sorted(combinations.items())}


if __name__ == '__main__':
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('input', nargs='+')
    ap.add_argument('--reduced', action='store_true',
                    help='Input rows give degrees in the fermionic quotient.')
    args = ap.parse_args()
    reports = {}
    for path in args.input:
        counts, both = audit(path, args.reduced)
        reports[Path(path).name] = dict(counts=counts, overlap=both)
    print(json.dumps(reports, indent=2, sort_keys=True))
