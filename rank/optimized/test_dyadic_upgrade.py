#!/usr/bin/env python3
"""Cross-check the two stronger kernels against fixed first-generation rank-31 receipts."""
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path
import json
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import frontier31 as f
import next_rank as nr


def exclusion(row):
    root = isqrt(sum(d*d for d in row))
    denominators = [root//d for d in row]
    n = denominators.count(4)
    if n >= 9 and root % 8:
        return 'height'
    if n == 14 and sum((Fraction(1,b*b) for b in denominators if b%4),Fraction()) >= Fraction(1,16):
        return 'fourteen mass'
    if n == 14:
        cofactor = root
        for p in (2,3,5):
            while cofactor%p == 0:
                cofactor //= p
        if cofactor != 1 or denominators.count(8)<5 or denominators.count(6)>1:
            return 'exceptional row'
    return None


def main():
    chosen = {5,221,327}
    with TemporaryDirectory() as tmp, ZipFile(ROOT/'data/rank31/rank31_seed.zip') as z:
        second, _ = f.build_second(Path(tmp))
        found = set()
        for line in z.open('receipts.jsonl'):
            i, old = json.loads(line)
            if i not in chosen:
                continue
            found.add(i)
            key = old['q'], tuple(old['prefix'])
            before = set(nr.check_receipt(old,31,key))
            new = nr.compute(second,31,key,use_theory=False)
            native = set(nr.check_receipt(new,31,key))
            _, raw = nr.exact_run(second['exact'],'shard',31,[str(key[0]),nr.prefix_text(key[1])])
            exact = set(nr.emit_rows(raw,31,key))
            nr.require(native==exact and native<=before,'kernel disagreement')
            excluded = Counter(exclusion(row) for row in before-native)
            nr.require(None not in excluded,'pruned row not justified by either theorem')
            nr.require(len(native)=={5:0,221:0,327:4}[i], 'unexpected replay count')
            print(f'shard {i}: {len(before)} old, {len(native)} retained, '+str(dict(sorted(excluded.items()))))
        nr.require(found == chosen,'missing test witness')


if __name__ == '__main__':
    main()
