#!/usr/bin/env python3
"""Regenerate and screen every necessary degree list through rank 30."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import json
from math import gcd, isqrt, lcm
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

from screen import exclusion_reason

HERE = Path(__file__).resolve().parent


def check_six():
    found = []
    for a, b, c in combinations(range(2, 36), 3):
        if a > 5 or b > 10 or sum(k % 3 == 0 for k in (a, b, c)) < 2:
            continue
        n = lcm(a, b, c)
        if n % 6 or n // c < 6 or abs(n//a - n//b - n//c) != 1:
            continue
        if (2*n + 1 + n//a + n//b + n//c)**2 >= 6*n*n:
            found.append((a, b, c))
    if found != [(3, 6, 7)]:
        raise ValueError('codegree-six subcase did not reproduce')


def screen_job(task):
    binary, r, shard = task
    command = [binary, str(r)]
    if shard is not None:
        command.append(str(shard))
    reasons = Counter()
    with tempfile.TemporaryFile(mode='w+') as diagnostics:
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=diagnostics,
                              text=True, encoding='ascii') as proc:
            try:
                for line in proc.stdout:
                    d = json.loads(line)
                    if len(d) != r or d[0] != 1 or d[1] < 6 or d != sorted(d):
                        raise ValueError(f'malformed degrees {d}')
                    D = sum(x*x for x in d)
                    square = isqrt(D)**2 == D
                    if (square != (binary.name == 'square') or
                            any(D % (x*x) for x in set(d)) or gcd(*d[1:]) != 1):
                        raise ValueError(f'inconsistent dimension list {d}')
                    answer = exclusion_reason(d)
                    if answer is None:
                        raise ValueError(f'surviving rank-{r} profile {d}')
                    reasons[answer] += 1
            except BaseException:
                proc.kill()
                raise
            status = proc.wait()
        diagnostics.seek(0)
        log = diagnostics.read().strip()
    receipt = re.fullmatch(r'rank (\d+), (?:fours (-?\d+), )?nodes (\d+), '
                           r'(?:profiles|candidates) (\d+)', log)
    if (status or not receipt or int(receipt[1]) != r or
            int(receipt[4]) != sum(reasons.values()) or
            (binary.name == 'square' and int(receipt[2]) !=
             (-1 if shard is None else shard))):
        raise RuntimeError(f'incomplete generation for {command}: {log}')
    return r, binary.name, reasons, int(receipt[3])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument('--rank', type=int, choices=range(9, 31),
                        help='check one rank instead of all ranks 9 through 30')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    start = time.monotonic()
    check_six()
    ranks = [args.rank] if args.rank else range(9, 31)
    totals, reasons, nodes = Counter(), Counter(), Counter()
    with tempfile.TemporaryDirectory() as build:
        jobs = []
        for case in ('square', 'nonsquare'):
            exe = Path(build)/case
            subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O3',
                            str(HERE/(case+'.cpp')), '-o', exe], check=True)
            for r in reversed(ranks):
                shards = reversed(range(14)) if case == 'square' and r == 30 else [None]
                jobs.extend((exe, r, shard) for shard in shards)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            pending = [pool.submit(screen_job, task) for task in jobs]
            for result in as_completed(pending):
                r, case, counts, n = result.result()
                totals[r, case] += sum(counts.values())
                reasons.update(counts)
                nodes[case] += n
    # Read the independent reference counts only after the search completes.
    expected = {}
    for line in (HERE/'counts.tsv').read_text().splitlines()[1:]:
        rank, sq, ns = map(int, line.split())
        expected[rank] = (sq, ns)
    print('rank     square     nonsquare    survivors')
    for r in ranks:
        sq, ns = totals[r, 'square'], totals[r, 'nonsquare']
        if (sq, ns) != expected[r]:
            raise ValueError(f'rank {r}: unexpected counts {sq}, {ns}')
        print(f'{r:4} {sq:10} {ns:13} {0:12}')
    print('First exclusions:', dict(sorted(reasons.items())))
    print('Nodes:', dict(nodes))
    print(f'All selected ranks verified in {time.monotonic()-start:.1f} seconds.')


if __name__ == '__main__':
    main()
