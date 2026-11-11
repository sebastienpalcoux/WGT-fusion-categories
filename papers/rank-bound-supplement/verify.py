#!/usr/bin/env python3
"""Regenerate and screen degree profiles in every rank through 31.

Ranks 9–30 use their small-rank generators. Rank 31 uses a prefix-partitioned
generator with exact-integer continuations and resumable receipts. Both routes
apply the same terminal necessary conditions in screen.py.
"""
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
import sys
import tempfile

from screen import exclusion_reason

HERE = Path(__file__).resolve().parent


def check_six():
    """Reproduce the small codegree-six calculation used before enumeration."""
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


def screen_small_job(task):
    binary, rank, shard = task
    command = [str(binary), str(rank)]
    if shard is not None:
        command.append(str(shard))
    reasons = Counter()
    with tempfile.TemporaryFile(mode='w+') as diagnostics:
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=diagnostics,
                              text=True, encoding='ascii') as proc:
            try:
                for line in proc.stdout:
                    degrees = json.loads(line)
                    if (len(degrees) != rank or degrees[0] != 1 or
                            degrees[1] < 6 or degrees != sorted(degrees)):
                        raise ValueError(f'malformed degree list {degrees}')
                    dimension = sum(d*d for d in degrees)
                    square = isqrt(dimension)**2 == dimension
                    if (square != (binary.name == 'square') or
                            any(dimension % (d*d) for d in set(degrees)) or
                            gcd(*degrees[1:]) != 1):
                        raise ValueError(f'inconsistent dimension list {degrees}')
                    reason = exclusion_reason(degrees, dimension=dimension)
                    if reason is None:
                        raise ValueError(f'surviving rank-{rank} profile {degrees}')
                    reasons[reason] += 1
            except BaseException:
                proc.kill()
                raise
            status = proc.wait()
        diagnostics.seek(0)
        receipt_text = diagnostics.read().strip()
    receipt = re.fullmatch(r'rank (\d+), (?:fours (-?\d+), )?nodes (\d+), '
                           r'(?:profiles|candidates) (\d+)', receipt_text)
    if (status or not receipt or int(receipt[1]) != rank or
            int(receipt[4]) != sum(reasons.values()) or
            (binary.name == 'square' and int(receipt[2]) !=
             (-1 if shard is None else shard))):
        raise RuntimeError(f'incomplete generation for {command}: {receipt_text}')
    return rank, binary.name, reasons, int(receipt[3])


def verify_small_ranks(ranks, jobs):
    totals, reasons, nodes = Counter(), Counter(), Counter()
    with tempfile.TemporaryDirectory() as build:
        tasks = []
        for case in ('square', 'nonsquare'):
            binary = Path(build)/case
            subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O3',
                            str(HERE/(case+'.cpp')), '-o', str(binary)], check=True)
            for rank in reversed(ranks):
                shards = reversed(range(14)) if case == 'square' and rank == 30 else [None]
                tasks.extend((binary, rank, shard) for shard in shards)
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            pending = [pool.submit(screen_small_job, task) for task in tasks]
            for result in as_completed(pending):
                rank, case, counts, node_count = result.result()
                totals[rank, case] += sum(counts.values())
                reasons.update(counts)
                nodes[case] += node_count

    expected = {}
    for line in (HERE/'counts.tsv').read_text().splitlines()[1:]:
        rank, square, nonsquare = map(int, line.split())
        expected[rank] = (square, nonsquare)
    print('rank     square     nonsquare    survivors', flush=True)
    for rank in ranks:
        square, nonsquare = totals[rank, 'square'], totals[rank, 'nonsquare']
        if (square, nonsquare) != expected[rank]:
            raise ValueError(f'rank {rank}: unexpected counts {square}, {nonsquare}')
        print(f'{rank:4} {square:10} {nonsquare:13} {0:12}', flush=True)
    print('First exclusions:', dict(sorted(reasons.items())), flush=True)
    print('Nodes:', dict(nodes), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument('--rank', type=int, choices=range(2, 32),
                        help='select one rank; the default checks the computational ranks 9–31')
    parser.add_argument('--work', type=Path, default=HERE/'rank31-results-fast',
                        help='resumable rank-31 result directory')
    parser.add_argument('--limit', type=int, help='run at most this many unfinished rank-31 tasks')
    parser.add_argument('--case', choices=('all', 'square', 'nonsquare'), default='all',
                        help='select a rank-31 dimension class')
    args = parser.parse_args()
    if args.jobs < 1 or (args.limit is not None and args.limit < 0):
        parser.error('--jobs must be positive and --limit nonnegative')
    if args.rank not in (None, 31) and (args.limit is not None or args.case != 'all'):
        parser.error('--limit and --case apply to rank 31')
    check_six()
    if args.rank is not None and args.rank < 9:
        print(f'Rank {args.rank} is excluded analytically in the manuscript; '
              'no degree enumeration is required.')
        return
    if args.rank != 31:
        verify_small_ranks(list(range(9, 31)) if args.rank is None else [args.rank],
                           args.jobs)
    if args.rank is None or args.rank == 31:
        command = [sys.executable, str(HERE/'run_rank31.py'), '--jobs', str(args.jobs),
                   '--work', str(args.work), '--case', args.case]
        if args.limit is not None:
            command.extend(('--limit', str(args.limit)))
        subprocess.run(command, check=True)
        work = args.work.resolve()
        with (work/'manifest.json').open() as file:
            manifest = json.load(file)
        if any(not (work/f'job-{index:05d}.jsonl').exists()
               for index in range(len(manifest['tasks']))):
            print('Rank 31 remains incomplete; continue with the same command and work directory.')
            return
        with (work/'summary.json').open() as file:
            summary = json.load(file)
        reference = {int(rank): (int(square), int(nonsquare))
                     for rank, square, nonsquare in
                     (line.split() for line in (HERE/'counts.tsv').read_text().splitlines()[1:])}
        expected = dict(zip(('square', 'nonsquare'), reference[31]))
        if (summary['rank'] != 31 or summary['complete'] is not True or
                len(manifest['tasks']) != 1167 or
                any(summary['by_case'][case]['profiles'] != count or
                    summary['by_case'][case]['survivors'] != 0
                    for case, count in expected.items()) or
                sum(v['tasks'] for v in summary['by_case'].values()) != 1167 or
                sum(summary['first_exclusions'].values()) != sum(expected.values()) or
                (work/'survivors.jsonl').read_text().strip()):
            raise ValueError('rank-31 result disagrees with counts.tsv or has survivors')
        print('Rank 31 verified against counts.tsv with zero survivors.')


if __name__ == '__main__':
    main()
