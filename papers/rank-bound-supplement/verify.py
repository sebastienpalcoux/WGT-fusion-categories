#!/usr/bin/env python3
"""Regenerate the degree profiles of a non-pointed simple integral modular category."""
import argparse
from collections import Counter
import fcntl
import gzip
import hashlib
import json
from math import gcd, isqrt
from multiprocessing import get_context
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile

from nonsquare_exact import complete as nonsquare_complete, complete_prefix
from primes import primes_up_to
from square_exact import complete as square_complete
from screen import exclusion_reason

HERE = Path(__file__).resolve().parent
RANKS = range(9, 32)
CASES = ('square', 'nonsquare')
FAMILIES = False
FAMILY_CERTIFICATE = None
FAMILY_SOURCES = ('catalogue.py', 'families.py', 'three.py', 'many_prime_families.py', 'support_certificates.py')
OBSTRUCTIONS = {'family', 'degree', 'crossed', 'activity', 'balance', 'orbit', 'dyadic', 'multiplication'}
# Capacities of the native position and prime-support masks, respectively.
NATIVE_MAX_RANK = {'square': 64, 'nonsquare': 136}
SOURCES = ('verify.py', 'square.cpp', 'nonsquare.cpp', 'square_exact.py',
           'nonsquare_exact.py', 'primes.py', 'screen.py', 'trace.py', 'multiplication.py')


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def atomic(path, data):
    tmp = path.with_suffix(path.suffix + '.part')
    tmp.write_bytes(data)
    tmp.replace(path)


def square_partition(square, rank):
    if rank > NATIVE_MAX_RANK['square']:
        return [([b], 1) for b in range(4, isqrt(rank-1)+1)]
    result = subprocess.run([square, '--rank', str(rank), '--tasks-cost', *(['--families'] if FAMILIES else [])],
                            capture_output=True, text=True, check=True)
    prefixes = sorted(json.loads(s) for s in result.stdout.splitlines())
    if result.stderr.strip() != f'TASKS {len(prefixes)}':
        raise ValueError('incomplete square partition')
    previous = None
    for prefix, weight in prefixes:
        if (not prefix or any(type(b) is not int for b in prefix) or prefix[0] < 4
                or prefix[0]**2 >= rank or prefix != sorted(prefix) or len(prefix) > rank-3
                or type(weight) is not int or weight < 1
                or (previous is not None and prefix[:len(previous)] == previous)):
            raise ValueError('invalid or overlapping square prefixes')
        previous = prefix
    return sorted(prefixes, key=lambda k: -k[1])


def task_list(square, rank):
    tasks = [[rank, 'square', key] for key, weight in square_partition(square, rank)]
    for q in range(2, rank):
        if any(q % (p*p) == 0 for p in range(2, isqrt(q)+1)):
            continue
        for b in range(1, isqrt((rank-1)//q)+1):
            f = q*b*b
            if f < 8 or (b == 1 and all(q % p for p in range(2, isqrt(q)+1))):
                continue
            # Split the two most expensive roots by maximal-degree multiplicity.
            for m in range(1, f) if (q, b) in ((2, 2), (3, 2)) else (0,):
                tasks.append([rank, 'nonsquare', [q, b, m]])
    return tasks


def checked(row, task):
    rank, case, key = task
    if (not isinstance(row, list) or len(row) != rank or any(type(d) is not int for d in row)
            or row != sorted(row) or row[0] != 1 or row[1] < 6 or gcd(*row[1:]) != 1):
        raise ValueError(f'malformed degree list in {task}: {row}')
    q, b = (1, key[0]) if case == 'square' else key[:2]
    M = b*row[-1]
    D = q*M*M
    if D != sum(d*d for d in row) or any(M % d for d in set(row)):
        raise ValueError('wrong square class or nonintegral denominator')
    matches = (all(d*z == M for d, z in zip(reversed(row), key)) if case == 'square'
               else not key[2] or row.count(row[-1]) == key[2])
    if not matches:
        raise ValueError('degree list outside its partition')
    if FAMILIES and FAMILY_CERTIFICATE(D):
        return 'family'
    return exclusion_reason(row, dimension=D)


def run_native(exe, task, accept):
    rank, case, key = task
    args = (['--rank', str(rank), '--prefix', ','.join(map(str, key))]
            if case == 'square' else [str(rank), *map(str, key)])
    if FAMILIES:
        args.append('--families')
    native_count = handoffs = 0
    with tempfile.TemporaryFile(mode='w+', encoding='utf-8', errors='replace') as err:
        with subprocess.Popen([exe, *args], stdout=subprocess.PIPE, stderr=err,
                              encoding='ascii', errors='surrogateescape') as proc:
            try:
                for line in proc.stdout:
                    if not line.endswith('\n'):
                        raise ValueError('truncated output')
                    item = json.loads(line)
                    if isinstance(item, list):
                        accept(item)
                        native_count += 1
                        continue
                    prefix = item['exact']
                    expected = key if case == 'square' else [key[1]]
                    if (not isinstance(prefix, list) or prefix[:len(expected)] != expected or
                            len(prefix) > rank-3 or prefix != sorted(prefix) or
                            any(type(b) is not int for b in prefix)):
                        raise ValueError('invalid exact-continuation prefix')
                    handoffs += 1
                    if case == 'square':
                        rows = square_complete(prefix, rank=rank, families=FAMILIES)
                    else:
                        primes = item['primes']
                        if (not isinstance(primes, list) or len(primes) < 3 or
                                primes != sorted(set(primes)) or any(p not in primes_up_to(rank) for p in primes)):
                            raise ValueError('invalid continuation prime support')
                        rows = complete_prefix(*key, rank, primes, prefix, families=FAMILIES)
                    for row in rows:
                        accept(row)
                status = proc.wait()
            finally:
                if proc.poll() is None:
                    proc.kill()
        err.seek(0)
        detail = err.read().strip()
    if status in (-signal.SIGINT, 128+signal.SIGINT):
        raise KeyboardInterrupt
    if case == 'square':
        receipt = rf'DONE nodes=(\d+) profiles={native_count} exact={handoffs}'
    else:
        q, b, m = key
        receipt = rf'rank {rank}, q {q}, B {b}, multiplicity {m}, nodes (\d+), candidates {native_count}, exact {handoffs}'
    match = re.fullmatch(receipt, detail)
    if status or not match:
        raise ValueError(f'native status {status}, receipt: {detail[-500:]}')
    return int(match[1]), handoffs


def worker(job):
    index, task, build, work = job
    rank, case, key = task
    reasons = Counter()
    survivors = []

    def accept(row):
        why = checked(row, task)
        if why is None:
            survivors.append(row)
        else:
            reasons[why] += 1

    nodes, handoffs = None, 0
    if rank <= NATIVE_MAX_RANK[case]:
        try:
            nodes, handoffs = run_native(build/case, task, accept)
        except (ValueError, KeyError, TypeError) as exc:
            print(f'{task}: {exc}; replaying this task with Python integers', file=sys.stderr, flush=True)
            reasons.clear()
            survivors.clear()
    if nodes is None:
        rows = (square_complete(key, rank=rank, families=FAMILIES) if case == 'square'
                else nonsquare_complete(*key, rank=rank, families=FAMILIES))
        for row in rows:
            accept(row)
    result = dict(task=task, nodes=nodes, handoffs=handoffs,
                  excluded=reasons, survivors=survivors, profiles=sum(reasons.values())+len(survivors))
    atomic(work/f'{index:05d}.json', encode(result).encode())
    return index, result['profiles'], len(survivors)


def summarize(ranks, tasks, receipts):
    if len(receipts) != len(tasks):
        raise ValueError('wrong number of task receipts')
    totals = {r: {k: 0 for k in (*CASES, 'survivors')} for r in ranks}
    for task, receipt in zip(tasks, receipts):
        excluded = receipt['excluded']
        survivors = receipt['survivors']
        if (task != receipt['task'] or not isinstance(survivors, list)
                or not isinstance(excluded, dict) or not set(excluded) <= OBSTRUCTIONS
                or any(type(v) is not int or v < 0 for v in excluded.values())
                or type(receipt['profiles']) is not int
                or receipt['profiles'] != len(survivors) + sum(excluded.values())
                or type(receipt['handoffs']) is not int or receipt['handoffs'] < 0
                or (receipt['nodes'] is None and receipt['handoffs'])
                or (receipt['nodes'] is not None and
                    (type(receipt['nodes']) is not int or receipt['nodes'] < 0))):
            raise ValueError(f'inconsistent receipt in {task}')
        if any(checked(row, task) is not None for row in survivors):
            raise ValueError(f'invalid survivor in {task}')
        rank, case, key = task
        totals[rank][case] += receipt['profiles']
        totals[rank]['survivors'] += len(survivors)
    return [dict(rank=rank, **counts) for rank, counts in totals.items()]


def display(summary):
    print('rank    square    nonsquare    survivors')
    for row in summary:
        print('{rank:4} {square:9} {nonsquare:12} {survivors:12}'.format(**row))
    if any(row['survivors'] for row in summary):
        message = 'Surviving degree profiles are recorded in the certificate.'
        if any(row['rank'] in RANKS and row['survivors'] for row in summary):
            raise SystemExit(message)
        print(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rank', type=int, help='select one positive rank (default: ranks 9 through 31)')
    parser.add_argument('--families', action='store_true', help='also use the full article dimension families (requires its repository checkout)')
    parser.add_argument('--jobs', type=int, default=4, help='parallel workers')
    parser.add_argument('--work', type=Path, default=Path('rank-results'), help='checkpoint directory')
    parser.add_argument('--limit', type=int, help='run at most this many unfinished tasks')
    parser.add_argument('--check', type=Path, metavar='CERTIFICATE', help='check archived receipts')
    args = parser.parse_args()
    global FAMILIES, FAMILY_CERTIFICATE
    FAMILIES = args.families
    family_root = HERE.parents[1]/'rank'/'code'
    if FAMILIES:
        if not all((family_root/name).is_file() for name in FAMILY_SOURCES):
            parser.error('--families requires the full repository, including rank/code')
        sys.path.insert(0, str(family_root))
        from catalogue import certificate
        FAMILY_CERTIFICATE = certificate
    if args.jobs < 1 or (args.limit is not None and args.limit < 0):
        parser.error('--jobs must be positive and --limit nonnegative')
    if args.rank is not None and args.rank < 1:
        parser.error('--rank must be positive')
    if args.rank is not None and args.rank < min(RANKS):
        if args.check:
            parser.error('ranks 1 through 8 are excluded analytically and have no computation certificate')
        print(f'Rank {args.rank} is excluded analytically.')
        return
    if args.rank == 9:
        print('Rank 9 is also excluded analytically; running the uniform verification.', flush=True)
    if args.rank is not None and args.rank > max(RANKS):
        print(f'Rank {args.rank} is exploratory and lies beyond the published verification.', flush=True)
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)  # Degree lists and receipts can contain arbitrarily large integers.
    ranks = [args.rank] if args.rank is not None else list(RANKS)
    # The denominator tree and terminal tests recurse over labels and prime supports.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), min(2**31-1, 8*max(ranks)+1024)))
    with tempfile.TemporaryDirectory() as td:
        build = Path(td)
        for case in CASES:
            if min(ranks) > NATIVE_MAX_RANK[case]:
                continue
            subprocess.run([os.getenv('CXX', 'c++'), '-O3', '-march=native', '-std=c++17',
                            HERE/(case+'.cpp'), '-o', build/case], check=True)
        tasks = [task for rank in reversed(ranks) for task in task_list(build/'square', rank)]
        hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in SOURCES}
        if FAMILIES:
            for name in FAMILY_SOURCES:
                hashes['../../rank/code/'+name] = hashlib.sha256((family_root/name).read_bytes()).hexdigest()
        spec = dict(ranks=ranks, sha256=hashes, tasks=tasks, families=FAMILIES)
        if args.check:
            with gzip.open(args.check, 'rt') as f:
                cert = json.load(f)
            if cert['manifest'] != spec:
                raise ValueError('certificate source hashes or task partition differ')
            summary = summarize(ranks, tasks, cert['receipts'])
            if cert['summary'] != summary:
                raise ValueError('archive summary disagrees with receipts')
            display(summary)
            return
        work = args.work.resolve()
        work.mkdir(parents=True, exist_ok=True)
        with (work/'run.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            manifest = work/'manifest.json'
            if manifest.exists():
                if json.loads(manifest.read_text()) != spec:
                    raise ValueError('sources or selected ranks changed; use a fresh work directory')
            elif any(work.glob('[0-9]*.json')):
                raise ValueError('checkpoint directory has no manifest')
            else:
                atomic(manifest, encode(spec).encode())
            missing = [i for i in range(len(tasks)) if not (work/f'{i:05d}.json').is_file()]
            selected = missing[:args.limit]
            done = len(tasks)-len(missing)
            print(f'{done}/{len(tasks)} tasks complete. Running {len(selected)}.', flush=True)
            if selected:
                jobs = ((i, tasks[i], build, work) for i in selected)
                with get_context('fork').Pool(min(args.jobs, len(selected))) as pool:
                    for i, profiles, kept in pool.imap_unordered(worker, jobs):
                        done += 1
                        print(f'[{done}/{len(tasks)}] rank {tasks[i][0]} task {i}: {profiles} profiles, {kept} survivors', flush=True)
            if done < len(tasks):
                return
            receipts = [json.loads((work/f'{i:05d}.json').read_text()) for i in range(len(tasks))]
            summary = summarize(ranks, tasks, receipts)
            cert = dict(manifest=spec, receipts=receipts, summary=summary)
            atomic(work/'certificate.json.gz', gzip.compress(encode(cert).encode(), mtime=0))
            display(summary)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit('Interrupted. Rerun with the same work directory to resume.')
