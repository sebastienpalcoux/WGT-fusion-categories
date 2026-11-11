#!/usr/bin/env python3
"""Resume the family-independent degree search at rank exactly 31 (Linux)."""
import argparse
from collections import Counter
import fcntl
import hashlib
import json
from math import gcd, isqrt
from multiprocessing import get_context
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time
import sys

from nonsquare_exact import InvalidPrefix, PRIMES, complete as nonsquare_complete, complete_prefix as nonsquare_prefix
from screen import exclusion_reason
from square_exact import complete as square_complete

HERE = Path(__file__).resolve().parent
SOURCES = ('run_rank31.py', 'square31.cpp', 'square_exact.py',
           'nonsquare31.cpp', 'nonsquare_exact.py', 'screen.py', 'trace.py', 'multiplication.py')


def packed(value):
    return json.dumps(value, separators=(',', ':'), sort_keys=True)


def job_list(square):
    listed = subprocess.run([square, '--tasks-cost'], capture_output=True, text=True, check=True)
    weighted = sorted(json.loads(line) for line in listed.stdout.splitlines())
    if listed.stderr.strip() != f'TASKS {len(weighted)}' or not weighted:
        raise ValueError('invalid square task manifest')
    previous = None
    for prefix, weight in weighted:
        if (not prefix or prefix[0] not in (4, 5) or prefix != sorted(prefix)
                or type(weight) is not int or weight < 1
                or (previous is not None and prefix[:len(previous)] == previous)):
            raise ValueError('invalid or overlapping square prefixes')
        previous = prefix
    tasks = []
    for q in range(2, 31):
        if any(q % (p*p) == 0 for p in range(2, isqrt(q)+1)):
            continue
        for b in range(1, isqrt(30//q)+1):
            f = q*b*b
            if f < 8 or (b == 1 and all(q%p for p in range(2, isqrt(q)+1))):
                continue
            multiplicities = range(1, f) if (q, b) in ((2, 2), (3, 2)) else (0,)
            tasks.extend(['nonsquare', [q, b, m]] for m in multiplicities)
    # Leaf counts affect scheduling only. Start expensive jobs on both sides early.
    light = [task for task in tasks if task[1][2] == 0]
    heavy = sorted((task for task in tasks if task[1][2]),
                   key=lambda task: (abs(task[1][2]-(task[1][0]*4-2)), task[1]))
    square_tasks = [['square', prefix] for prefix, _ in sorted(weighted, key=lambda v: (-v[1], v[0]))]
    smoke = [task for task in square_tasks if task[1] in ([5], [4, 5], [4]*8+[5])]
    rest = [task for task in square_tasks if task not in smoke]
    scheduled = light + smoke
    for i, task in enumerate(heavy):
        scheduled.append(task)
        if i < len(rest):
            scheduled.append(rest[i])
    return scheduled + rest[len(heavy):]


class InvalidProfile(ValueError):
    pass


def checked(row, case, key):
    if (not isinstance(row, list) or len(row) != 31 or any(type(d) is not int for d in row)
            or row != sorted(row) or row[0] != 1 or row[1] < 6 or gcd(*row[1:]) != 1):
        raise InvalidProfile(f'malformed degree list: {row}')
    q, first = (1, key[0]) if case == 'square' else key[:2]
    M = first*row[-1]
    D = q*M*M
    if D != sum(d*d for d in row) or any(M % d for d in set(row)):
        raise InvalidProfile('wrong square class or nonintegral denominator')
    if case == 'square':
        matches = all(d*b == M for d, b in zip(reversed(row), key))
    else:
        matches = not key[2] or row.count(row[-1]) == key[2]
    if not matches:
        raise InvalidProfile('degree list in wrong task')
    return exclusion_reason(row, dimension=D)


class GenerationFailure(Exception):
    """The native process did not certify a complete stream."""

    def __init__(self, message, *, retry=True):
        super().__init__(message)
        self.retry = retry


def run_native(exe, case, key, accept):
    args = ['--prefix', ','.join(map(str, key))] if case == 'square' else ['31', *map(str, key)]
    native_count = handoffs = 0
    stream_error = None
    with tempfile.TemporaryFile(mode='w+', encoding='utf-8', errors='replace') as err:
        with subprocess.Popen([exe, *args], stdout=subprocess.PIPE, stderr=err,
                              encoding='ascii', errors='surrogateescape') as child:
            try:
                for line_number, line in enumerate(child.stdout, start=1):
                    if stream_error is not None:
                        continue
                    try:
                        if not line.endswith('\n'):
                            raise ValueError('truncated final record')
                        item = json.loads(line)
                    except ValueError as exc:
                        # Even on abort the pipe may contain a fragment of the last row.
                        # Read to EOF so that the real exit status and stderr are available.
                        stream_error = f'output record {line_number}: {exc}'
                        continue
                    if isinstance(item, dict):
                        prefix = item.get('exact')
                        expected = key if case == 'square' else [key[1]]
                        if (not isinstance(prefix, list) or prefix[:len(expected)] != expected
                                or len(prefix) > 28 or any(type(b) is not int for b in prefix)
                                or prefix != sorted(prefix)):
                            stream_error = f'invalid continuation at output record {line_number}'
                            continue
                        handoffs += 1
                        if case == 'square':
                            rows = square_complete(prefix)
                        else:
                            primes = item.get('primes')
                            if (not isinstance(primes, list) or len(primes) < 3
                                    or any(type(p) is not int or p not in PRIMES for p in primes)
                                    or primes != sorted(set(primes))):
                                stream_error = f'invalid prime support at output record {line_number}'
                                continue
                            rows = nonsquare_prefix(*key, 31, primes, prefix)
                        try:
                            for row in rows:
                                accept(row)
                        except InvalidPrefix as exc:
                            stream_error = f'invalid continuation at output record {line_number}: {exc}'
                    else:
                        try:
                            accept(item)
                        except InvalidProfile as exc:
                            stream_error = f'output record {line_number}: {exc}'
                            continue
                        native_count += 1
                code = child.wait()
            finally:
                if child.poll() is None:
                    child.kill()
        err.seek(0)
        diag = err.read().strip()
    if code in (-signal.SIGINT, 128+signal.SIGINT):
        raise KeyboardInterrupt
    if case == 'square':
        pattern = rf'DONE nodes=(\d+) profiles={native_count} exact={handoffs}'
    else:
        q, b, m = key
        pattern = rf'rank 31, q {q}, B {b}, multiplicity {m}, nodes (\d+), candidates {native_count}, exact {handoffs}'
    match = re.fullmatch(pattern, diag)
    if code or stream_error or not match:
        detail = ['native exit '+str(code)]
        if stream_error:
            detail.append(stream_error)
        if diag:
            detail.append(diag[-800:])
        else:
            detail.append('missing completion receipt')
        overflow = 'overflow_error' in diag or 'OVERFLOW' in diag
        raise GenerationFailure(', '.join(detail), retry=not overflow)
    return int(match[1]), handoffs


def checkpoint_path(work, i):
    return work / f'job-{i:05d}.jsonl'


def worker(task):
    i, (case, key), build, work = task
    path = checkpoint_path(work, i)
    tmp = path.with_suffix('.part')
    reasons = Counter()
    survivors = 0
    exact_replay = False
    started = time.monotonic()
    with tmp.open('w+') as out:
        def accept(row):
            nonlocal survivors
            reason = checked(row, case, key)
            if reason is None:
                out.write(packed(row)+'\n')
                survivors += 1
            else:
                reasons[reason] += 1

        # A failed native pass is provisional; only a complete pass may commit.
        for attempt in (1, 2):
            try:
                nodes, handoffs = run_native(build / case, case, key, accept)
                break
            except GenerationFailure as exc:
                out.seek(0)
                out.truncate()
                reasons.clear()
                survivors = 0
                retry = attempt == 1 and exc.retry
                action = 'retrying native generator' if retry else 'replaying with exact integers'
                print(f'Task {i}, {case} {key}: {exc}; {action}.', file=sys.stderr, flush=True)
                if retry:
                    continue
                exact_replay = True
                rows = square_complete(key) if case == 'square' else nonsquare_complete(*key)
                for row in rows:
                    accept(row)
                nodes, handoffs = None, 0
                break
        receipt = dict(job=i, case=case, key=key, nodes=nodes, handoffs=handoffs,
                       exact_replay=exact_replay, excluded=dict(reasons), survivors=survivors,
                       profiles=sum(reasons.values())+survivors)
        out.write(packed(receipt)+'\n')
    tmp.replace(path)
    return i, receipt['profiles'], survivors, time.monotonic()-started


def interrupt_worker(signum, frame):
    raise SystemExit(128+signum)


def worker_signals():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, interrupt_worker)


def summarize(work, tasks):
    totals = {case: Counter() for case in ('square', 'nonsquare')}
    reasons = Counter()
    tmp = work/'survivors.part'
    with tmp.open('w') as out:
        for i, (case, key) in enumerate(tasks):
            n = 0
            receipt = None
            with checkpoint_path(work, i).open() as data:
                for line in data:
                    item = json.loads(line)
                    if receipt is not None:
                        raise ValueError(f'extra data after receipt {i}')
                    if isinstance(item, list):
                        if checked(item, case, key) is not None:
                            raise ValueError(f'invalid survivor in task {i}')
                        n += 1
                        out.write(packed(item)+'\n')
                    else:
                        receipt = item
            if (not receipt or (receipt['job'], receipt['case'], receipt['key']) != (i, case, key)
                    or receipt['survivors'] != n or
                    receipt['profiles'] != n+sum(receipt['excluded'].values())):
                raise ValueError(f'inconsistent receipt {i}')
            totals[case].update(tasks=1, profiles=receipt['profiles'], survivors=n)
            reasons.update(receipt['excluded'])
    tmp.replace(work/'survivors.jsonl')
    result = dict(rank=31, complete=True, by_case=totals, first_exclusions=reasons)
    tmp = work/'summary.part'
    tmp.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    tmp.replace(work/'summary.json')
    print(json.dumps(result, indent=2), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--jobs', type=int, default=4, help='number of parallel workers')
    ap.add_argument('--work', type=Path, default=HERE/'rank31-results-fast')
    ap.add_argument('--limit', type=int, help='run at most this many unfinished tasks')
    ap.add_argument('--case', choices=['all', 'square', 'nonsquare'], default='all')
    args = ap.parse_args()
    if args.jobs < 1 or (args.limit is not None and args.limit < 0):
        ap.error('--jobs must be positive and --limit nonnegative')
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    with (work/'run.lock').open('a') as lock, tempfile.TemporaryDirectory() as td:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            ap.exit(1, 'Another run already owns this work directory.\n')
        build = Path(td)
        for case in ('square', 'nonsquare'):
            subprocess.run([os.getenv('CXX', 'c++'), '-O3', '-march=native', '-std=c++17',
                            HERE/(case+'31.cpp'), '-o', build/case], check=True)
        tasks = job_list(build/'square')
        hashes = {s: hashlib.sha256((HERE/s).read_bytes()).hexdigest() for s in SOURCES}
        spec = packed(dict(rank=31, sha256=hashes, tasks=tasks))+'\n'
        manifest = work/'manifest.json'
        if manifest.exists() and manifest.read_text() != spec:
            ap.exit(1, 'Sources or tasks have changed. Supply a fresh --work directory.\n')
        if not manifest.exists():
            if any(work.glob('job-*.jsonl')):
                ap.exit(1, 'Missing manifest for existing checkpoints. Use a fresh --work directory.\n')
            tmp = work/'manifest.part'
            tmp.write_text(spec)
            tmp.replace(manifest)
        missing = [i for i in range(len(tasks)) if not checkpoint_path(work, i).is_file()]
        done = len(tasks)-len(missing)
        selected = [i for i in missing if args.case == 'all' or tasks[i][0] == args.case]
        if args.limit is not None:
            selected = selected[:args.limit]
        print(f'{done}/{len(tasks)} tasks already complete. Running {len(selected)}.', flush=True)
        if selected:
            # On Linux this command starts its workers before creating any threads.
            with get_context("fork").Pool(min(args.jobs, len(selected)), initializer=worker_signals) as pool:
                jobs = ((i, tasks[i], build, work) for i in selected)
                for i, rows, kept, seconds in pool.imap_unordered(worker, jobs):
                    done += 1
                    print(f'[{done}/{len(tasks)}] task {i} complete: {rows} profiles, '
                          f'{kept} survivors ({seconds:.1f}s)', flush=True)
        if done == len(tasks):
            summarize(work, tasks)
        else:
            print(f'Partial run; {len(tasks)-done} tasks remain. Checkpoints are in {work}.')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit('Interrupted. Rerun the command to resume completed tasks.')
