#!/usr/bin/env python3
"""Fresh, source-bound degree enumeration for every rank from 1 through 44."""
import argparse
from collections import Counter
from contextlib import suppress
import ctypes
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
from math import isfinite
import multiprocessing as mp
from multiprocessing.connection import wait
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parent
for directory in (ROOT, ROOT/'engine', ROOT/'engine_src', ROOT/'phase-mass'):
    sys.path.insert(0, str(directory))
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
import square
import nonsquare
import terminal as checker
import preflight

RANK = 44
FORMAT = 'uniform-rank-census-v2'
SOURCE_DIRS = ('engine', 'engine_src', 'arithmetic', 'phase-mass', 'tests')


def now():
    return datetime.now(timezone.utc).isoformat(timespec='microseconds')


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode('utf-8')


def digest(value):
    return hashlib.sha256(encode(value)).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def save_json(path, data):
    path = Path(path)
    tmp = path.with_suffix('.tmp')
    with tmp.open('wb') as stream:
        stream.write(encode(data))
        stream.flush()
        os.fsync(stream.fileno())
    tmp.replace(path)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def sources():
    files = list(ROOT.glob('*.py'))
    for name in SOURCE_DIRS:
        files.extend(p for p in (ROOT/name).rglob('*')
                     if p.is_file() and p.suffix in ('.py', '.cpp', '.hpp', '.json'))
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(files)}


def build(directory):
    directory.mkdir(exist_ok=True)
    compiler = shutil.which(os.environ.get('CXX', 'g++'))
    if not compiler:
        raise RuntimeError('Install a C++17 compiler: sudo apt install g++')
    flags = ['-std=c++17', '-O3', '-march=native']
    tick = time.monotonic()
    for name, source in [('square', 'square44.cpp'), ('nonsquare', 'nonsquare.cpp')]:
        subprocess.run([compiler, *flags, '-I', str(ROOT/'engine'),
                        str(ROOT/'engine_src'/source), '-o', str(directory/name)], check=True)
    return dict(compiler=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                flags=flags, seconds=round(time.monotonic()-tick, 6),
                binaries={name: sha(directory/name) for name in ('square', 'nonsquare')})


def plan(directory):
    sq = square.tasks(RANK, directory)
    ns = nonsquare.tasks(RANK)
    # Submit expensive square prefixes early and advance nonsquare cases too.
    tasks = []
    for i, task in enumerate(ns):
        tasks.append(task)
        tasks.extend(sq[4*i:4*i+4])
    tasks.extend(sq[4*len(ns):])
    if len({digest(t) for t in tasks}) != len(tasks):
        raise ValueError('duplicate task in the regenerated partition')
    groups = {}
    for task in tasks:
        q, b = task['q'], task['B']
        prefix = task.get('prefix', [b])
        valid = ((b, q) in square.start_pairs(RANK)) or ((b, q) in nonsquare.start_pairs(RANK))
        if (not valid or not prefix or prefix[0] != b or prefix != sorted(prefix)
                or len(prefix) > RANK-3 or any(type(z) is not int or z < b for z in prefix)):
            raise ValueError('invalid task in the regenerated partition')
        if q > 1 and task['m'] not in nonsquare.multiplicities(RANK, b, q):
            raise ValueError('invalid top multiplicity')
        groups.setdefault((q, b, task.get('m')), []).append(prefix)
    for prefixes in groups.values():
        ordered = sorted(prefixes)
        if any(b[:len(a)] == a for a, b in zip(ordered, ordered[1:])):
            raise ValueError('overlapping task prefixes')
    return tasks


def check_manifest(spec):
    if (spec.get('format') != FORMAT or spec.get('rank') != RANK
            or digest({k: v for k, v in spec.items() if k != 'sha256'}) != spec.get('sha256')
            or digest(spec['tasks']) != spec.get('partition_sha256')):
        raise ValueError('damaged or incompatible work manifest')
    if spec['source_sha256'] != sources():
        raise ValueError('Source files changed. Use the original package to resume, or a fresh --work directory.')


def reason_for(row, task):
    return checker.reason_for(row, task, RANK,
                              adapter_failure=nonsquare.row_failure if task['q'] > 1 else None)


def command(directory, task):
    if task['q'] == 1:
        return square.command(directory, task, RANK)
    args = [str(directory/'nonsquare'), str(RANK), str(task['q']),
            str(task['B']), str(task['m'])]
    if 'prefix' in task:
        args += ['--prefix', ','.join(map(str, task['prefix']))]
    return args


def run_task(work, directory, index, task):
    started, tick = now(), time.monotonic()
    ledger_path = work/'rows'/f'{index:06d}.jsonl.gz'
    tmp = ledger_path.with_suffix('.tmp')
    native = handoffs = extra = count = 0
    excluded, survivors = Counter(), []
    records_hash = hashlib.sha256()
    exact_rows = square.exact_rows if task['q'] == 1 else nonsquare.exact_rows
    try:
        with tmp.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as ledger:
                def accept(row, backend):
                    nonlocal count
                    reason = reason_for(row, task)
                    if reason is None:
                        survivors.append(row)
                    elif reason in checker.OBSTRUCTIONS:
                        excluded[reason] += 1
                    else:
                        raise ValueError('unknown terminal obstruction: '+str(reason))
                    record = encode([backend, reason, row])+b'\n'
                    ledger.write(record)
                    records_hash.update(record)
                    count += 1

                with tempfile.TemporaryFile(mode='w+b') as err:
                    with subprocess.Popen(command(directory, task), stdout=subprocess.PIPE,
                                          stderr=err, bufsize=65536) as proc:
                        try:
                            for line in proc.stdout:
                                if not line.endswith(b'\n'):
                                    raise ValueError('truncated native output')
                                item = json.loads(line)
                                if isinstance(item, list):
                                    accept(item, 'native')
                                    native += 1
                                elif isinstance(item, dict):
                                    handoffs += 1
                                    # Exhaustion, including zero output, is required for completion.
                                    for row in exact_rows(item, task, RANK):
                                        accept(row, 'exact')
                                        extra += 1
                                else:
                                    raise ValueError('invalid native output object')
                            status = proc.wait()
                        finally:
                            if proc.poll() is None:
                                proc.kill()
                    err.seek(0)
                    diagnostic = err.read().decode('utf-8', errors='strict').strip()
                if task['q'] == 1:
                    pattern = rf'DONE nodes=(\d+) profiles={native} exact={handoffs}'
                else:
                    pattern = (rf'rank {RANK}, q {task["q"]}, B {task["B"]}, multiplicity {task["m"]}, '
                               rf'nodes (\d+), candidates {native}, exact {handoffs}')
                match = re.fullmatch(pattern, diagnostic)
                if status != 0 or match is None:
                    raise RuntimeError(f'native task incomplete (exit {status}): {diagnostic!r}')
            raw.flush()
            os.fsync(raw.fileno())
        if count != native+extra:
            raise ValueError('raw ledger count mismatch')
        tmp.replace(ledger_path)
        return dict(nodes=int(match.group(1)), native=native, handoffs=handoffs,
                    exact_rows=extra, profiles=count, excluded=dict(excluded), survivors=survivors,
                    started_utc=started, finished_utc=now(), seconds=round(time.monotonic()-tick, 6),
                    native_completion=diagnostic, native_exit=0, exact_continuations_exhausted=True,
                    raw_ledger=str(ledger_path.relative_to(work)), raw_ledger_sha256=sha(ledger_path),
                    raw_records_sha256=records_hash.hexdigest())
    finally:
        with suppress(FileNotFoundError):
            tmp.unlink()


def worker(channel, work, directory, coordinator):
    os.setsid()
    def stop_group(*_):
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        os.killpg(os.getpgrp(), signal.SIGTERM)
    signal.signal(signal.SIGTERM, stop_group)
    libc = ctypes.CDLL(None, use_errno=True)
    libc.prctl.argtypes = [ctypes.c_int]+[ctypes.c_ulong]*4
    if libc.prctl(1, signal.SIGTERM, 0, 0, 0) != 0:
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err))
    if os.getppid() != coordinator:
        stop_group()
    try:
        while True:
            job = channel.recv()
            if job is None:
                return
            index, task = job
            try:
                channel.send((index, run_task(work, directory, index, task), None))
            except Exception:
                channel.send((index, None, traceback.format_exc()))
                return
    finally:
        channel.close()


def receipt_path(work, index):
    return work/'receipts'/f'{index:06d}.json'


def check_receipt(rec, spec):
    i = rec['index']
    if (type(i) is not int or not 0 <= i < len(spec['tasks'])
            or rec['manifest_sha256'] != spec['sha256']
            or any(type(rec[k]) is not int or rec[k] < 0
                   for k in ('nodes', 'native', 'handoffs', 'exact_rows', 'profiles'))
            or type(rec['seconds']) not in (int, float) or not isfinite(rec['seconds'])
            or rec['seconds'] < 0 or rec['native_exit'] != 0
            or rec['exact_continuations_exhausted'] is not True
            or rec['native']+rec['exact_rows'] != rec['profiles']
            or (rec['exact_rows'] and not rec['handoffs'])
            or not isinstance(rec['survivors'], list) or not isinstance(rec['excluded'], dict)
            or not rec['excluded'].keys() <= checker.OBSTRUCTIONS
            or any(type(n) is not int or n < 0 for n in rec['excluded'].values())
            or sum(rec['excluded'].values())+len(rec['survivors']) != rec['profiles']
            or rec['raw_ledger'] != f'rows/{i:06d}.jsonl.gz'
            or not re.fullmatch('[0-9a-f]{64}', rec.get('raw_records_sha256', ''))):
        raise ValueError(f'invalid receipt for task {i!r}')
    task = spec['tasks'][i]
    expected = (f'DONE nodes={rec["nodes"]} profiles={rec["native"]} exact={rec["handoffs"]}'
                if task['q'] == 1 else
                f'rank {spec["rank"]}, q {task["q"]}, B {task["B"]}, multiplicity {task["m"]}, '
                f'nodes {rec["nodes"]}, candidates {rec["native"]}, exact {rec["handoffs"]}')
    if rec['native_completion'] != expected:
        raise ValueError('native completion receipt disagrees with task accounting')
    for row in rec['survivors']:
        checker.validate_row(row, task, spec["rank"])
    return i


def receipts(work, spec):
    done = {}
    for path in sorted((work/'receipts').glob('*.json')):
        packet = read_json(path)
        rec = packet['receipt']
        if packet.get('sha256') != digest(rec):
            raise ValueError('damaged task receipt: '+str(path))
        index = check_receipt(rec, spec)
        if path != receipt_path(work, index) or index in done:
            raise ValueError('duplicate or misnamed task receipt')
        if sha(work/rec['raw_ledger']) != rec['raw_ledger_sha256']:
            raise ValueError('raw ledger checksum mismatch')
        done[index] = rec
    return done


def compute(work, directory, spec, selected, jobs, done, running):
    it, procs, busy = iter(selected), [], {}
    ctx = mp.get_context('fork')
    def assign(channel):
        index = next(it, None)
        if index is None:
            channel.send(None)
        else:
            busy[channel] = index
            running[index] = time.monotonic()
            channel.send((index, spec['tasks'][index]))
    def progress():
        oldest = sorted(running, key=running.get)[:3]
        ages = ', '.join(f'#{i} {int(time.monotonic()-running[i])}s' for i in oldest)
        print(f'{len(done):,}/{len(spec["tasks"]):,} tasks complete; {len(running)} running'
              + ('; longest '+ages if ages else ''), flush=True)
    try:
        for _ in range(min(jobs, len(selected))):
            parent, child = ctx.Pipe()
            proc = ctx.Process(target=worker, args=(child, work, directory, os.getpid()))
            proc.start()
            child.close()
            procs.append((proc, parent))
            assign(parent)
        last = time.monotonic()
        while busy:
            for channel in wait(list(busy), timeout=2):
                index = busy.pop(channel)
                received, data, error = channel.recv()
                if received != index or error:
                    raise RuntimeError(f'task {index} {spec["tasks"][index]} failed\n{error}')
                rec = dict(index=index, manifest_sha256=spec['sha256'], **data)
                check_receipt(rec, spec)
                save_json(receipt_path(work, index), dict(receipt=rec, sha256=digest(rec)))
                done[index] = rec
                del running[index]
                assign(channel)
            # A forked worker can inherit another worker's pipe descriptor.
            # Detect a dead busy worker even if that inheritance delays EOF.
            for proc, channel in procs:
                if channel in busy and proc.exitcode is not None:
                    raise RuntimeError(f'worker for task {busy[channel]} exited unexpectedly '
                                       f'with status {proc.exitcode}; no completion is certified')
            if time.monotonic()-last >= 30:
                progress()
                last = time.monotonic()
        progress()
    finally:
        for proc, channel in procs:
            if running:
                with suppress(ProcessLookupError):
                    os.killpg(proc.pid, signal.SIGTERM)
            channel.close()
        for proc, _ in procs:
            proc.join(2)
            if proc.is_alive():
                with suppress(ProcessLookupError):
                    os.killpg(proc.pid, signal.SIGKILL)
                proc.kill()
                proc.join()


def verify_rows(work, spec, done, recheck=False):
    """Bind every row and its classification to a completed task receipt."""
    seen = set()
    for index, rec in sorted(done.items()):
        counts, reasons, retained = Counter(), Counter(), []
        records_hash = hashlib.sha256()
        with gzip.open(work/rec['raw_ledger'], 'rt') as stream:
            for line in stream:
                backend, reason, row = json.loads(line)
                if backend not in ('native', 'exact') or (reason is not None and reason not in checker.OBSTRUCTIONS):
                    raise ValueError('invalid raw row record')
                checker.validate_row(row, spec['tasks'][index], RANK)
                key = tuple(row)
                if key in seen:
                    raise ValueError('duplicate degree list in the native/exact task partition')
                seen.add(key)
                records_hash.update(encode([backend, reason, row])+b'\n')
                if recheck and reason_for(row, spec['tasks'][index]) != reason:
                    raise ValueError('terminal recheck differs from the recorded classification')
                counts[backend] += 1
                if reason is None:
                    retained.append(row)
                else:
                    reasons[reason] += 1
        if (counts['native'] != rec['native'] or counts['exact'] != rec['exact_rows']
                or sum(counts.values()) != rec['profiles'] or dict(reasons) != rec['excluded']
                or retained != rec['survivors'] or records_hash.hexdigest() != rec['raw_records_sha256']):
            raise ValueError('raw row ledger disagrees with task receipt')


def summary_text(summary):
    lines = [f'RANK {summary["rank"]} | {summary["status"].upper()} | '
             f'{summary["completed"]}/{summary["tasks"]} tasks complete', '',
             f'Raw degree profiles: {summary["profiles"]}',
             f'Arbitrary-precision handoffs fully completed: {summary["handoffs"]}',
             f'Distinct retained degree profiles: {summary["survivor_count"]}',
             f'Total compute-session wall seconds: {summary["compute_session_seconds"]}',
             f'Sum of completed worker-task seconds: {summary["worker_task_seconds"]}', '']
    if summary['status'] != 'complete':
        lines += ['This is a PARTIAL computation. Unfinished tasks carry no conclusion.',
                  'Resume the same rank selection with: python3 run.py --resume --jobs 8', '']
    elif summary['survivor_count']:
        lines += ['The retained degree profiles passed necessary conditions only.',
                  'Their presence is not a construction of modular data.', '']
    else:
        lines += ['Every task completed and no degree profile passed all implemented conditions.', '']
    lines += ['See survivors.jsonl for all retained profiles and manifest.json for the full task partition.',
              'The ZIP includes raw-census.jsonl.gz, every completed task receipt, fresh structural',
              'checks, source files, source hashes, timings, and any error report.']
    if summary['error']:
        lines += ['', summary['error']]
    return '\n'.join(lines)+'\n'


def report(work, spec, done, running, session, error=None, recheck=False):
    check_manifest(spec)
    structural = read_json(work/'preflight.json')
    if digest(structural) != spec['preflight_sha256'] or structural.get('status') != 'complete':
        raise ValueError('fresh structural check is missing or damaged')
    verify_rows(work, spec, done, recheck=recheck)
    history = read_json(work/'sessions.json')
    totals, reasons, unique = {}, Counter(), {}
    for index, rec in sorted(done.items()):
        task = spec['tasks'][index]
        reasons.update(rec['excluded'])
        for row in rec['survivors']:
            key = (task['q'], task['B'], tuple(row))
            entry = unique.setdefault(key, dict(rank=RANK, B=task['B'], q=task['q'],
                                               dimension=task['q']*(task['B']*row[-1])**2,
                                               degrees=row, tasks=[]))
            entry['tasks'].append(index)
    for b, q in sorted({(t['B'], t['q']) for t in spec['tasks']}):
        selected = [i for i in done if (spec['tasks'][i]['B'], spec['tasks'][i]['q']) == (b, q)]
        totals[f'{b},{q}'] = dict(B=b, q=q, completed=len(selected),
                                 tasks=sum((t['B'], t['q']) == (b, q) for t in spec['tasks']),
                                 **{k: sum(done[i][k] for i in selected)
                                    for k in ('nodes', 'native', 'exact_rows', 'profiles', 'handoffs')})
    summary = dict(format=FORMAT, rank=RANK,
                   status='complete' if len(done) == len(spec['tasks']) and error is None else 'partial',
                   tasks=len(spec['tasks']), completed=len(done),
                   missing_tasks=[i for i in range(len(spec['tasks'])) if i not in done],
                   profiles=sum(r['profiles'] for r in done.values()),
                   native=sum(r['native'] for r in done.values()),
                   exact_rows=sum(r['exact_rows'] for r in done.values()),
                   handoffs=sum(r['handoffs'] for r in done.values()),
                   survivor_count=len(unique), survivor_records=sum(len(r['survivors']) for r in done.values()),
                   excluded_by=dict(reasons), cases=totals, manifest_sha256=spec['sha256'],
                   source_sha256=digest(spec['source_sha256']),
                   interrupted_tasks=sorted(running), session=session, error=error,
                   compute_session_seconds=round(sum(s['seconds'] for s in history
                                                     if s['operation'] == 'compute'), 6),
                   worker_task_seconds=round(sum(r['seconds'] for r in done.values()), 6),
                   row_accounting_verified=True, all_terminal_rechecks=bool(recheck))
    census = work/'raw-census.jsonl.gz'
    with census.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as output:
            for index, rec in sorted(done.items()):
                with gzip.open(work/rec['raw_ledger'], 'rt') as stream:
                    for line in stream:
                        backend, reason, row = json.loads(line)
                        output.write(encode([RANK, index, backend, reason, row])+b'\n')
    summary['raw_census_sha256'] = sha(census)
    message = summary_text(summary)
    save_json(work/'summary.json', summary)
    (work/'summary.txt').write_text(message)
    print(f'Rank {RANK}: {summary["status"]}; '
          f'{summary["completed"]}/{summary["tasks"]} tasks; '
          f'{summary["profiles"]} degree profiles; '
          f'{summary["survivor_count"]} retained.', flush=True)
    return summary


def run_rank(args):
    global RANK
    RANK = args.rank
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    with (work/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        existing = (work/'manifest.json').exists()
        if existing and not (args.resume or args.report or args.verify):
            raise ValueError('This work directory already has a manifest. Use --resume, or --work with a NEW directory for a fresh computation.')
        if not existing and (args.report or args.verify):
            raise ValueError('No completed preparation in this work directory; first run the computation.')
        for name in ('receipts', 'rows'):
            (work/name).mkdir(exist_ok=True)
        if not existing and any((work/'receipts').glob('*.json')):
            raise ValueError('Receipts without a manifest cannot be reused. Choose a NEW --work directory.')
        available = len(os.sched_getaffinity(0))
        jobs = min(args.jobs, available)
        tick = time.monotonic()
        initial_sources = sources()
        session = dict(started_utc=now(), python=platform.python_version(),
                       system=platform.platform(), jobs=jobs, available_logical_cpus=available,
                       fresh=not existing, operation='verify' if args.verify else 'report' if args.report else 'compute')
        print(f'Rank {RANK} | {jobs} workers maximum | '+('resuming local work' if existing else 'fresh computation'), flush=True)
        if existing:
            spec = read_json(work/'manifest.json')
            check_manifest(spec)
        session['build'] = args.preparation['build']
        print('Regenerating the complete square and nonsquare task partition...', flush=True)
        tasks = plan(args.build_directory)
        if existing:
            if tasks != spec['tasks']:
                raise ValueError('regenerated task partition differs from the saved manifest')
        else:
            structural = args.preparation['structural']
            save_json(work/'preflight.json', structural)
            if sources() != initial_sources:
                raise ValueError('Sources changed during preparation. Restart with unchanged sources.')
            spec = dict(format=FORMAT, rank=RANK, created_utc=now(), source_sha256=initial_sources,
                        preflight_sha256=digest(structural), tasks=tasks, partition_sha256=digest(tasks),
                        scope=f'Exact rank {RANK}; general simple integral rotation-admissible modular data.',
                        initial_state='Fresh: no earlier census, result or task receipt was loaded.',
                        theoretical_reductions=['Odd dimension excluded below rank209.',
                                                'Largest-degree codegree8 requires rank at least46.',
                                                'Square and nonsquare support, packet and phase restrictions at the requested rank.',
                                                'Uniform actual/effective prime-square trace support and gcd obstruction.'],
                        planner=dict(square='complete native --tasks-cost 24 32',
                                     nonsquare='all admissible supports and top multiplicities; disjoint exact prefixes'),
                        receipt_semantics='A task completes only after exit0, matched native counts and full exhaustion of all exact continuations.')
            spec['sha256'] = digest(spec)
            save_json(work/'manifest.json', spec)
        done = receipts(work, spec)
        missing = [i for i in range(len(tasks)) if i not in done]
        print(f'{len(tasks):,} tasks: {sum(t["q"] == 1 for t in tasks):,} square, '
              f'{sum(t["q"] > 1 for t in tasks):,} nonsquare. {len(done):,} already complete.', flush=True)
        selected = missing if args.limit is None else missing[:args.limit]
        running, error, status = {}, None, 0
        session['completed_before'] = len(done)
        try:
            if not (args.report or args.verify):
                compute(work, args.build_directory, spec, selected, jobs, done, running)
        except KeyboardInterrupt:
            status = 130
            print('Interrupted. Saving completed tasks and writing a PARTIAL results ZIP...', flush=True)
            done = receipts(work, spec)
            for index in done:
                running.pop(index, None)
        except Exception:
            status, error = 1, traceback.format_exc()
        session.update(finished_utc=now(), seconds=round(time.monotonic()-tick, 6),
                       completed_after=len(done), exit_status=status)
        history_path = work/'sessions.json'
        history = read_json(history_path) if history_path.exists() else []
        history.append(session)
        save_json(history_path, history)
        print('Checking row accounting and packaging the results...', flush=True)
        report(work, spec, done, running, session, error, recheck=args.recheck)
        return status


def batch_report(work, ranks, session, error=None):
    """One compact result archive, containing source once and fresh rank records."""
    summaries = []
    for rank in ranks:
        path = work/f'rank-{rank:02d}'/'summary.json'
        if path.exists():
            item = read_json(path)
            if item['rank'] != rank or item['format'] != FORMAT:
                raise ValueError('rank summary disagrees with the requested batch')
            summaries.append(item)
    complete = (len(summaries) == len(ranks)
                and all(s['status'] == 'complete' for s in summaries)
                and error is None)
    summary = dict(format=FORMAT, requested_ranks=ranks,
                   status='complete' if complete else 'partial',
                   ranks=[dict(rank=s['rank'], status=s['status'], tasks=s['tasks'],
                               completed=s['completed'], profiles=s['profiles'],
                               square=sum(c['profiles'] for c in s['cases'].values() if c['q'] == 1),
                               nonsquare=sum(c['profiles'] for c in s['cases'].values() if c['q'] != 1),
                               survivors=s['survivor_count'],
                               compute_seconds=s['compute_session_seconds'],
                               worker_seconds=s['worker_task_seconds'],
                               raw_census_sha256=s['raw_census_sha256'],
                               manifest_sha256=s['manifest_sha256']) for s in summaries],
                   tasks=sum(s['tasks'] for s in summaries),
                   completed=sum(s['completed'] for s in summaries),
                   survivors=sum(s['survivor_count'] for s in summaries),
                   source_sha256=sources(), session=session, error=error,
                   scope='Nonpointed based-simple integral rotation-admissible modular data; necessary conditions only.',
                   computation='All requested rank partitions and degree profiles are generated from zero. No earlier census or receipts are loaded.')
    if complete and summary['survivors'] == 0 and ranks == list(range(1, max(ranks)+1)):
        summary['rank_lower_bound_in_stated_scope'] = max(ranks)+1
    save_json(work/'summary.json', summary)
    title = 'RANKS '+(','.join(map(str,ranks)) if len(ranks)<4 else f'{ranks[0]}–{ranks[-1]}')
    lines = [title+' | '+summary['status'].upper(), '',
             'Rank  Tasks complete/total  Square  Nonsquare  Retained  Wall seconds']
    for s in summary['ranks']:
        lines.append(f'{s["rank"]:4d}  {s["completed"]:7d}/{s["tasks"]:<7d}  '
                     f'{s["square"]:6d}  {s["nonsquare"]:9d}  {s["survivors"]:8d}  '
                     f'{s["compute_seconds"]:12.6f}')
    lines += ['', f'Total completed tasks: {summary["completed"]}/{summary["tasks"]}',
              f'Total retained profiles: {summary["survivors"]}',
              'Scope: nonpointed based-simple integral rotation-admissible modular data.',
              'Partial tasks and resource-limited necessary tests never certify exclusions.']
    if error:
        lines += ['', error]
    message = '\n'.join(lines)+'\n'
    (work/'summary.txt').write_text(message)
    path = work/'rank-results.zip'
    tmp = path.with_suffix('.tmp')
    with ZipFile(tmp,'w',compression=ZIP_DEFLATED,compresslevel=9) as archive:
        archive.write(work/'summary.json','summary.json')
        archive.write(work/'summary.txt','summary.txt')
        for name in ('batch.json','preparation.json'):
            if (work/name).exists():
                archive.write(work/name,name)
        for s in summaries:
            rank = s['rank']
            directory = work/f'rank-{rank:02d}'
            base = f'ranks/{rank:02d}/'
            for name in ('summary.json','summary.txt','manifest.json','sessions.json','raw-census.jsonl.gz'):
                archive.write(directory/name,base+name)
            spec = read_json(directory/'manifest.json')
            done = receipts(directory,spec)
            with archive.open(base+'receipts.jsonl','w',force_zip64=True) as stream:
                for index in sorted(done):
                    stream.write(encode(done[index])+b'\n')
            with archive.open(base+'survivors.jsonl','w',force_zip64=True) as stream:
                for index in sorted(done):
                    for row in done[index]['survivors']:
                        stream.write(encode(dict(task=index,degrees=row))+b'\n')
        for name in summary['source_sha256']:
            archive.write(ROOT/name,'source/'+name)
        for name in ('README.txt','Prime-square-criterion.tex'):
            if (ROOT/name).exists():
                archive.write(ROOT/name,'source/'+name)
    tmp.replace(path)
    print('\n'+message,flush=True)
    print(f'Results ZIP: {path}\nZIP size: {path.stat().st_size:,} bytes',flush=True)
    return summary


def run(args):
    work = args.work.resolve()
    ranks = [args.rank] if args.rank is not None else list(range(1,args.through+1))
    work.mkdir(parents=True,exist_ok=True)
    with (work/'batch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        existing = (work/'batch.json').exists()
        if existing and not (args.resume or args.report or args.verify):
            raise ValueError('This work directory already contains a batch. Use --resume or a new --work directory.')
        if not existing and (args.report or args.verify):
            raise ValueError('No prepared batch exists in this work directory.')
        hashes = sources()
        if existing:
            batch = read_json(work/'batch.json')
            if batch['format'] != FORMAT or batch['ranks'] != ranks or batch['source_sha256'] != hashes:
                raise ValueError('Source or rank selection differs from the saved batch; use an unchanged package or a fresh work directory.')
        else:
            batch = dict(format=FORMAT,ranks=ranks,source_sha256=hashes,created_utc=now(),
                         initial_state='Fresh computation; no historical results loaded.')
            save_json(work/'batch.json',batch)
        tick = time.monotonic()
        session = dict(started_utc=now(),fresh=not existing,requested_ranks=ranks)
        args.build_directory = work/'build'
        build_info = build(args.build_directory)
        if existing and (work/'preparation.json').exists():
            preparation = read_json(work/'preparation.json')
            if preparation['source_sha256'] != hashes or preparation['structural'].get('status') != 'complete':
                raise ValueError('Invalid or incompatible shared preparation')
            preparation['build'] = build_info
        else:
            print('Recomputing shared structural prerequisites from zero...',flush=True)
            preparation = dict(source_sha256=hashes,build=build_info,structural=preflight.run())
            if sources() != hashes:
                raise ValueError('Sources changed during preparation')
        save_json(work/'preparation.json',preparation)
        args.preparation = preparation
        status,error = 0,None
        try:
            for rank in ranks:
                rank_args = argparse.Namespace(**vars(args))
                rank_args.rank = rank
                rank_args.work = work/f'rank-{rank:02d}'
                status = run_rank(rank_args)
                if status:
                    break
        except KeyboardInterrupt:
            status,error = 130,'Interrupted; incomplete ranks carry no mathematical conclusion.'
        except Exception:
            status,error = 1,traceback.format_exc()
        session.update(finished_utc=now(),seconds=round(time.monotonic()-tick,6),exit_status=status)
        batch_report(work,ranks,session,error)
        return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--rank',type=int,help='compute one exact rank')
    group.add_argument('--through',type=int,default=44,help='compute every rank from 1 to this rank (default: 44)')
    parser.add_argument('--jobs',type=int,default=8,help='maximum worker processes (default: 8)')
    parser.add_argument('--work',type=Path,default=Path('rank-work'),help='fresh work/output directory')
    parser.add_argument('--resume',action='store_true',help='resume this package\'s own unfinished batch')
    parser.add_argument('--report',action='store_true',help='verify saved accounting and rebuild the result ZIP')
    parser.add_argument('--verify',action='store_true',help='regenerate partitions and verify saved receipts and ledgers')
    parser.add_argument('--recheck',action='store_true',help='also repeat all terminal necessary-condition tests')
    parser.add_argument('--limit',type=int,help='pilot only: maximum new tasks per rank; omitted for a complete run')
    args = parser.parse_args()
    if args.jobs<1 or (args.limit is not None and args.limit<0):
        parser.error('jobs must be positive and limit nonnegative')
    if not 1 <= (args.rank if args.rank is not None else args.through) <= 44:
        parser.error('the certified search supports ranks 1 through 44')
    if sys.version_info<(3,9) or not sys.platform.startswith('linux'):
        parser.error('use Linux (including WSL2) and Python 3.9 or later')
    if sys.flags.optimize:
        parser.error('run without Python -O: structural proof assertions must remain enabled')
    try:
        return run(args)
    except BlockingIOError:
        print('This work directory is already being used by another run.',file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print('Interrupted during preparation. Resume the same rank selection; no unfinished task is certified.',file=sys.stderr)
        return 130
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())


