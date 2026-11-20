#!/usr/bin/env python3
"""Independently audit the portable fresh ranks1--44 result archive.

Standard library only. No driver or archived result is imported or executed.
This verifies coverage of the source-bound DECLARED task partitions; it does
not regenerate native partitions or replay the searches. Optional independent
prime-square trace checks repeat that mathematical necessary condition only.

Exit0: complete, internally consistent archive (survivors may be nonempty).
Exit2: internally consistent partial archive; no empty-census conclusion.
Exit1: damaged, inconsistent, unsupported, or otherwise invalid archive.
"""
import argparse
import ast
from collections import Counter, defaultdict
from functools import lru_cache
import gzip
import hashlib
from io import BytesIO
import json
from math import gcd, isfinite, isqrt
from pathlib import Path, PurePosixPath
import re
import sys
import time
from zipfile import ZipFile


FORMAT = 'uniform-rank-census-v2'
RANKS = list(range(1, 45))
SOURCE_DIRS = ('engine', 'engine_src', 'arithmetic', 'phase-mass', 'tests')
SOURCE_SUFFIXES = {'.py', '.cpp', '.hpp', '.json'}
HASH = re.compile(r'[0-9a-f]{64}')
PRIMES = tuple(p for p in range(2, 45)
               if all(p % d for d in range(2, isqrt(p)+1)))
COUNTERS = ('nodes', 'native', 'handoffs', 'exact_rows', 'profiles')
# Optional documentation and validation attachments, never substitutes for
# manifest-bound sources or census records. No general prefix is permitted.
SUPPLEMENTAL_NAMES = frozenset({
    'README.txt',
    'verification/payload-audit.json',
    'verification/verify_results.py',
    'verification/compare_rank44_censuses.py',
    'verification/census-comparison.json',
    'verification/auditor-tests.json',
    'verification/controlled-benchmark.json',
    'verification/native-pruning-validation.json',
})

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


class AuditError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise AuditError(message)


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode('utf-8')


def digest(value):
    return hashlib.sha256(encode(value)).hexdigest()


def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha_file(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            value.update(block)
    return value.hexdigest()


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def finite_seconds(value):
    return type(value) in (int, float) and isfinite(value) and value >= 0


def _object(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, 'Duplicate JSON object key: '+key)
        out[key] = value
    return out


def _bad_constant(value):
    raise AuditError('Nonfinite JSON number: '+value)


def loads(data):
    return json.loads(data, object_pairs_hook=_object, parse_constant=_bad_constant)


def json_member(archive, name):
    return loads(archive.read(name))


def json_lines(archive, name):
    with archive.open(name) as stream:
        for number, line in enumerate(stream, 1):
            require(line.endswith(b'\n') and line.strip(),
                    f'Truncated or empty JSONL record at {name}:{number}')
            yield loads(line)


@lru_cache(maxsize=4096)
def factors_small(n):
    """Factor a number whose support is contained among primes below45."""
    require(integer(n, 1), 'Expected a positive integer to factor')
    result = []
    for p in PRIMES:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            result.append((p, exponent))
    require(n == 1, 'Dimension root has a prime outside the supported rank range')
    return tuple(result)


def validate_task(task, rank):
    require(isinstance(task, dict), f'Rank{rank}: task is not an object')
    q, b = task.get('q'), task.get('B')
    require(integer(q, 1) and integer(b, 2), f'Rank{rank}: invalid q or B')
    f = q*b*b
    require(8 <= f <= rank-2 and f not in (8, 9),
            f'Rank{rank}: task has an excluded starting codegree')
    require(all(e == 1 for _, e in factors_small(q)), 'Task kernel is not squarefree')
    if q == 1:
        require(b >= 4 and 'm' not in task and integer(task.get('cost'), 1),
                'Invalid square task schema')
    else:
        require(integer(task.get('m'), 1) and task['m'] < rank,
                'Invalid nonsquare multiplicity')
    prefix = task.get('prefix', [b])
    require(isinstance(prefix, list) and 1 <= len(prefix) <= rank-3
            and all(integer(z, b) for z in prefix) and prefix[0] == b
            and prefix == sorted(prefix), 'Invalid denominator task prefix')
    return q, b, task.get('m'), prefix


def validate_partition(tasks, rank):
    require(isinstance(tasks, list), 'Task partition is not a list')
    seen, groups = set(), defaultdict(list)
    for task in tasks:
        token = digest(task)
        require(token not in seen, f'Rank{rank}: duplicate task')
        seen.add(token)
        q, b, m, prefix = validate_task(task, rank)
        groups[q, b, m].append(prefix)
    for prefixes in groups.values():
        ordered = sorted(prefixes)
        require(not any(b[:len(a)] == a for a, b in zip(ordered, ordered[1:])),
                f'Rank{rank}: overlapping task prefixes')


def validate_degree_row(row, task, rank):
    require(isinstance(row, list) and len(row) == rank and rank >= 2
            and all(integer(d, 1) for d in row) and row[0] == 1
            and row[1] >= 6 and row == sorted(row) and gcd(*row[1:]) == 1,
            f'Rank{rank}: invalid degree-list shape or nonunit gcd')
    q, b = task['q'], task['B']
    M = b*row[-1]
    require(sum(d*d for d in row) == q*M*M and all(M % d == 0 for d in row),
            f'Rank{rank}: dimension or denominator identity failed')
    factors = factors_small(M)
    require(all(p < rank for p, _ in factors), 'Degree row escapes native prime support')
    prefix = task.get('prefix', [b])
    require(all(d*z == M for d, z in zip(reversed(row), prefix)),
            'Degree row escapes its task prefix')
    if q != 1:
        require(row.count(row[-1]) == task['m'], 'Degree row escapes its multiplicity task')
    return M, factors


def obstruction_names(source_text):
    """Read only string literals in OBSTRUCTIONS assignments; execute no code."""
    names = set()
    for text in source_text:
        tree = ast.parse(text)
        for statement in tree.body:
            if isinstance(statement, ast.Assign) and any(
                    isinstance(target, ast.Name) and target.id == 'OBSTRUCTIONS'
                    for target in statement.targets):
                names.update(node.value for node in ast.walk(statement.value)
                             if isinstance(node, ast.Constant) and isinstance(node.value, str))
    require('prime_square_trace' in names, 'Source does not declare the new trace obstruction')
    return names


def check_sources(archive, expected, source_dir=None):
    require(isinstance(expected, dict) and expected, 'Missing source manifest')
    for name, claimed in expected.items():
        path = PurePosixPath(name)
        require(not path.is_absolute() and '..' not in path.parts
                and isinstance(claimed, str) and HASH.fullmatch(claimed),
                'Invalid source entry')
        require(sha_bytes(archive.read('source/'+name)) == claimed,
                'Embedded source checksum mismatch: '+name)
    embedded = {name[len('source/'):] for name in archive.namelist()
                if name.startswith('source/') and PurePosixPath(name).suffix in SOURCE_SUFFIXES}
    require(embedded == set(expected), 'Unlisted or missing code/data source file')
    if source_dir is not None:
        root = Path(source_dir).resolve()
        files = list(root.glob('*.py'))
        for directory in SOURCE_DIRS:
            files.extend(path for path in (root/directory).rglob('*')
                         if path.is_file() and path.suffix in SOURCE_SUFFIXES)
        external = {str(path.relative_to(root)): sha_file(path) for path in files}
        require(external == expected, 'Archive source manifest differs from the supplied source directory')
    return obstruction_names(archive.read('source/'+name).decode('utf-8')
                             for name in ('screen.py', 'terminal.py', 'engine_src/nonsquare.py'))


def check_preparation(preparation, sources):
    require(preparation.get('source_sha256') == sources, 'Preparation source binding differs')
    structural = preparation['structural']
    require(structural.get('status') == 'complete' and structural.get('fresh') is True
            and structural.get('ranks') == RANKS and structural.get('rank') == 44,
            'Shared structural preparation does not cover ranks1--44')
    require(structural.get('all_primary_survivors') == 0, 'Primary structural survivors remain')
    entries = structural['primary_profiles']
    require(len(entries) == 46
            and [entry['index'] for entry in entries] == list(range(46)),
            'Duplicate or missing primary-profile index')
    require(len(structural['mixed_signatures']) == 39
            and len(structural['mixed_two_five_excluded']) == 2,
            'Incomplete structural signature census')
    roots = Counter((entry['index'], root) for entry in entries for root in entry['roots'])
    closures = structural['root_closures']
    require(all(entry.get('status') == 'excluded' for entry in closures)
            and roots == Counter((entry['profile'], entry['root']) for entry in closures),
            'Primary structural root lacks an exact closure record')
    require(structural['mixed_two_five_check'].get('status') == 'complete'
            and structural['mixed_two_five_check'].get('fresh') is True,
            'Mixed two-five structural check is unfinished')
    lower = structural['lower_rank_structural_checks']
    require(lower.get('status') == 'complete' and lower.get('fresh') is True
            and lower.get('ranks') == RANKS
            and [entry['rank'] for entry in lower['mixed_checks']] == [40, 41, 42, 43]
            and all(entry.get('native_agreement') is True for entry in lower['mixed_checks']),
            'Lower-rank structural/native comparison is incomplete')
    return digest(structural)


def independent_prime_cut_failure(row, M, factors):
    """Direct degree-gcd/distinct-value cut check, importing no criterion."""
    primes = tuple(p for p, _ in factors)
    for mask in range(1, 1 << len(primes)):
        if mask == (1 << len(primes))-1:
            continue
        left = tuple(p for j, p in enumerate(primes) if mask >> j & 1)
        right = tuple(p for j, p in enumerate(primes) if not (mask >> j & 1))
        crossing = {d for d in row if d > 1
                    and any((M//d) % p == 0 for p in left)
                    and any((M//d) % p == 0 for p in right)}
        if len(crossing) < 3 or gcd(*crossing) != 1:
            return True
    return False


def independent_prime_square_failure(row, q, M, factors):
    """Repeat the mathematical trace test without importing the saved code."""
    if q != 1 and q not in PRIMES:
        return False
    counts = Counter(row)
    support_primes = [p for p, _ in factors]
    for p in support_primes:
        if q == p:
            continue
        n = M//p
        if n <= 1 or not counts[n]:
            continue
        if q == 1 and p >= 5:
            allowed = [j for j in range(1, p-1) if (p-1) % j == 0]
            if p == 5:
                allowed = [1]
            elif p == 7:
                allowed = [1, 2]
            if M % 2:
                odd_part = p-1
                while odd_part % 2 == 0:
                    odd_part //= 2
                allowed = [j for j in allowed if ((p-1)//j) % 2 == 1 and odd_part >= 5]
            if counts[n] not in {p*j for j in allowed}:
                return True
        target = n
        while target % p == 0:
            target //= p
        hits = [d for d in counts if d > 1 and target % d == 0]
        hit_support = {s for d in hits for s in support_primes if d % s == 0}
        if (p == 2 or n % p != 0) and gcd(*hits) == 1 and len(hit_support) >= 3:
            continue
        if q == 1:
            return True
        options = defaultdict(set)
        if q > 2:
            # Enumerate all tame inactive/active orbit divisor solutions.
            divisors = [a for a in range(1, q) if (q-1) % a == 0]
            for a in divisors:
                for b in divisors:
                    if b % ((q-1)//2):
                        continue
                    numerator = (q-1)*a
                    if numerator % b:
                        continue
                    k = isqrt(numerator//b)
                    if k < 2 or k*k*b != numerator:
                        continue
                    if counts[k*n]:
                        options[k*n].add((b, a))
                    if n % k == 0 and counts[n//k]:
                        options[n//k].add((a, b))
        capacity = 0
        for partner, choices in options.items():
            choices = tuple(sorted(choices))
            @lru_cache(None)
            def best(index, left):
                if index == len(choices):
                    return 0
                cost, gain = choices[index]
                return max(j*gain+best(index+1, left-j*cost)
                           for j in range(left//cost+1))
            capacity += best(0, counts[partner])
        if counts[n] > capacity:
            return True
    return False


def audit_rank(archive, rank, sources, preparation_sha, reasons, global_seen, recheck_prime):
    base = f'ranks/{rank:02d}/'
    spec = json_member(archive, base+'manifest.json')
    require(spec.get('format') == FORMAT and type(spec.get('rank')) is int
            and spec['rank'] == rank, 'Rank manifest has an incompatible format or rank')
    require(digest({k: v for k, v in spec.items() if k != 'sha256'}) == spec.get('sha256'),
            f'Rank{rank}: manifest checksum mismatch')
    require(spec['source_sha256'] == sources and spec['preflight_sha256'] == preparation_sha,
            f'Rank{rank}: source or preparation binding differs')
    tasks = spec['tasks']
    require(digest(tasks) == spec['partition_sha256'], f'Rank{rank}: partition checksum mismatch')
    validate_partition(tasks, rank)
    done = {}
    for rec in json_lines(archive, base+'receipts.jsonl'):
        index = rec.get('index')
        require(integer(index) and index < len(tasks) and index not in done,
                f'Rank{rank}: duplicate or out-of-range receipt index')
        require(rec.get('manifest_sha256') == spec['sha256'], 'Receipt manifest binding differs')
        require(all(integer(rec.get(key)) for key in COUNTERS)
                and finite_seconds(rec.get('seconds')), 'Invalid receipt counter or timing')
        require(type(rec.get('native_exit')) is int and rec['native_exit'] == 0
                and rec.get('exact_continuations_exhausted') is True,
                'Unfinished native task or exact continuation was counted complete')
        require(rec['profiles'] == rec['native']+rec['exact_rows']
                and (not rec['exact_rows'] or rec['handoffs']), 'Receipt backend count mismatch')
        excluded, retained = rec.get('excluded'), rec.get('survivors')
        require(isinstance(excluded, dict) and set(excluded) <= reasons
                and all(integer(n) for n in excluded.values()) and isinstance(retained, list)
                and sum(excluded.values())+len(retained) == rec['profiles'],
                'Receipt disposition accounting mismatch')
        require(rec.get('raw_ledger') == f'rows/{index:06d}.jsonl.gz'
                and isinstance(rec.get('raw_records_sha256'), str)
                and HASH.fullmatch(rec['raw_records_sha256'])
                and isinstance(rec.get('raw_ledger_sha256'), str)
                and HASH.fullmatch(rec['raw_ledger_sha256']), 'Invalid raw-ledger binding')
        task = tasks[index]
        completion = (f'DONE nodes={rec["nodes"]} profiles={rec["native"]} exact={rec["handoffs"]}'
                      if task['q'] == 1 else
                      f'rank {rank}, q {task["q"]}, B {task["B"]}, multiplicity {task["m"]}, '
                      f'nodes {rec["nodes"]}, candidates {rec["native"]}, exact {rec["handoffs"]}')
        require(rec.get('native_completion') == completion, 'Native completion trailer disagrees')
        done[index] = rec
    require(list(done) == sorted(done), f'Rank{rank}: receipts are not in task-index order')
    hashes = {index: hashlib.sha256() for index in done}
    row_counts = {index: Counter() for index in done}
    rejected = {index: Counter() for index in done}
    retained = {index: [] for index in done}
    compressed = archive.read(base+'raw-census.jsonl.gz')
    last_index = -1
    prime_checks = 0
    with gzip.GzipFile(fileobj=BytesIO(compressed), mode='rb') as stream:
        for number, line in enumerate(stream, 1):
            require(line.endswith(b'\n') and line.strip(), 'Invalid merged raw JSONL record')
            record = loads(line)
            require(isinstance(record, list) and len(record) == 5, 'Invalid merged raw record schema')
            r, index, backend, reason, row = record
            require(type(r) is int and r == rank and integer(index) and index in done
                    and index >= last_index, f'Rank{rank}: raw row lacks an ordered complete task')
            last_index = index
            require(backend in ('native', 'exact') and (reason is None or reason in reasons),
                    'Unknown raw backend or terminal obstruction')
            M, factors = validate_degree_row(row, tasks[index], rank)
            token = tuple(row)
            require(token not in global_seen, f'Rank{rank}: duplicate degree row across raw ledgers')
            global_seen.add(token)
            hashes[index].update(encode([backend, reason, row])+b'\n')
            row_counts[index][backend] += 1
            if reason is None:
                retained[index].append(row)
            else:
                rejected[index][reason] += 1
            if recheck_prime and reason == 'prime_square_trace':
                require(independent_prime_square_failure(row, tasks[index]['q'], M, factors),
                        f'Rank{rank}, task{index}: independent prime-square trace check did not exclude')
                prime_checks += 1
            if recheck_prime and reason == 'prime_cut':
                require(independent_prime_cut_failure(row, M, factors),
                        f'Rank{rank}, task{index}: independent prime-cut check did not exclude')
    for index, rec in done.items():
        require(hashes[index].hexdigest() == rec['raw_records_sha256'],
                f'Rank{rank}, task{index}: semantic raw-row digest mismatch')
        require(row_counts[index]['native'] == rec['native']
                and row_counts[index]['exact'] == rec['exact_rows']
                and dict(rejected[index]) == rec['excluded'] and retained[index] == rec['survivors'],
                f'Rank{rank}, task{index}: raw records disagree with receipt')
    expected_survivors = [dict(task=index, degrees=row) for index in sorted(done)
                          for row in retained[index]]
    require(list(json_lines(archive, base+'survivors.jsonl')) == expected_survivors,
            f'Rank{rank}: survivor file disagrees with raw records')
    summary = json_member(archive, base+'summary.json')
    sessions = json_member(archive, base+'sessions.json')
    require(isinstance(sessions, list) and sessions, 'Missing per-rank session history')
    require(sessions[0].get('fresh') is True and sessions[0].get('completed_before') == 0,
            f'Rank{rank}: first recorded session is not a fresh start')
    for session in sessions:
        require(finite_seconds(session.get('seconds'))
                and integer(session.get('completed_before')) and integer(session.get('completed_after'))
                and session['completed_before'] <= session['completed_after'] <= len(tasks),
                'Invalid session timing or task accounting')
    missing = [index for index in range(len(tasks)) if index not in done]
    errors = summary.get('error')
    status = 'complete' if not missing and errors is None else 'partial'
    totals = {key: sum(rec[key] for rec in done.values()) for key in COUNTERS}
    excluded = Counter()
    for rec in done.values():
        excluded.update(rec['excluded'])
    cases = {}
    for b, q in sorted({(task['B'], task['q']) for task in tasks}):
        indices = [i for i in done if (tasks[i]['B'], tasks[i]['q']) == (b, q)]
        cases[f'{b},{q}'] = dict(B=b, q=q, completed=len(indices),
                                tasks=sum((t['B'], t['q']) == (b, q) for t in tasks),
                                **{key: sum(done[i][key] for i in indices) for key in COUNTERS})
    expected = dict(format=FORMAT, rank=rank, status=status, tasks=len(tasks), completed=len(done),
                    missing_tasks=missing,
                    **{key: value for key, value in totals.items() if key != 'nodes'},
                    survivor_count=len(expected_survivors), survivor_records=len(expected_survivors),
                    excluded_by=dict(excluded), cases=cases, manifest_sha256=spec['sha256'],
                    source_sha256=digest(sources), raw_census_sha256=sha_bytes(compressed),
                    row_accounting_verified=True)
    require(all(summary.get(key) == value for key, value in expected.items()),
            f'Rank{rank}: summary disagrees with independently counted records')
    interrupted = summary.get('interrupted_tasks')
    require(isinstance(interrupted, list) and interrupted == sorted(set(interrupted))
            and all(integer(i) and i in missing for i in interrupted),
            'Interrupted tasks were also certified complete')
    worker_seconds = round(sum(rec['seconds'] for rec in done.values()), 6)
    compute_seconds = round(sum(s['seconds'] for s in sessions if s['operation'] == 'compute'), 6)
    require(abs(summary['worker_task_seconds']-worker_seconds) <= 0.000001
            and abs(summary['compute_session_seconds']-compute_seconds) <= 0.000001,
            'Summary timing totals disagree with records')
    projected = dict(rank=rank, status=status, tasks=len(tasks), completed=len(done),
                     profiles=totals['profiles'],
                     square=sum(c['profiles'] for c in cases.values() if c['q'] == 1),
                     nonsquare=sum(c['profiles'] for c in cases.values() if c['q'] != 1),
                     survivors=len(expected_survivors), compute_seconds=summary['compute_session_seconds'],
                     worker_seconds=summary['worker_task_seconds'], raw_census_sha256=sha_bytes(compressed),
                     manifest_sha256=spec['sha256'])
    report = dict(projected, native=totals['native'], exact_rows=totals['exact_rows'],
                  handoffs=totals['handoffs'], nodes=totals['nodes'], excluded_by=dict(excluded),
                  prime_square_rechecks=prime_checks, missing_tasks=missing,
                  partition_sha256=spec['partition_sha256'])
    return projected, report


def audit(path, source_dir=None, recheck_prime=False):
    started = time.monotonic()
    with ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), 'Duplicate ZIP member name')
        require(all(not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts
                    for name in names), 'Unsafe ZIP member name')
        require(archive.testzip() is None, 'ZIP member CRC failure')
        summary = json_member(archive, 'summary.json')
        batch = json_member(archive, 'batch.json')
        preparation = json_member(archive, 'preparation.json')
        require(summary.get('format') == FORMAT and batch.get('format') == FORMAT,
                'Unsupported archive format')
        requested = summary.get('requested_ranks')
        require(isinstance(requested, list) and all(type(r) is int for r in requested)
                and requested == RANKS and batch.get('ranks') == RANKS,
                'This auditor requires requested ranks1 through44, in order')
        sources = summary['source_sha256']
        require(batch.get('source_sha256') == sources, 'Batch source binding differs')
        reasons = check_sources(archive, sources, source_dir)
        preparation_sha = check_preparation(preparation, sources)
        ranks = sorted(int(match.group(1)) for name in names
                       if (match := re.fullmatch(r'ranks/(\d{2})/summary\.json', name)))
        require(len(ranks) == len(set(ranks)) and set(ranks) <= set(RANKS), 'Invalid rank coverage')
        expected_names = {'summary.json', 'summary.txt', 'batch.json', 'preparation.json'}
        expected_names.update('source/'+name for name in sources)
        expected_names.update(name for name in ('source/README.txt', 'source/Prime-square-criterion.tex')
                              if name in names)
        expected_names.update(SUPPLEMENTAL_NAMES.intersection(names))
        projections, reports, seen = [], [], set()
        for rank in ranks:
            base = f'ranks/{rank:02d}/'
            expected_names.update(base+name for name in (
                'summary.json', 'summary.txt', 'manifest.json', 'sessions.json',
                'raw-census.jsonl.gz', 'receipts.jsonl', 'survivors.jsonl'))
            projected, report = audit_rank(archive, rank, sources, preparation_sha,
                                           reasons, seen, recheck_prime)
            projections.append(projected)
            reports.append(report)
        require(set(names) == expected_names, 'Unexpected, duplicate, or missing archive member')
        complete = ranks == RANKS and all(r['status'] == 'complete' for r in reports)
        complete = complete and summary.get('error') is None
        status = 'complete' if complete else 'partial'
        require(summary.get('status') == status and summary.get('ranks') == projections,
                'Aggregate rank coverage or per-rank summary differs')
        totals = {key: sum(r[key] for r in reports)
                  for key in ('tasks', 'completed', 'profiles', 'native', 'exact_rows',
                              'handoffs', 'nodes', 'survivors', 'prime_square_rechecks')}
        require(all(summary.get(key) == totals[key] for key in ('tasks', 'completed', 'survivors')),
                'Aggregate task or survivor totals differ')
        excluded = Counter()
        for r in reports:
            excluded.update(r['excluded_by'])
        require(sum(excluded.values())+totals['survivors'] == totals['profiles'] == len(seen),
                'Global raw/disposition/uniqueness accounting differs')
        if 'rank_lower_bound_in_stated_scope' in summary:
            require(complete and totals['survivors'] == 0
                    and summary['rank_lower_bound_in_stated_scope'] == 45,
                    'A partial or nonempty census was reported as a rank exclusion')
    output = dict(audit_format='independent-uniform-rank-audit-v1', integrity='passed',
                  archive=str(Path(path).resolve()), archive_sha256=sha_file(path),
                  verifier_sha256=sha_file(__file__), requested_ranks=RANKS,
                  coverage=status, covered_ranks=ranks,
                  missing_ranks=[r for r in RANKS if r not in ranks],
                  source_files=len(sources), source_manifest_sha256=digest(sources),
                  source_binding='matches supplied directory' if source_dir else 'embedded manifest and content',
                  preflight_sha256=preparation_sha, totals=totals,
                  excluded_by=dict(sorted(excluded.items())), ranks=reports,
                  empty_complete_census=complete and totals['survivors'] == 0,
                  prime_square_rechecks_requested=bool(recheck_prime),
                  seconds=round(time.monotonic()-started, 6),
                  not_replayed=['native task-partition generation', 'degree enumeration',
                                'structural arithmetic calculations',
                                'other terminal obstruction proofs'],
                  compressed_task_ledgers='Not embedded; every task is bound by its semantic row digest and the merged census hash.',
                  scope='Nonpointed based-simple integral rotation-admissible modular data.')
    if output['empty_complete_census']:
        output['rank_lower_bound_supported_in_stated_scope'] = 45
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--source-dir', type=Path,
                        help='also compare every source file with this unchanged directory')
    parser.add_argument('--recheck-prime-square', action='store_true',
                        help='independently repeat each recorded prime_square_trace exclusion')
    parser.add_argument('--output', type=Path, help='save the JSON audit outside the frozen source tree')
    args = parser.parse_args()
    try:
        result = audit(args.archive, args.source_dir, args.recheck_prime_square)
        status = 0 if result['coverage'] == 'complete' else 2
    except Exception as error:
        result = dict(audit_format='independent-uniform-rank-audit-v1', integrity='failed',
                      archive=str(args.archive.resolve()), error=f'{type(error).__name__}: {error}')
        status = 1
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix+'.tmp')
        temporary.write_text(text)
        temporary.replace(args.output)
    print(text, end='')
    return status


if __name__ == '__main__':
    sys.exit(main())
