#!/usr/bin/env python3
"""Short differential checks for partner-capacity and prime-cut pruning.

This never starts a rank census. Each selected native task and each exact
continuation has a five-second limit. The reference rows identify small
test inputs only; they are never used by the production enumerator.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import subprocess
import sys
import time
from tempfile import TemporaryDirectory
from zipfile import ZipFile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runner', type=Path, default=Path(__file__).resolve().parents[2]/'census')
    args = parser.parse_args()
    runner = args.runner.resolve()
    temporary = TemporaryDirectory(prefix='prefix-check-')
    here = Path(temporary.name)
    reference_results = here/'reference'
    with ZipFile(runner/'reference-census.zip') as archive:
        for suffix in ('py', 'cpp'):
            (here/f'nonsquare-before.{suffix}').write_bytes(
                archive.read(f'source/engine_src/nonsquare.{suffix}'))
        (here/'nonsquare-exact-before.py').write_bytes(
            archive.read('source/engine/nonsquare_exact.py'))
        for rank in (41, 43, 44):
            folder = reference_results/f'rank-{rank:02d}'
            folder.mkdir(parents=True)
            (folder/'manifest.json').write_bytes(archive.read(f'ranks/{rank:02d}/manifest.json'))
    for stage, source in (('before', here/'nonsquare-before.cpp'),
                          ('after', runner/'engine_src/nonsquare.cpp')):
        subprocess.run(['g++','-O2','-std=c++17',str(source),'-o',str(here/f'nonsquare-{stage}')],
                       check=True, capture_output=True, text=True, timeout=30)
    sys.path[:0] = [str(runner), str(runner/'engine_src'),
                   str(runner/'engine'), str(runner/'phase-mass')]
    new = load('capacity_after', runner/'engine_src/nonsquare.py')
    old = load('capacity_before', here/'nonsquare-before.py')
    old_helper = load('mass_before', here/'nonsquare-exact-before.py')
    old.remaining_mass_possible = old_helper.remaining_mass_possible
    from prime_square import prime_square_failure
    from prime_cut import prime_cut_failure
    from screen import factor

    # Open-class, passed-class, support, and p_b=2 boundary cases.
    cases = [
        ((5, 5), 4, (2, 3, 5), 3, True),
        ((5, 5), 3, (2, 3, 5), 3, False),
        ((5, 5, 10, 10), 2, (2, 3, 5), 3, True),
        ((5, 5, 10, 10, 12), 9, (2, 3, 5), 3, False),
        ((5, 5, 10, 10, 10, 10), 0, (2, 3, 5), 3, True),
        ((2, 2, 4), 3, (2, 3, 5), 3, True),
        ((2, 2, 4), 2, (2, 3, 5), 3, False),
        ((2, 2, 4, 4, 4, 4), 0, (2, 3, 5), 3, True),
        ((2, 2, 4, 4, 4, 5), 9, (2, 3, 5), 3, False),
        ((5, 5), 0, (2, 3, 5, 7), 3, True),
        ((5, 5), 0, (2, 3, 5), 5, True),
    ]
    for prefix, slots, ps, q, expected in cases:
        assert new._three_prime_block_capacity(
            slots, prefix, Counter(prefix), ps, q) is expected

    # For every finite partner completion in a small universe, rejection
    # of the prefix implies rejection of the completion. This checks
    # that an open or not-yet-reached class is never undercounted.
    from itertools import combinations_with_replacement
    projection_checks = 0
    universe = (2, 3, 4, 5, 6, 9, 10, 12)
    for length in range(1, 6):
        for prefix in combinations_with_replacement(universe, length):
            for slots in range(3):
                if new._three_prime_block_capacity(
                        slots, prefix, Counter(prefix), (2, 3, 5), 3):
                    continue
                future = tuple(z for z in universe if z >= prefix[-1])
                for suffix in combinations_with_replacement(future, slots):
                    full = prefix + suffix
                    assert not new._three_prime_block_capacity(
                        0, full, Counter(full), (2, 3, 5), 3)
                    projection_checks += 1

    native_results = []
    exact_inputs = []
    for rank, index in ((41, 395), (43, 90), (44, 0), (44, 860)):
        source = reference_results/f'rank-{rank:02d}'
        manifest = json.loads((source/'manifest.json').read_text())
        task = manifest['tasks'][index]
        tail = [str(rank), str(task['q']), str(task['B']), str(task['m'])]
        if 'prefix' in task:
            tail += ['--prefix', ','.join(map(str, task['prefix']))]
        observations = []
        rowsets = []
        for stage in ('before', 'after'):
            started = time.monotonic()
            result = subprocess.run([str(here/f'nonsquare-{stage}'), *tail],
                                    capture_output=True, text=True, timeout=5,
                                    check=True)
            rows = [json.loads(line) for line in result.stdout.splitlines()]
            assert all(isinstance(row, list) for row in rows), 'unexpected handoff'
            assert len(rows) == len({tuple(row) for row in rows})
            rowsets.append({tuple(row) for row in rows})
            observations.append(dict(stage=stage, rows=len(rows),
                elapsed_seconds=time.monotonic()-started,
                completion=result.stderr.strip(),
                output_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()))
        before, after = rowsets
        assert after <= before
        removed = before-after
        certificates = [prime_square_failure(row, task['q'], task['B']*row[-1])
                        for row in sorted(removed)]
        assert all(certificates)
        native_results.append(dict(rank=rank, task_index=index, task=task,
            manifest_sha256=sha(source/'manifest.json'),
            observations=observations, removed=len(removed),
            removed_all_certified_by_prime_free_theorem=True))
        # One retained and one newly removed positive input are sufficient
        # to exercise both outcomes of the exact terminal continuation.
        for label, candidates in (('retained', after), ('removed', removed)):
            if candidates and not any(x[0] == label for x in exact_inputs):
                exact_inputs.append((label, rank, task, min(candidates)))

    exact_results = []
    def alarm(*_):
        raise TimeoutError('five-second exact-continuation limit')
    signal.signal(signal.SIGALRM, alarm)
    for label, rank, task, original in exact_inputs:
        M = task['B']*original[-1]
        denominators = sorted(M//d for d in original[1:])
        item = {'exact': denominators[:-2], 'primes': [p for p, _ in factor(M)]}
        rowsets = []
        timings = []
        for module in (old, new):
            signal.alarm(5)
            started = time.process_time()
            result = {tuple(row) for row in module.exact_rows(item, task, rank)}
            timings.append(time.process_time()-started)
            signal.alarm(0)
            rowsets.append(result)
        before, after = rowsets
        assert original in before
        assert after <= before
        def current_certificate(row):
            M = task['B']*row[-1]
            return (prime_square_failure(row, task['q'], M)
                    or prime_cut_failure(row, M, tuple(p for p, _ in factor(M))))
        assert all(current_certificate(row) for row in before-after)
        if label == 'removed':
            assert original not in after
        elif original not in after:
            # A former partner-capacity survivor may now fail the cut.
            assert prime_cut_failure(original, M, tuple(p for p, _ in factor(M)))
        exact_results.append(dict(case=label, rank=rank, prefix=item,
            before_rows=len(before), after_rows=len(after), cpu_seconds=timings,
            removed_all_certified_by_current_terminal=True,
            original_row=list(original)))

    output = dict(scope='targeted differential validation, not a fresh census',
        reference_archive_sha256=sha(runner/'reference-census.zip'),
        source_sha256={
            'native_before': sha(here/'nonsquare-before.cpp'),
            'native_after': sha(runner/'engine_src/nonsquare.cpp'),
            'python_before': sha(here/'nonsquare-before.py'),
            'python_after': sha(runner/'engine_src/nonsquare.py'),
            'mass_before': sha(here/'nonsquare-exact-before.py'),
            'mass_after': sha(runner/'engine/nonsquare_exact.py'),
            'prime_cut': sha(runner/'prime_cut.py'),
            'prime_cut_prefix': sha(runner/'prime_cut_prefix.py'),
            'terminal_criterion': sha(runner/'prime_square.py'),
            'check_script': sha(__file__)},
        boundary_cases=len(cases), projection_checks=projection_checks,
        native=native_results, exact=exact_results, status='passed')
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
