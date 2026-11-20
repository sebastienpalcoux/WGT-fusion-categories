#!/usr/bin/env python3
"""Reproduce the computations for Theorem A.

python3 verify.py                 Replay terminal tests and both odd censuses.
python3 verify.py --partitions    Also regenerate every native task partition.
python3 verify.py --full --jobs 8 Recompute the general census from its sources.

The default needs Python 3.10 or newer. The two native options additionally
need Linux and a C++17 compiler. The full run compares the sorted raw records
at every rank, including all native and arbitrary-precision completions.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'census/reference-census.zip'
PYTHON = [sys.executable, '-E', '-B']


def run(command, cwd=ROOT):
    print(' '.join(str(part) for part in command), flush=True)
    subprocess.run([*PYTHON, *command], cwd=cwd, check=True)


@contextmanager
def producer_sources():
    with TemporaryDirectory(prefix='rank-verification-') as temporary:
        directory = Path(temporary)
        with ZipFile(ARCHIVE) as archive:
            manifest = json.loads(archive.read('preparation.json'))['source_sha256']
            for name, expected in manifest.items():
                relative = PurePosixPath(name)
                if relative.is_absolute() or '..' in relative.parts:
                    raise ValueError('Invalid producer source path')
                data = archive.read('source/' + name)
                if sha256(data).hexdigest() != expected:
                    raise ValueError(f'Source checksum mismatch at {name}')
                target = directory / 'source' / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
        yield directory


def raw_counts(archive, rank):
    prefix = f'ranks/{rank:02d}/'
    data = gzip.decompress(archive.read(prefix + 'raw-census.jsonl.gz'))
    return Counter(data.splitlines())


def compare_censuses(fresh):
    sys.path.insert(0, str(ROOT / 'census'))
    import audit
    reference_summary = audit.audit(ARCHIVE, recheck_prime=False)
    fresh_summary = audit.audit(fresh, recheck_prime=False)
    if not (reference_summary['empty_complete_census'] and fresh_summary['empty_complete_census']):
        raise ValueError('Incomplete census')
    if reference_summary['totals'] != fresh_summary['totals']:
        raise ValueError('Task or exact continuation accounting differs')
    with ZipFile(ARCHIVE) as reference, ZipFile(fresh) as rebuilt:
        for rank in range(1, 45):
            if raw_counts(reference, rank) != raw_counts(rebuilt, rank):
                raise ValueError(f'Raw general census differs at rank {rank}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--partitions', action='store_true', help='compile the native generators and compare their partitions')
    parser.add_argument('--full', action='store_true', help='regenerate the general census through rank 44 from zero')
    parser.add_argument('--jobs', type=int, default=min(3, os.cpu_count() or 1), help='workers for the odd-square replay and full general generation')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    run(['census/recheck.py'])
    run(['odd_nonsquare/check.py'])
    run(['odd_square/check.py', '--jobs', str(args.jobs)])
    for script in ('three_prime_cutoff', 'support_five', 'isolation44'):
        run([f'arithmetic/{script}.py'])
    with producer_sources() as directory:
        source = directory / 'source'
        run(['arithmetic/packet_checks.py', str(source / 'arithmetic')])
        if args.full:
            run(['run.py', '--through', '44', '--jobs', str(args.jobs),
                 '--work', str(directory / 'fresh')], cwd=source)
            compare_censuses(directory / 'fresh/rank-results.zip')
        else:
            run(['-c', 'import preflight; preflight.run()'], cwd=source)
    if args.partitions and not args.full:
        run(['regenerate_tasks.py', str(ARCHIVE)])
    print('All selected checks passed.', flush=True)


if __name__ == '__main__':
    main()
