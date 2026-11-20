#!/usr/bin/env python3
"""Verify the general and odd rank computations using the saved complete runs.

The default checks every covering-census Table 1 row and the bounded
arithmetic tests, and replays both odd-dimensional certificates. --quick
skips the two odd certificate recursions. --independent adds the standalone
finite verifier (SymPy required). Denominator generation is separate.
"""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ARCHIVE = 'census/reference-census.zip'
# Each tuple is a script path followed by its arguments, relative to the
# supplement root. Add any new receipt checker after its prerequisites.
QUICK_CHECKS = (
    ('census/five_tests.py',),
    ('arithmetic/three_prime_cutoff.py',),
    ('arithmetic/support_five_moments.py',),
    ('arithmetic/odd_edge_mass.py',),
    ('arithmetic/gram_packets.py',),
    ('census/recheck.py',),
    ('verification/rank_tables/audit_slices.py', ARCHIVE),
    ('verification/rank_tables/audit_overlap.py', ARCHIVE),
    ('verification/rank_tables/audit_fusion.py', ARCHIVE),
    ('verification/odd_theory/check_edge37.py',),
    ('arithmetic/three_prime_rank.py',),
    ('verification/check_global_weyl_independent.py',),
    ('verification/sector_frontier/check.py',),
    ('verification/sector_frontier/verify_neutral.py',),
    ('verification/sector_frontier/verify_component.py',),
    ('verification/frontier217/check.py', '--graphs'),
    ('arithmetic/perfect_counts.py',),
    ('verification/odd_receipts/check.py',),
)
COMPARISON_ODD = (
    ('verify_odd_receipts.py',),
    ('verification/odd_uniform/check.py',),
    ('verification/odd_uniform/test_conditions.py',),
    ('verification/sparse_charge/check.py',),
    ('verification/nonsquare/check.py',),
    ('verification/isolated/inputs.py',),
    ('verification/profile_census/check.py', '--regenerate-graphs'),
    ('verification/isolated/check.py',),
    ('verification/mixed/test_fiber_sign.py',),
    ('verification/mixed/check_201_dense.py',),
    ('verification/mixed/check.py',),
)
FULL_CHECKS = (
    'arithmetic/support_five.py', 'arithmetic/large_four.py',
    'arithmetic/isolation44.py', 'arithmetic/two_five44.py',
    'arithmetic/high_four44.py', 'arithmetic/frontier_cases.py',
    'arithmetic/gauss_divisors.py', 'arithmetic/odd_nonsquare.py',
    'arithmetic/odd_square.py',
)


def run_checks(root, args):
    checks = list(QUICK_CHECKS)
    if not args.quick:
        checks.extend([('odd_nonsquare/check.py',), ('odd_square/check.py', '--jobs', str(args.jobs))])
    if args.independent:
        checks.append(('verification/independent/check.py',))
    if args.independent_odd:
        checks.extend(COMPARISON_ODD)
        checks.append(('verification/mixed/check.py', '--independent-replay'))
    if args.native_prefix:
        checks.extend([('verification/prefix/check.py',), ('verification/native_boundary/check.py',)])
    if args.full:
        if args.native_odd:
            checks.append(('build_odd_native.py',))
        for script in FULL_CHECKS:
            checks.append((script, '--jobs', str(args.jobs)) if script.endswith('/odd_square.py') else (script,))
    for command in checks:
        print('Checking', *command, flush=True)
        env = os.environ.copy()
        env['ODD_NATIVE'] = '1' if args.native_odd and command[0].endswith('/odd_square.py') else '0'
        subprocess.run([sys.executable, '-E', '-B', *command], cwd=root, env=env, check=True)
    for directory in (root, root / 'census'):
        print('Checking unit tests in', directory.relative_to(root) or '.', flush=True)
        subprocess.run([sys.executable, '-E', '-B', '-m', 'unittest', 'discover', '-s', 'tests'],
                       cwd=directory, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the two odd certificate recursions; retain bounded finite checks and accounting tests')
    parser.add_argument('--independent', action='store_true', help='also verify the independent finite receipt and all terminal witnesses; requires SymPy')
    parser.add_argument('--full', action='store_true', help='also replay the original final odd searches and full ancillary arithmetic')
    parser.add_argument('--jobs', type=int, default=1, help='workers for the odd-square replays')
    parser.add_argument('--native-odd', action='store_true', help='use the C++17/GMP continuation with --full')
    parser.add_argument('--independent-odd', action='store_true', help='also repeat the rank-209 mixed proof without the neutral-fiber theorem')
    parser.add_argument('--native-prefix', action='store_true', help='also compare a small native-prefix sample and force width handoffs, requiring g++ with UBSan')
    args = parser.parse_args()
    if args.jobs < 1 or (args.native_odd and not args.full):
        parser.error('--jobs must be positive and --native-odd requires --full')
    if args.quick and (args.full or args.independent_odd):
        parser.error('--quick cannot be combined with --full or --independent-odd')
    start = time.monotonic()
    # Some inherited arithmetic scripts write timing receipts next to their
    # sources. Replaying in a temporary copy preserves the release files.
    with tempfile.TemporaryDirectory(prefix='rank-bound-replay-') as temporary:
        root = Path(temporary) / 'rank-bound-supplement'
        shutil.copytree(Path(__file__).resolve().parent, root,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        run_checks(root, args)
    print(f'All selected checks passed ({time.monotonic()-start:.1f} s).', flush=True)
    print('Denominator enumeration was not repeated.')
    if args.quick:
        print('The odd-square and odd-nonsquare production certificate recursions were not repeated.')
    if not args.full:
        print('The auxiliary rank-201 Gauss replay was not repeated.')


if __name__ == '__main__':
    main()
