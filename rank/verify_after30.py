#!/usr/bin/env python3
"""Verify the supplied rank-30 return, new projections, rank-32 data, and higher square slices."""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
DATA = Path('data/revision')
ARCHIVE = DATA / 'rank30_completed.zip'


def run(*args):
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


def main():
    run('verify_rank30.py')
    run('verify_nonsquare32.py')
    run('square_slices.py', '--stored')
    run('rank31_fourteen.py')
    run('code/check_dyadic_zero.py', ARCHIVE)
    run('code/check_dyadic_concentration.py', ARCHIVE, '--expect', DATA / 'dyadic_concentration.json')
    run('code/rank30_mass_union.py', '--verify', DATA / 'rank30_mass_union.json')
    septic = subprocess.check_output([sys.executable, 'code/check_septic.py', ARCHIVE], cwd=ROOT)
    if json.loads(septic) != json.loads((ROOT / DATA / 'septic_projection.json').read_bytes()):
        raise ValueError('septic projection disagrees')
    print('All completed-rank and projection checks passed.')


if __name__ == '__main__':
    main()
