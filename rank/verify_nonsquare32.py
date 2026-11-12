#!/usr/bin/env python3
"""Check the exhaustive nonsquare rank-32 enumeration and its independent arithmetic screens."""
import argparse
from hashlib import sha256
from collections import defaultdict
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import tempfile
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'code'))
from current_screen import screen
from rank_audit import validate


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(blob):
    return sha256(blob).hexdigest()


def canonical(rows):
    return "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in sorted(rows)).encode("ascii")


def root_map(rank):
    """Map a possible largest-object codegree to its unique nonsquare root (q,b)."""
    roots = {}
    for q in range(2, rank):
        divisors = range(2, isqrt(q) + 1)
        if any(q % (i*i) == 0 for i in divisors):
            continue
        prime = all(q % i != 0 for i in divisors)
        for b in range(1, isqrt((rank-1)//q) + 1):
            f = q*b*b
            if f >= 8 and not (b == 1 and prime):
                roots[f] = (q, b)
    return roots


DATA = ROOT / 'data/revision/nonsquare'
SRC = ROOT / 'optimized/nonsquare'
RANK = 32


def verify(replay=False):
    manifest = json.loads((DATA / 'rank32.json').read_bytes())
    for name in ('nonsquare_joint.cpp', 'joint_capacity.hpp'):
        require(digest((SRC / name).read_bytes()) == manifest['sources'][name], 'wrong source ' + name)
    roots = root_map(RANK)
    raw = (DATA / 'rank32.jsonl').read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    require(canonical(rows) == raw and digest(raw) == manifest['sha256'], 'changed data')
    require(len(rows) == len(set(map(tuple, rows))) == manifest['types'], 'duplicate or missing list')
    jobs = manifest['roots']
    require([(r['q'], r['b']) for r in jobs] == sorted(roots.values()), 'incomplete root partition')
    grouped = defaultdict(list)
    for row in rows:
        D = validate(row, RANK)
        q, b = roots[D // row[-1]**2]
        M = isqrt(D // q)
        require(D == q*M*M and row[-1]*b == M and all(M % d == 0 for d in row), 'wrong root')
        grouped[q, b].append(row)
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'search'
        if replay:
            subprocess.run(['c++', '-O3', '-std=c++17', SRC / 'nonsquare_joint.cpp', '-o', exe], check=True)
        for item in jobs:
            q, b = item['q'], item['b']
            expected = canonical(grouped[q, b])
            require(digest(expected) == item['sha256'] and len(grouped[q,b]) == item['types'], 'wrong root output')
            if replay:
                job = subprocess.run([exe, str(RANK), '8', str(q), str(b)], capture_output=True, check=True)
                match = re.fullmatch(rb'Complete rank (\d+), nodes (\d+), types (\d+), seconds \S+', job.stderr.splitlines()[-1])
                require(match and list(map(int, match.groups())) == [RANK,item['nodes'],item['types']],
                        f'incomplete or different enumeration at {(q,b)}')
                require(canonical([json.loads(x) for x in job.stdout.splitlines()]) == expected,
                        f'different degree list at {(q,b)}')
    report = screen(rows, RANK)
    require(report == manifest['screen'] and report['survivors'] == [], 'unresolved or different exclusion')
    print(f"Rank 32, nonsquare: {len(roots)} roots, {len(rows):,} distinct lists, zero survivors.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true', help='recompile and enumerate every root')
    verify(parser.parse_args().replay)
