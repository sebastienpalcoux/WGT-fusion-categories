#!/usr/bin/env python3
"""Finish rank 28 from the 13 completed roots and the two outstanding GMP roots."""
import argparse
from collections import defaultdict
import fcntl
import hashlib
import io
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
sys.path.insert(0, str(ROOT / 'continuation'))
from continue_ranks import package, sha256
from rank_audit import fast_report
from verify_result import degree_vector, require, verify_one

MISSING = [[1, 3], [1, 4]]
FRONTIER = ROOT / 'data/frontier'
SEED = FRONTIER / 'rank28_partial.jsonl'
META = FRONTIER / 'rank28_partial.json'
GMP = ROOT / 'continuation/search.cpp'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encode(rows):
    return b''.join((json.dumps(t, separators=(',', ':'))+'\n').encode('ascii') for t in sorted(rows))


def canonical_rows(data):
    rows = [degree_vector(line, 28) for line in data.splitlines()]
    require(len(set(rows)) == len(rows) and encode(rows) == data, 'malformed or duplicate canonical type')
    return rows


def build_gmp(binary):
    subprocess.run(['c++','-O3','-std=c++17',str(GMP),'-lgmpxx','-lgmp','-o',str(binary)], check=True)


def check_seed(binary, replay=False):
    data = SEED.read_bytes()
    meta = json.loads(META.read_bytes())
    require(meta['rank'] == 28 and meta['complete'] is False and meta['missing_roots'] == MISSING,
            'incorrect frontier scope')
    require(meta['sha256'] == digest(data) and meta['gmp_sha256'] == sha256(GMP), 'seed/source digest mismatch')
    native = next((s for s in sorted((ROOT/'code').glob('*.cpp')) if sha256(s) == meta['canonical_native_sha256']), None)
    require(native is not None, 'archived native generator digest unavailable')
    found = subprocess.run([binary,'--manifest','28','1'], capture_output=True, check=True)
    roots = [(int(q), int(b)) for q,b in (line.split(b'\t') for line in found.stdout.splitlines())]
    marker = re.fullmatch(rb'COMPLETE manifest rank=28 nodes=\d+ types=0 shards=(\d+)\n', found.stderr)
    require(marker and int(marker[1]) == len(roots) == len(set(roots)), 'invalid GMP root manifest')
    require(digest(found.stdout) == meta['root_list_sha256'], 'depth-one manifest hash disagrees')
    listed = meta['roots']
    require(roots == [(v['q'],v['b']) for v in listed] and
            [[v['q'],v['b']] for v in listed if not v['complete']] == MISSING, 'frontier misses a root')
    complete = {key for key in roots if list(key) not in MISSING}
    rows = canonical_rows(data)
    require(len(rows) == meta['candidates'] == 414, 'expected 414 archived types')
    groups = defaultdict(list)
    for t in rows:
        f = sum(d*d for d in t)//t[-1]**2
        matching = [key for key in roots if key[0]*key[1]**2 == f]
        require(len(matching) == 1 and matching[0] in complete, 'seed contains an unknown or unfinished root')
        groups[matching[0]].append(t)
    for v in listed:
        if v['complete']:
            part = groups[v['q'],v['b']]
            require(len(part) == v['types'] and digest(encode(part)) == v['sha256'], 'seed root digest mismatch')
    if replay:
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp)/'native'
            subprocess.run(['c++','-O3','-std=c++17',str(native),'-o',str(exe)], check=True)
            for q,b in sorted(complete):
                r = subprocess.run([exe,'28','8',str(q),str(b)], capture_output=True, check=True)
                marks = re.findall(rb'^Complete rank 28, nodes \d+, types (\d+), seconds [\d.e+-]+$',r.stderr,re.M)
                again = [degree_vector(line,28) for line in r.stdout.splitlines()]
                require(len(marks) == 1 and int(marks[0]) == len(again) and sorted(again) == groups[q,b],
                        f'native root {q},{b} replay disagrees')
    return rows


def verify_archive(result, binary, seed, replay):
    with ZipFile(result) as archive:
        names = ['gmp_roots.zip','rank28_partial.json','rank28_partial.jsonl','candidates.jsonl','family_result.json']
        require(sorted(archive.namelist()) == sorted(names), 'unexpected or repeated return archive member')
        require(archive.read('rank28_partial.json') == META.read_bytes() and
                archive.read('rank28_partial.jsonl') == SEED.read_bytes(), 'return seed disagrees with published frontier')
        with ZipFile(io.BytesIO(archive.read('gmp_roots.zip'))) as gmp_result:
            def compiled(key):
                require(key == sha256(GMP), 'unknown GMP generator')
                return binary, False
            header = json.loads(gmp_result.read('search.json'))
            require(header['rank'] == 28, 'wrong rank in two-root certificate')
            verify_one(gmp_result, compiled, replay, roots=MISSING)
            extra = canonical_rows(gmp_result.read('candidates.jsonl'))
        merged = canonical_rows(archive.read('candidates.jsonl'))
        require(sorted(seed+extra) == merged, 'combined types do not match the disjoint roots')
        report = fast_report([list(t) for t in merged], 28)
        require(report == json.loads(archive.read('family_result.json')), 'joint arithmetic audit disagrees')
        print(f"Rank 28: {len(merged)} distinct profiles across all roots, {len(report['survivors'])} arithmetic survivors.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=12)
    parser.add_argument('--depth', type=int, default=11)
    parser.add_argument('--work', type=Path, default=Path('rank28-work'))
    parser.add_argument('--prepare-only', action='store_true', help='validate the published seed and generate the GMP shard manifest')
    parser.add_argument('--verify', type=Path, metavar='RETURN_ZIP', help='verify a completed independent return archive')
    parser.add_argument('--replay', action='store_true', help='with --verify, regenerate the 13 native and two GMP roots')
    args = parser.parse_args()
    if args.jobs < 1 or not 1 <= args.depth <= 25 or args.replay and not args.verify or args.prepare_only and args.verify:
        parser.error('invalid job count, split depth, or combination of modes')
    if args.verify:
        with tempfile.TemporaryDirectory() as temp:
            binary = Path(temp)/'gmp'; build_gmp(binary)
            verify_archive(args.verify, binary, check_seed(binary, args.replay), args.replay)
        return
    work = args.work.resolve(); work.mkdir(exist_ok=True,parents=True)
    with (work/'.finish.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        part = work/'gmp'
        cmd = [sys.executable,str(ROOT/'continuation/run_shards.py'),'28','--jobs',str(args.jobs),
               '--depth',str(args.depth),'--work',str(part),'--root','1','3','--root','1','4']
        subprocess.run([*cmd,'--prepare-only'],check=True)
        binary = part/'search'
        seed = check_seed(binary)
        print(f'{len(seed)} published profiles validated; two GMP roots selected.',flush=True)
        if args.prepare_only: return
        subprocess.run(cmd,check=True)
        with (part/'.lock').open('a') as part_lock:
            fcntl.flock(part_lock,fcntl.LOCK_EX)
            extra = canonical_rows((part/'candidates.jsonl').read_bytes())
            local = fast_report([list(t) for t in extra], 28)
            (part/'family_result.json').write_text(json.dumps(local,indent=2)+'\n')
            inner = package(28,part)
            merged = sorted(seed+extra)
            full = fast_report([list(t) for t in merged], 28)
            target = work/'rank28_completed.zip'; tmp = target.with_suffix('.zip.partial')
            with ZipFile(tmp,'w',compression=ZIP_DEFLATED) as archive:
                archive.write(inner, 'gmp_roots.zip')
                archive.write(SEED, SEED.name); archive.write(META, META.name)
                archive.writestr('candidates.jsonl',encode(merged))
                archive.writestr('family_result.json',json.dumps(full,indent=2)+'\n')
            verify_archive(tmp,binary,seed,False)
            tmp.replace(target)
            print(f'Return this file: {target}')

if __name__ == '__main__':
    main()
