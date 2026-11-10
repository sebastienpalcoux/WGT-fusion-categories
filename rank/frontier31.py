#!/usr/bin/env python3
"""Adaptive rank-31 continuation with exact prefix coverage and bounded-memory verification."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import groupby
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

import next_rank as search
from dyadic_upgrade import strengthened
import rank31_fourteen as fourteen

RANK = 31
ROOT = Path(__file__).resolve().parent
ENGINE = ROOT / 'optimized'
DEPTH = 20
GROUP = 12
GENERATION = "dyadic-4"


def build(tmp):
    exe = {}
    fingerprints = {}
    for name, libs in [('fast', []), ('exact', ['-lgmpxx', '-lgmp'])]:
        source = ENGINE / (name + '.cpp')
        fingerprints[name] = search.digest(source.read_bytes())
        exe[name] = tmp / name
        subprocess.run(['c++', '-O3', '-std=c++17', source, *libs, '-o', exe[name]], check=True)
    # Include the identical GMP engine when the driver starts inside a parent.
    wrapper = ROOT / 'current' / 'child_manifest.cpp'
    fingerprints['wrapper'] = search.digest(wrapper.read_bytes())
    local = tmp / 'child_manifest.cpp'
    local.write_bytes(wrapper.read_bytes())
    exe['split'] = tmp / 'split'
    subprocess.run(['c++', '-O3', '-std=c++17', '-I', ENGINE, local,
                    '-lgmpxx', '-lgmp', '-o', exe['split']], check=True)
    return exe, fingerprints


def build_second(tmp):
    exe = {}
    fingerprints = {}
    header = (ENGINE/'dyadic_mass.hpp').read_bytes()
    fingerprints['rules_header'] = search.digest(header)
    fingerprints['rule_integration'] = search.digest((ROOT/'dyadic_upgrade.py').read_bytes())
    for name, libs in [('fast', []), ('exact', ['-lgmpxx', '-lgmp'])]:
        revised = strengthened(name, (ENGINE/(name+'.cpp')).read_text(), header.decode('ascii'))
        source = tmp/('second_' + name + '.cpp')
        source.write_text(revised)
        fingerprints[name] = search.digest(source.read_bytes())
        exe[name] = tmp/('second_' + name)
        subprocess.run(['c++', '-O3', '-std=c++17', source, *libs, '-o', exe[name]], check=True)
    return exe, fingerprints


def read_seed(first_spec, keys):
    path = ROOT/'data'/'rank31'/'rank31_seed.zip'
    supplied = {}
    with ZipFile(path) as z:
        search.require(sorted(z.namelist()) == ['receipts.jsonl','spec.json'] and
                       json.loads(z.read('spec.json')) == first_spec, 'incompatible original seed')
        with z.open('receipts.jsonl') as records:
            for line in records:
                i, item = json.loads(line)
                search.require(type(i) is int and 0 <= i < len(keys) and i not in supplied, 'invalid seed key')
                search.check_receipt(item, RANK, keys[i])
                supplied[i] = item
    return supplied, search.digest(path.read_bytes())


def in_fourteen_root(key):
    q,b = key
    return q==1 and b.count(4)==14 and b[-1]>4


def closed_receipt(key):
    item = search.receipt(key,0,0,0,[])
    del item['sha256']
    item.update(generation=GENERATION, closure='rank-31 complete fourteen-fours root')
    return {**item,'sha256':search.digest(search.encoded(item))}


def make_receipt(i, key, first, second, supplied):
    if i in supplied:
        return search.compute(first, RANK, key)
    if in_fourteen_root(key):
        return closed_receipt(key)
    item = search.compute(second, RANK, key)
    del item['sha256']
    item['generation'] = GENERATION
    return {**item, 'sha256':search.digest(search.encoded(item))}


def check_generation(i, key, item, supplied):
    if i in supplied:
        search.require(item == supplied[i], 'original receipt differs from the frozen seed')
    elif in_fourteen_root(key):
        search.require(item==closed_receipt(key), 'incorrect root-closure receipt')
    else:
        search.require(item.get('generation')==GENERATION and 'closure' not in item, 'receipt from wrong kernel generation')


def grouping(denominators, prefix):
    # At most GROUP descendants at the probe depth are assigned to each job.
    if len(denominators) <= GROUP:
        yield prefix
        return
    level = len(prefix)
    for b, items in groupby(denominators, key=lambda t: t[level]):
        yield from grouping(list(items), (*prefix, b))


def plan(exe):
    roots = search.uniform_manifest(exe['exact'], RANK, 1)
    keys = []
    for q, root in roots:
        if q != 1 or root != (4,):
            keys.append((q, root))
            continue
        _, lines = search.exact_run(exe['split'], 'manifest', RANK,
                                    [str(DEPTH), '1', '4'])
        probes = []
        for line in lines:
            mult, text = line.split(b'\t')
            denoms = tuple(map(int, text.split(b',')))
            search.require(mult == b'1' and len(denoms) == DEPTH and denoms[0] == 4
                           and denoms == tuple(sorted(denoms)), 'invalid probe')
            probes.append(denoms)
        search.require(len(probes) == len(set(probes)), 'duplicate probe')
        # Every completion has a depth-20 prefix. Branches that terminate sooner
        # have no completed row because the last two terms lie at depth 28.
        if probes:
            keys.extend((q, p) for p in grouping(sorted(probes), root))
    keys.sort()
    for (q, a), (s, b) in zip(keys, keys[1:]):
        search.require(q != s or b[:len(a)] != a, 'overlapping terminal jobs')
    return keys


def summary(items):
    counts = Counter()
    survivors = []
    total = 0
    for key, item in items:
        rows = search.check_receipt(item, RANK, key)
        if not rows:
            continue
        part = search.screen(rows, RANK)
        counts.update(part['first_exclusions'])
        total += len(rows)
        survivors.extend(part['survivors'])
    return {'rank': RANK, 'profiles': total, 'first_exclusions': dict(sorted(counts.items())),
            'survivors': sorted(survivors)}


def spec(hashes, keys):
    return {'format': 'adaptive-31-v1', 'rank': RANK, 'depth': DEPTH, 'group': GROUP,
            'sources': hashes, 'manifest_sha256': search.digest(search.encoded(keys))}


def verify(path, first, second, keys, expected, supplied, replay=False):
    with ZipFile(path) as z:
        search.require(sorted(z.namelist()) == ['receipts.jsonl', 'screen.json', 'spec.json'],
                       'unexpected archive contents')
        search.require(json.loads(z.read('spec.json')) == expected, 'source or partition changed')
        with z.open('receipts.jsonl') as records:
            def stream():
                for i, key in enumerate(keys):
                    line = records.readline()
                    search.require(bool(line), 'missing receipt')
                    item = json.loads(line)
                    check_generation(i, key, item, supplied)
                    if replay:
                        new = make_receipt(i,key,first,second,supplied)
                        search.require(search.encoded(new) == search.encoded(item), 'search replay disagrees')
                    yield key, item
                search.require(not records.readline(), 'extra receipt')
            report = summary(stream())
        search.require(report == json.loads(z.read('screen.json')), 'screen replay disagrees')
    return report


def run(args, first, second, keys, first_spec, expected, supplied):
    work = args.work.resolve()
    with search.locked(work):
        metadata = work/'spec.json'
        def path(i):
            return work/f'shard-{i:06}.json'
        if metadata.exists():
            old = json.loads(metadata.read_bytes())
            if old == first_spec or (old.get('format') in ('adaptive-31-v2','adaptive-31-v3') and
                    all(old.get(k)==expected.get(k) for k in ('rank','depth','group','sources','manifest_sha256','seed_sha256'))
                    and old.get('remaining_sources',{}).get('rules_header') in (
                    '85d9c0094396d7edc4970afecc0732f9c03833f7496763b29a33c62e0a595111',
                    'e97e2cf58337ba453a31ec0d4426821a2ea388aacad9351e3e4b0adcc814916a')):
                # The original seed remains valid. Unseeded first-generation work
                # is retained separately because its issued proof rules differ.
                archive = work/('archived-'+old['format'])
                archive.mkdir(exist_ok=True)
                for i in range(len(keys)):
                    if i not in supplied and path(i).exists():
                        dest = archive/path(i).name
                        search.check_receipt(json.loads(path(i).read_bytes()),RANK,keys[i])
                        if dest.exists():
                            search.require(dest.read_bytes()==path(i).read_bytes(),'incompatible saved first-generation work')
                        path(i).replace(dest)
                search.save(metadata, search.encoded(expected))
            else:
                search.require(old == expected, 'incompatible checkpoint directory')
        else:
            search.save(metadata, search.encoded(expected))
        for i,item in supplied.items():
            if not path(i).exists():
                search.save(path(i), search.encoded(item))
        closed = 0
        for i,key in enumerate(keys):
            if i not in supplied and in_fourteen_root(key):
                closed += 1
                if not path(i).exists():
                    search.save(path(i),search.encoded(closed_receipt(key)))
        missing = []
        for i, key in enumerate(keys):
            if path(i).exists():
                item = json.loads(path(i).read_bytes())
                check_generation(i,key,item,supplied)
                search.check_receipt(item,RANK,key)
            else:
                missing.append(i)
        print(f'{len(keys)-len(missing)}/{len(keys)} checked shards ({len(supplied)} frozen, '
              f'{closed} covered by the complete fourteen-fours root). {len(missing)} outstanding.',flush=True)
        if args.prepare_only:
            return
        if args.max_jobs is not None:
            missing = missing[:args.max_jobs]
        def one(i):
            item = make_receipt(i,keys[i],first,second,supplied)
            search.check_receipt(item, RANK, keys[i])
            search.save(path(i), search.encoded(item))
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending = [pool.submit(one, i) for i in missing]
            for n, future in enumerate(as_completed(pending), 1):
                future.result()
                if n % 100 == 0 or n == len(pending):
                    print(f'Finished {n}/{len(pending)} scheduled jobs.', flush=True)
        if any(not path(i).exists() for i in range(len(keys))):
            print('Partial checkpoint saved. Run the same command to continue.', flush=True)
            return
        target = work / 'rank31_completed.zip'
        part = target.with_suffix('.partial')
        with ZipFile(part, 'w', ZIP_DEFLATED) as z:
            z.writestr('spec.json', search.encoded(expected))
            with z.open('receipts.jsonl', 'w', force_zip64=True) as records:
                def stream():
                    for i,key in enumerate(keys):
                        raw = path(i).read_bytes()
                        item = json.loads(raw)
                        check_generation(i,key,item,supplied)
                        records.write(search.encoded(item) + b'\n')
                        yield key, item
                report = summary(stream())
            z.writestr('screen.json', search.encoded(report))
        part.replace(target)
        print(f"{report['profiles']} profiles, {len(report['survivors'])} survivors. Return {target}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--work', type=Path, default=Path('rank-work/31'))
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--max-jobs', type=int, help='process at most this many outstanding shards and save progress')
    p.add_argument('--verify', type=Path, metavar='ZIP')
    p.add_argument('--replay', action='store_true', help='repeat the searches while verifying a completed archive')
    args = p.parse_args()
    if args.jobs < 1 or args.max_jobs is not None and args.max_jobs < 1 or args.replay and not args.verify:
        p.error('positive job counts are required; replay requires --verify')
    with TemporaryDirectory() as tmp:
        first, hashes = build(Path(tmp))
        second, second_hashes = build_second(Path(tmp))
        keys = plan(first)
        first_spec = spec(hashes, keys)
        supplied, seed_digest = read_seed(first_spec, keys)
        expected = {**first_spec, 'format':'adaptive-31-v4', 'seed_sha256':seed_digest,
                    'remaining_sources':second_hashes,
                    'fourteen_ledger_sha256':search.digest((fourteen.DATA/'fourteen_31.json').read_bytes())}
        if args.replay:
            fourteen.replay(second,second_hashes)
        else:
            fourteen.stored(second_hashes)
        search.check_nonsquare(RANK)
        if args.verify:
            report = verify(args.verify, first, second, keys, expected, supplied, args.replay)
            print(f"Verified rank 31: {len(keys)} shards, {report['profiles']} profiles, {len(report['survivors'])} survivors.")
        else:
            run(args, first, second, keys, first_spec, expected, supplied)


if __name__ == '__main__':
    main()
