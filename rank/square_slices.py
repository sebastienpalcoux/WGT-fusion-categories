#!/usr/bin/env python3
"""Replay three finite frontiers for square dimension D=M²."""
import argparse
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import next_rank as nr

ROOT = Path(__file__).resolve().parent
LEDGER = ROOT / 'data' / 'square-slices' / 'square_slices.json'


def replace_once(text, needle, addition):
    nr.require(text.count(needle) == 1, 'unexpected frozen source')
    return text.replace(needle, needle + '\n' + addition)


def valuation_engines(tmp):
    # Restrict the frozen q=1 recursion to M == 2 mod 4. The prefix lcm L
    # divides M. At a terminal pair both kernels explicitly compute M=k*L.
    result = {}
    for name in ('fast', 'exact'):
        text = (ROOT/'current'/(name+'.cpp')).read_text()
        if name == 'fast':
            text = replace_once(text, 'const PacketState &st,Mask bmask, const Ternary& tern) {',
                                '    if (L % 4 == 0) return;')
            text = replace_once(text, 'Wide root=Wide(k)*L;', '            if (root % 4 != 2) return;')
        else:
            text = replace_once(text, 'const Packets& st, const Ternary& ternary) {',
                                '        if (L % 4 == 0) return;')
            text = replace_once(text, 'Z s = k*L, b = s/e, c = s/d;', '        if (s % 4 != 2) return;')
        source = tmp / ('v2_' + name + '.cpp')
        source.write_text(text)
        exe = tmp / ('v2_' + name)
        libs = [] if name == 'fast' else ['-lgmpxx', '-lgmp']
        subprocess.run(['c++', '-O3', '-std=c++17', source, *libs, '-o', exe], check=True)
        result[name] = exe
    return result


def record(rank, rows, counters, valuation):
    rows = sorted(rows)
    nr.require(len(rows) == len(set(rows)), 'duplicate row across roots')
    if valuation:
        nr.require(all(isqrt(sum(x*x for x in row)) % 4 == 2 for row in rows), 'wrong square class')
    screened = nr.screen(rows, rank)
    return {**screened, **counters, 'rows_sha256': sha256(nr.encoded(rows)).hexdigest()}


def calculate(tmp):
    regular, _ = nr.executables(tmp)
    v2 = valuation_engines(tmp)
    report = {}
    for name,lo,hi,least in [('least_5',30,35,5), ('root_2_mod_4',30,38,5), ('least_6',37,45,6)]:
        result = []
        for rank in range(lo, hi+1):
            rows = []
            counters = {'native_nodes':0,'gmp_nodes':0,'handoffs':0}
            for b in range(least, isqrt(rank-1)+1):
                key = (1,(b,))
                if name == 'least_6':
                    nodes, raw = nr.exact_run(regular['exact'], 'shard', rank, ['1', str(b)])
                    rows.extend(nr.emit_rows(raw, rank, key))
                    counters['gmp_nodes'] += nodes
                else:
                    engine = v2 if name == 'root_2_mod_4' else regular
                    result_b = nr.compute(engine, rank, key)
                    rows.extend(nr.check_receipt(result_b, rank, key))
                    for key, field in [('native_nodes','native_nodes'),('gmp_nodes','gmp_nodes'),('handoffs','fallbacks')]:
                        counters[key] += result_b[field]
            entry = record(rank, rows, counters, name == 'root_2_mod_4')
            nr.require(not entry['survivors'], f'unexpected survivor for {name} at {rank}')
            result.append(entry)
        report[name] = result
    return report


def check_stored(ledger):
    limits = {'least_5': (30, 35, 5), 'root_2_mod_4': (30, 38, 5), 'least_6': (37, 45, 6)}
    nr.require(set(ledger) == set(limits), 'unexpected slice names')
    with ZipFile(LEDGER.with_name('square_slice_rows.zip')) as saved:
        expected = []
        for name, (first, last, minimum) in limits.items():
            items = ledger[name]
            nr.require([item['rank'] for item in items] == list(range(first, last+1)), 'rank coverage mismatch')
            for item in items:
                rank = item['rank']
                filename = f'{name}/{rank}.jsonl'
                expected.append(filename)
                with saved.open(filename) as member:
                    rows = [nr.degree_vector(line, rank) for line in member]
                for row in rows:
                    dim = sum(d*d for d in row)
                    root = isqrt(dim)
                    nr.require(root*root == dim and root//row[-1] >= minimum, 'wrong square slice')
                counters = {key:item[key] for key in ('native_nodes', 'gmp_nodes', 'handoffs')}
                nr.require(rows == sorted(rows) and record(rank, rows, counters, name=='root_2_mod_4') == item,
                           'supplied row list or screen disagrees')
        nr.require(sorted(saved.namelist()) == sorted(expected), 'unexpected list archive contents')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',type=Path,help='write a new ledger instead of comparing to the bundled one')
    parser.add_argument('--stored', action='store_true', help='check stored lists and exclusions without re-enumerating')
    args = parser.parse_args()
    if args.stored:
        if args.write:
            parser.error('--stored cannot produce a new ledger')
        report = json.loads(LEDGER.read_text())
        check_stored(report)
    else:
        with TemporaryDirectory() as tmp:
            report = calculate(Path(tmp))
    if args.write:
        args.write.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    else:
        nr.require(report == json.loads(LEDGER.read_text()), 'square-slice ledger disagrees')
    for name,records in report.items():
        print(name, [(x['rank'], x['profiles']) for x in records], 'all excluded')


if __name__ == '__main__':
    main()
