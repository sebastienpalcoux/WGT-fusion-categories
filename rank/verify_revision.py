#!/usr/bin/env python3
"""Validate the returned archive, revised range certificates and the optimized replay."""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data/revision'
sys.path.insert(0,str(ROOT/'code'))
from check_family28 import check
from compare28 import admitted
from fiber_packets import inspect
from mixed_row import test as mixed_row
from nine_two import solutions
from nine_signed_check import verify as nine_crosscheck
from regular_indicator import obstructions
from next_rank import executables, manifest, check_receipt, verify as verify_current
from current_screen import screen


def require(condition, name):
    if not condition:
        raise ValueError(name)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replay',action='store_true',help='also regenerate optimized and nonsquare shards')
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args()
    original=ROOT/'data/ordinary/rank28_completed.zip'
    subprocess.run([sys.executable,ROOT/'finish_rank28.py','--verify',original],check=True)
    check(original)
    known=json.loads((DATA/'rank28_proof.json').read_bytes())
    with ZipFile(original) as z:
        rows=[json.loads(l) for l in z.read('candidates.jsonl').splitlines()]
        survivors=json.loads(z.read('family_result.json'))['survivors']
    require([v['degrees'] for v in known['proof']]==survivors,'residual lists disagree')
    require((len(rows),len(survivors))==(9458,66),'historical counts disagree')
    records=[inspect(ds) for ds in survivors]
    require(records==json.loads((DATA/'rank28_fiber_packets.json').read_bytes()),'fiber records disagree')
    pairs=Counter()
    for ds,rec,v in zip(survivors,records,known['proof']):
        D=sum(d*d for d in ds)
        g=min(d for d in ds if d>1)
        cap=sum(d*d for d in ds if d>1 and D//(d*d)%5==0)
        eligible=rec['primes']['5']
        require((D,g,cap,eligible['available'])==(v['D'],v['g'],v['U5'],v['pairs5']),'mass record')
        require(10*g*cap<4*D and eligible['required']==6,'five-active proof')
        pairs[eligible['available']]+=1
    require(pairs=={2:4,3:59,4:3},'packet capacity histogram')
    require(sum(bool(obstructions(x)) for x in rows)==known['regular_indicator_all']==1494,'regular indicator')
    require(sum(bool(obstructions(x)) for x in survivors)==known['regular_indicator_residual']==35,'indicator residuals')
    require(Counter(map(mixed_row,survivors))=={'orthogonality':3,'modular inequality':63},'mixed rational rows')
    require(sorted({n for n,_ in solutions()})==nine_crosscheck()==[70,110,140,220,380],'signed Egyptian sums')
    subprocess.run([sys.executable,ROOT/'code/check_adjacent.py'],check=True)
    with ZipFile(DATA/'rank28_optimized.zip') as z:
        optimized=sorted(row for item in z.read('receipts.jsonl').splitlines() for row in json.loads(item)['rows'])
    require(optimized==[r for r in rows if admitted(r)=='kept'] and len(optimized)==1940,
            'optimized generator and independent historical projection disagree')
    with TemporaryDirectory() as temp:
        tmp=Path(temp)
        exe,hashes=executables(tmp)
        for path in sorted(DATA.glob('rank*_optimized.zip')):
            verify_current(path,exe,hashes,a.replay,a.jobs)
        with ZipFile(DATA/'rank30_checkpoint.zip') as z:
            info=json.loads(z.read('spec.json'))
            require(info=={'depth':11,'format':1,'rank':30,'sources':hashes},'checkpoint sources')
            keys=set(manifest(exe['exact'],30,11))
            results=[]
            for line in z.read('receipts.jsonl').splitlines():
                item=json.loads(line)
                key=(item['q'],tuple(item['prefix']))
                require(key in keys,'duplicate checkpoint shard or wrong coordinate')
                keys.remove(key)
                results.extend(check_receipt(item,30,key))
            pending={(q,tuple(b)) for q,b in json.loads(z.read('pending.json'))}
            require(keys==pending and all(q==1 and b[0]==4 for q,b in keys),'wrong remaining search')
            require(len(results)==len(set(results)),'duplicate rank-30 types')
            partial=screen(sorted(results),30)
            require(partial==json.loads(z.read('screen.json')) and not partial['survivors'],'rank-30 tests')
            print(f'Rank 30 checkpoint: {len(results)} excluded degree types, {len(keys)} pending shards.')
        odd_jobs=subprocess.check_output([sys.executable,ROOT/'code/odd_jobs.py','57'])
        require(json.loads(odd_jobs)==json.loads((DATA/'odd57_jobs.json').read_bytes()),'rank-57 support list')
        for ref in sorted(DATA.glob('*.json')):
            if ref.name.startswith('odd') and ref.stem[3:].isdigit():
                args=[str(int(ref.stem[3:])),str(tmp/'odd.json')]
            elif ref.name.startswith('nonsquare') and ref.stem[9:].isdigit():
                args=[str(int(ref.stem[9:])),str(tmp/'odd.json'),'--omit-square-class-one']
            else:
                continue
            subprocess.run([sys.executable,ROOT/'code/odd_frontier.py',*args],check=True)
            require(json.loads((tmp/'odd.json').read_bytes())==json.loads(ref.read_bytes()),ref.name)
    subprocess.run([sys.executable,ROOT/'verify_nonsquare.py',*(['--replay'] if a.replay else [])],check=True)
    print('All completed rank and local arithmetic certificates passed.')


if __name__=='__main__':
    main()
