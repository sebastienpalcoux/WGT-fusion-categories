#!/usr/bin/env python3
"""Resume rank 30, splitting only the outstanding checkpoint parents across workers."""
import argparse
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

import next_rank as nr

RANK = 30
BASE_DEPTH = 11
SPLITTER = nr.ROOT/'current/child_manifest.cpp'


def build_splitter(directory):
    binary = directory/'child_manifest'
    subprocess.run(['c++','-O3','-std=c++17',str(SPLITTER),'-lgmpxx','-lgmp','-o',str(binary)],check=True)
    return binary


def partition(exe, rank, key, depth):
    q,prefix = key
    nr.require(type(depth) is int and len(prefix) < depth <= rank-3,'invalid child depth')
    _, lines = nr.exact_run(exe,'manifest',rank,[str(depth),str(q),nr.prefix_text(prefix)])
    children=[]
    for line in lines:
        child_q,denoms=line.split(b'\t')
        b=tuple(map(int,denoms.split(b',')))
        nr.require(int(child_q)==q and len(b)==depth and b[:len(prefix)]==prefix and
                   tuple(sorted(b))==b, 'invalid child manifest coordinate')
        children.append((q,b))
    nr.require(len(set(children))==len(children),'duplicate child manifest coordinate')
    return children


def check_plain(item, rank, key):
    nr.require(not any(k in item for k in ('children','split_depth','splitter_sha256')),'nested split is not permitted')
    return nr.check_receipt(item,rank,key)


def aggregate(rank, key, depth, children, items):
    nr.require(len(children)==len(items), 'missing child receipt')
    rows=[]
    for child,item in zip(children,items):
        rows.extend(check_plain(item,rank,child))
    merged=nr.receipt(key,sum(x['native_nodes'] for x in items),sum(x['gmp_nodes'] for x in items),
                      sum(x['fallbacks'] for x in items),rows)
    del merged['sha256']
    merged.update(split_depth=depth,splitter_sha256=nr.digest(SPLITTER.read_bytes()),children=items)
    merged['sha256']=nr.digest(nr.encoded(merged))
    return merged


def audit_split(splitter, exe, rank, key, item, replay=False, jobs=1):
    nr.require('rule' not in item and item.get('splitter_sha256')==nr.digest(SPLITTER.read_bytes()),
               'unrecognized split receipt')
    children=partition(splitter,rank,key,item['split_depth'])
    items=item['children']
    nr.require(isinstance(items,list),'invalid children list')
    expected=aggregate(rank,key,item['split_depth'],children,items)
    nr.require(nr.encoded(expected)==nr.encoded(item),'subshard receipts disagree with their parent')
    if replay:
        def repeat(pair):
            child,old=pair
            new=nr.compute(exe,rank,child,use_theory='rule' in old)
            nr.require(nr.encoded(new)==nr.encoded(old),'subshard replay disagrees')
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(repeat,zip(children,items)))


def rolled_tasks(groups):
    # Interleave parents, so a single large parent cannot occupy the whole queue.
    queues=deque((p,iter(range(len(keys)))) for p,keys in groups.items())
    while queues:
        p,queue=queues.popleft()
        i=next(queue,None)
        if i is not None:
            yield p,i
            queues.append((p,queue))


def run(work, depth, jobs, prepare):
    with TemporaryDirectory() as tmp, nr.locked(work):
        tmp=Path(tmp)
        exe,hashes=nr.executables(tmp)
        splitter=build_splitter(tmp)
        spec={'format':1,'rank':RANK,'depth':BASE_DEPTH,'sources':hashes}
        config=work/'spec.json'
        if config.exists():
            nr.require(config.read_bytes()==nr.encoded(spec),'different rank or frozen sources in work directory')
        else:
            nr.save(config,nr.encoded(spec))
        keys=nr.manifest(exe['exact'],RANK,BASE_DEPTH)
        with ZipFile(nr.ROOT/'data/revision/rank30_checkpoint.zip') as seed:
            nr.require(json.loads(seed.read('spec.json'))==spec,'checkpoint sources do not match')
            supplied=[json.loads(row) for row in seed.read('receipts.jsonl').splitlines()]
            given={(x['q'],tuple(x['prefix'])) for x in supplied}
            pending={(q,tuple(b)) for q,b in json.loads(seed.read('pending.json'))}
        nr.require(len(keys)==524 and len(given)==len(supplied) and given.issubset(keys)
                   and pending==set(keys)-given, 'checkpoint partition disagrees')
        nr.check_nonsquare(RANK)
        nr.seed_receipts(work,RANK,BASE_DEPTH,keys,spec)
        groups={}; receipts={}; subdir=work/'subshards'/f'depth-{depth}'
        for i,key in enumerate(keys):
            file=work/f'shard-{i:07}.json'
            if file.exists():
                item=json.loads(file.read_bytes())
                if 'children' in item:
                    audit_split(splitter,exe,RANK,key,item)
                else:
                    check_plain(item,RANK,key)
                receipts[i]=item
            elif nr.gauss_shortcut(RANK,key):
                receipts[i]=nr.compute(exe,RANK,key)
                nr.save(file,nr.encoded(receipts[i]))
            else:
                groups[i]=partition(splitter,RANK,key,depth)
                (subdir/f'{i:07}').mkdir(parents=True,exist_ok=True)
        work_items=[]; children_by_parent=defaultdict(dict)
        for p,i in rolled_tasks(groups):
            dest=subdir/f'{p:07}'/f'{i:07}.json'
            key=groups[p][i]
            if dest.exists():
                item=json.loads(dest.read_bytes());check_plain(item,RANK,key)
                children_by_parent[p][i]=item
            else:
                work_items.append((p,i,dest,key))
        n=sum(map(len,groups.values()))
        remaining_theory=sum(nr.gauss_shortcut(RANK,k) is not None for _,_,_,k in work_items)
        print(f'{len(receipts)}/524 base receipts, {len(groups)} parents split into {n} child shards; '
              f'{n-len(work_items)} child receipts cached. {remaining_theory} pending children have a theorem shortcut.',flush=True)
        if prepare:return
        def calculate(x):
            p,i,dest,key=x
            item=nr.compute(exe,RANK,key)
            check_plain(item,RANK,key)
            nr.save(dest,nr.encoded(item))
            return p,i,item
        def finish_parent(p):
            values=[children_by_parent[p][i] for i in range(len(groups[p]))]
            item=aggregate(RANK,keys[p],depth,groups[p],values)
            nr.save(work/f'shard-{p:07}.json',nr.encoded(item))
            receipts[p]=item
            print(f'{len(receipts)}/524 base receipts complete.',flush=True)
        for p in groups:
            if len(children_by_parent[p])==len(groups[p]):finish_parent(p)
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            tasks=[pool.submit(calculate,x) for x in work_items]
            for future in as_completed(tasks):
                p,i,item=future.result();children_by_parent[p][i]=item
                if len(children_by_parent[p])==len(groups[p]):finish_parent(p)
        ordered=[receipts[i] for i in range(len(keys))]
        result=nr.evaluate(ordered,RANK,keys)
        out=work/'rank30_completed.zip'
        part=out.with_suffix('.partial')
        with ZipFile(part,'w',ZIP_DEFLATED) as archive:
            archive.writestr('spec.json',nr.encoded(spec))
            archive.writestr('receipts.jsonl',b''.join(nr.encoded(item)+b'\n' for item in ordered))
            archive.writestr('screen.json',nr.encoded(result))
        part.replace(out)
        print(f"{result['profiles']} profiles, {len(result['survivors'])} survivors. Return {out}",flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs',type=int,default=12)
    p.add_argument('--child-depth',type=int,default=20,help='absolute prefix depth for unfinished parents (18–27)')
    p.add_argument('--work',type=Path,default=Path('rank-work/30'),help='rank-30 work directory')
    p.add_argument('--prepare-only',action='store_true',help='check the seed and subshard partition without running searches')
    a=p.parse_args()
    if a.jobs<1 or not 18<=a.child_depth<=27:p.error('jobs must be positive and child depth must be 18..27')
    run(a.work.resolve(),a.child_depth,a.jobs,a.prepare_only)
