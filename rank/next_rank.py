#!/usr/bin/env python3
"""Resumable exact rank search. For example, python3 next_rank.py 32 --jobs 12."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import fcntl
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'continuation'))
from verify_result import degree_vector, require
from current_screen import screen


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return json.dumps(value, separators=(',', ':'), sort_keys=True).encode('ascii')


def save(path, raw):
    temporary = path.with_suffix('.partial')
    temporary.write_bytes(raw)
    temporary.replace(path)


@contextmanager
def locked(work):
    work.mkdir(parents=True, exist_ok=True)
    with (work/'run.lock').open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def executables(tmp):
    result = {}
    fingerprints = {}
    for name, libs in [('fast', []), ('exact', ['-lgmpxx', '-lgmp'])]:
        source = ROOT/'current'/(name+'.cpp')
        fingerprints[name] = digest(source.read_bytes())
        exe = tmp/name
        subprocess.run(['c++','-O3','-std=c++17',str(source), *libs, '-o',str(exe)],check=True)
        result[name] = exe
    return result, fingerprints


def prefix_text(prefix):
    return ','.join(map(str,prefix))


def exact_run(exe, mode, rank, args):
    task = subprocess.run([exe,'--'+mode,str(rank),*args],capture_output=True,check=True)
    tag = 'manifest' if mode=='manifest' else 'search'
    result = re.fullmatch(rf'COMPLETE {tag} rank={rank} nodes=(\d+) types=(\d+) shards=(\d+)\n',
                          task.stderr.decode('ascii'))
    require(result is not None, 'missing GMP completion marker')
    nodes,nrows,nshards = map(int,result.groups())
    lines = task.stdout.splitlines()
    if mode=='manifest':
        require(nrows==0 and nshards==len(lines),'GMP manifest count mismatch')
    else:
        require(nrows==len(lines) and nshards==0,'GMP row count mismatch')
    return nodes, lines


def uniform_manifest(exe,rank,depth):
    _, lines = exact_run(exe,'manifest',rank,[str(depth)])
    keys=[]
    for line in lines:
        q,values=line.split(b'\t')
        key=(int(q),tuple(map(int,values.split(b','))))
        require(len(key[1])==depth and sorted(key[1])==list(key[1]),'malformed partition')
        keys.append(key)
    require(len(set(keys))==len(keys),'duplicate manifest shard')
    return keys


def manifest(exe, rank, depth):
    base = uniform_manifest(exe, rank, depth)
    hot = (1, (4,)*depth)
    if depth >= 17 or hot not in base:
        return base
    # Refine the expensive, all-four parent. The GMP manifest defines both partitions.
    children = uniform_manifest(exe, rank, min(17, rank-3))
    return [key for key in base if key != hot] + [
        key for key in children if key[0] == 1 and key[1][:depth] == hot[1]]


def emit_rows(raw,rank,key):
    q,prefix = key
    rows=[]
    for line in raw:
        row=degree_vector(line,rank)
        dimension=sum(d*d for d in row)
        M=isqrt(dimension//q)
        require(M*M*q==dimension and all(M%d==0 for d in row),'wrong square class')
        require(tuple(M//d for d in reversed(row[1:]))[:len(prefix)]==prefix,'row outside shard')
        rows.append(row)
    return rows


def gauss_shortcut(rank, key):
    q, prefix = key
    if q > 1 and rank in (30,31):
        return 'independent complete nonsquare enumeration'
    if q != 1 or prefix[0] != 4:
        return None
    if prefix[-1] > 5 and prefix.count(5) == 1 and rank < 23 + prefix.count(4):
        return 'unique quintic denominator rank bound'
    for i,j in [(13,4),(11,7),(11,6)]:
        if prefix.count(4) == i and prefix.count(5) >= j:
            return f'quintic Gauss row {i},{j}'
    # Once the prefix has passed five the multiplicity at five is final.
    j=prefix.count(5)
    if prefix[-1] > 5 and j in (2,3,6,7):
        cutoff=Fraction(3 if j==3 else 2,25)
        complement=1-Fraction(j,25)-sum(Fraction(1,b*b) for b in prefix if b%5)
        known=sum(Fraction(1,b*b) for b in prefix if b%5==0 and b!=5)
        next_multiple=5*((prefix[-1]+4)//5)
        slots=rank-1-len(prefix)
        from_tail=known+Fraction(slots,next_multiple**2)+Fraction(1,5400**2)
        if min(complement,from_tail)<cutoff:
            return f'quintic support mass, multiplicity {j}'
    return None


def compute(exe, rank, key, use_theory=True):
    q,prefix = key
    rule = gauss_shortcut(rank,key) if use_theory else None
    if rule is not None:
        body = {'q':q,'prefix':list(prefix),'native_nodes':0,'gmp_nodes':0,'fallbacks':0,'rows':[],'rule':rule}
        return {**body,'sha256':digest(encoded(body))}
    # Exceeding the native range is allowed only through an explicit handoff or a full GMP retry.
    try:
        if any(b >= 1 << 64 for b in prefix):
            raise OverflowError
        native=subprocess.run([exe['fast'],str(rank),'8',str(q),prefix_text(prefix)],
                              capture_output=True,check=True)
    except (OverflowError,subprocess.CalledProcessError):
        nodes,raw=exact_run(exe['exact'],'shard',rank,[str(q),prefix_text(prefix)])
        rows=emit_rows(raw,rank,key)
        return receipt(key,0,nodes,1,rows)
    report=native.stderr.decode('ascii').splitlines()
    marker=re.fullmatch(rf'Complete rank {rank}, nodes (\d+), types (\d+), seconds [\d.eE+\-]+',report[-1])
    require(marker is not None,'native completion marker missing')
    nodes,count=map(int,marker.groups())
    rows=emit_rows(native.stdout.splitlines(),rank,key)
    require(len(rows)==count,'native row count disagrees')
    delegated=[]
    for line in report[:-1]:
        if re.fullmatch(rf'q={q} cumulative nodes \d+ types \d+',line):
            continue
        match=re.fullmatch(rf'GMP {q} (\d+(?:,\d+)*)',line)
        require(match is not None,'unexpected native diagnostic: '+line)
        child=tuple(map(int,match[1].split(',')))
        require(child[:len(prefix)]==prefix,'GMP handoff outside shard')
        delegated.append(child)
    ordered=sorted(delegated)
    for first,second in zip(ordered,ordered[1:]):
        require(second[:len(first)]!=first,'overlapping GMP handoffs')
    exact_nodes=0
    for child in ordered:
        work,raw=exact_run(exe['exact'],'shard',rank,[str(q),prefix_text(child)])
        rows.extend(emit_rows(raw,rank,(q,child)))
        exact_nodes+=work
    return receipt(key,nodes,exact_nodes,len(ordered),rows)


def receipt(key,native,gmp,jobs,rows):
    require(len(set(rows))==len(rows),'duplicate degree vector')
    body={'q':key[0],'prefix':list(key[1]),'native_nodes':native,'gmp_nodes':gmp,
          'fallbacks':jobs,'rows':sorted(rows)}
    return {**body,'sha256':digest(encoded(body))}


def check_receipt(item,rank,key):
    body={k:v for k,v in item.items() if k!='sha256'}
    require(digest(encoded(body))==item['sha256'],'shard hash disagrees')
    require((body['q'],tuple(body['prefix']))==key,'wrong shard coordinate')
    require(all(type(body[k]) is int and body[k]>=0 for k in ('gmp_nodes','native_nodes','fallbacks')),
            'bad search counts')
    if 'rule' in body:
        require(body['rule']==gauss_shortcut(rank,key) and not body['rows'] and
                all(body[k]==0 for k in ('gmp_nodes','native_nodes','fallbacks')), 'invalid theoretical receipt')
    rows=emit_rows(map(encoded,body['rows']),rank,key)
    require(rows==sorted(set(rows)),'duplicate or unordered shard output')
    return rows


def evaluate(receipts,rank,keys):
    require(len(receipts)==len(keys),'incomplete receipt set')
    rows=sorted(row for item,key in zip(receipts,keys) for row in check_receipt(item,rank,key))
    require(len(set(rows))==len(rows),'overlapping shards')
    return screen(rows,rank)


def check_nonsquare(rank):
    if rank in (30,31):
        subprocess.run([sys.executable, ROOT/'verify_nonsquare.py','--ranks',str(rank)],check=True)


def seed_receipts(work, rank, depth, keys, spec):
    path=ROOT/'data/revision'/f'rank{rank}_checkpoint.zip'
    if not path.exists():
        return
    with ZipFile(path) as z:
        old=json.loads(z.read('spec.json'))
        if old != spec or old['depth']!=depth:
            return
        available={key:i for i,key in enumerate(keys)}
        count=0
        for line in z.read('receipts.jsonl').splitlines():
            item=json.loads(line)
            key=(item['q'],tuple(item['prefix']))
            require(key in available,'seed outside manifest')
            index=available.pop(key)
            check_receipt(item,rank,key)
            dest=work/f'shard-{index:07}.json'
            if not dest.exists():
                save(dest,encoded(item))
                count+=1
        if count:
            print(f'Restored {count} completed shards.',flush=True)


def run_rank(rank,depth,jobs,work,exe,fingerprints):
    spec={'format':1,'rank':rank,'depth':depth,'sources':fingerprints}
    check_nonsquare(rank)
    with locked(work):
        if (work/'spec.json').exists():
            require((work/'spec.json').read_bytes()==encoded(spec),'different rank, depth or source in work directory')
        else:
            save(work/'spec.json',encoded(spec))
        keys=manifest(exe['exact'],rank,depth)
        seed_receipts(work,rank,depth,keys,spec)
        def one(pair):
            i,key=pair
            path=work/f'shard-{i:07}.json'
            if path.exists():
                item=json.loads(path.read_bytes())
                check_receipt(item,rank,key)
            else:
                item=compute(exe,rank,key)
                save(path,encoded(item))
            return item
        print(f'Rank {rank}, {len(keys)} disjoint shards, {jobs} workers.',flush=True)
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            receipts=list(pool.map(one,enumerate(keys)))
        result=evaluate(receipts,rank,keys)
        target=work/f'rank{rank}_completed.zip'
        temporary=target.with_suffix('.partial')
        with ZipFile(temporary,'w',ZIP_DEFLATED) as z:
            z.writestr('spec.json',encoded(spec))
            z.writestr('receipts.jsonl',b''.join(encoded(r)+b'\n' for r in receipts))
            z.writestr('screen.json',encoded(result))
        temporary.replace(target)
        print(f"{len(keys)} shards, {result['profiles']} profiles, {len(result['survivors'])} survivors. {target}",flush=True)


def verify(path,exe,fingerprints,replay=False,jobs=1):
    with ZipFile(path) as z:
        require(sorted(z.namelist())==['receipts.jsonl','screen.json','spec.json'],'unexpected archive members')
        spec=json.loads(z.read('spec.json'))
        require(spec['sources']==fingerprints and spec['format']==1,'unrecognized generator source or format')
        rank=spec['rank']
        check_nonsquare(rank)
        keys=manifest(exe['exact'],rank,spec['depth'])
        receipts=[json.loads(line) for line in z.read('receipts.jsonl').splitlines()]
        result=evaluate(receipts,rank,keys)
        require(result==json.loads(z.read('screen.json')),'current screening disagrees')
        plain=[]
        split=[]
        for key,item in zip(keys,receipts):
            if any(k in item for k in ('children','split_depth','splitter_sha256')):
                split.append((key,item))
            else:
                plain.append((key,item))
        if split:
            from finish_frontier import audit_split, build_splitter
            splitter=build_splitter(exe['exact'].parent)
            for key,item in split:
                audit_split(splitter,exe,rank,key,item,replay,jobs)
        if replay:
            def redo(pair):
                key,old=pair
                new=compute(exe,rank,key,use_theory='rule' in old)
                require(encoded(old)==encoded(new),'shard replay disagrees')
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                list(pool.map(redo,plain))
        print(f"Verified rank {rank}: {len(keys)} shards, {result['profiles']} profiles, {len(result['survivors'])} survivors.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('ranks',nargs='*',type=int)
    group.add_argument('--verify',type=Path)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--depth',type=int,default=11)
    p.add_argument('--work',type=Path,default=Path('rank-work'))
    p.add_argument('--replay',action='store_true',help='regenerate each shard with --verify')
    a=p.parse_args()
    if a.jobs<1 or a.replay and not a.verify or any(not 24<=r<=40 or not 1<=a.depth<=r-3 for r in a.ranks):
        p.error('ranks must be 24..40, depth must be 1..rank-3, and job count positive')
    with TemporaryDirectory() as scratch:
        exe,fingerprints=executables(Path(scratch))
        if a.verify:
            verify(a.verify,exe,fingerprints,a.replay,a.jobs)
        else:
            for r in a.ranks:
                run_rank(r,a.depth,a.jobs,a.work.resolve()/str(r),exe,fingerprints)

if __name__=='__main__':
    main()
