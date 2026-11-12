"""Shared exact arithmetic for the conditional square searches (R04)."""
import hashlib
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'code'))
from current_screen import screen


def require(condition, message):
    if not condition:
        raise ValueError(message)


def degree_vector(line, rank):
    d = json.loads(line)
    require(isinstance(d, list) and len(d) == rank and all(type(x) is int and x > 0 for x in d),
            "malformed dimension vector")
    require(d == sorted(d) and d.count(1) == 1, "invalid type order or multiplicity of the unit")
    total = sum(x*x for x in d)
    require(all(total % (x*x) == 0 for x in d), "square divisibility failed")
    return tuple(d)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return json.dumps(value, separators=(',', ':'), sort_keys=True).encode('ascii')


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


def compute(exe, rank, key):
    q,prefix = key
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
    require('rule' not in body, 'conditional square receipts must be computed')
    rows=emit_rows(map(encoded,body['rows']),rank,key)
    require(rows==sorted(set(rows)),'duplicate or unordered shard output')
    return rows


if __name__ == "__main__":
    raise SystemExit(
        "Use python3 papers/rank-bound-supplement/verify.py --rank R for a new search; "
        "this module supplies the conditional square-suite helpers."
    )
