"""Compare legacy native, revised native/GMP and an independent projection of archived raw rows."""
import argparse
from collections import Counter
from math import isqrt
from pathlib import Path
import subprocess
import json
import tempfile
import time
import zipfile

# These recorded rank-30 coordinates span the three most expensive old patterns,
# ones eliminated during generation, and two large mixed-prime terminal lists.
COORDINATES = [
    [4]*14 + [7]*3 + [7,7,9],
    [4]*14 + [6]*3 + [8,12,12],
    [4]*12 + [5]*4 + [6,6,8,15],
    [4]*12 + [5]*4 + [6,6,6,14],
    [4]*12 + [5]*4 + [6,6,6,18],
    [4]*12 + [5]*5 + [9,10,10],
]


def prime_divisors(n):
    p = 2
    out = []
    while p*p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def prefix_condition(ds, prefix, h):
    """Independent integer projection for these square-class-one rank-30 fixtures."""
    counts = Counter(prefix)
    if counts[4] == 14 and counts[7]:
        return False
    if counts[4] == 13 and counts[5] >= 4 or counts[4] == 11 and counts[5] >= 6:
        return False
    last=prefix[-1]
    if last>5:
        j=counts[5]
        if j==1 and len(ds)<23+counts[4]:
            return False
        if j in (2,3,6,7):
            scale=1<<20
            need=3 if j==3 else 2
            outside=sum(scale//(b*b) for b in prefix if b%5)
            known=sum((scale+b*b-1)//(b*b) for b in prefix if b%5==0 and b!=5)
            next_b=5*((last+4)//5)
            cap=known+h*((scale+next_b**2-1)//next_b**2)+1
            if 25*outside+(j+need)*scale>25*scale or 25*cap<need*scale:
                return False
    required=0
    for p in {p for b in counts for p in prime_divisors(b) if p>2}:
        step=(p-1)//2
        done=sum((m//step)*step for b,m in counts.items() if b%p==0)
        if last%p==0:
            done+=counts[last]%step
        required+=max(0,p-1-done)
    return required<=h


def accepted(row):
    M=isqrt(sum(d*d for d in row))
    bs=[M//d for d in row[:0:-1]]
    return all(prefix_condition(row, bs[:j], len(bs)-j)
               for j in list(range(1,len(bs)-1))+[len(bs)])


def run(args):
    started=time.monotonic()
    p=subprocess.run(args,check=True,capture_output=True,text=True,timeout=60)
    rows=sorted(tuple(json.loads(x)) for x in p.stdout.splitlines())
    if len(rows)!=len(set(rows)) or 'GMP ' in p.stderr or not ('Complete rank ' in p.stderr or 'COMPLETE search' in p.stderr):
        raise ValueError(p.stderr)
    return rows,time.monotonic()-started


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    default_archive=Path(__file__).resolve().parents[1]/'data'/'revision'/'rank30_completed.zip'
    parser.add_argument('archive',nargs='?',type=Path,default=default_archive)
    parser.add_argument('--quick',action='store_true',help='test the first two shard coordinates')
    parser.add_argument('--baseline',type=Path,help='optional legacy fast.cpp executable')
    parser.add_argument('--include-slow-baseline',action='store_true')
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    receipts=list(map(json.loads,zipfile.ZipFile(args.archive).read('receipts.jsonl').splitlines()))
    shards={tuple(z['prefix']):z for r in receipts for z in r.get('children',[]) if z['q']==1}
    results=[]
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        for name,macro,libs in [('fast',False,[]),('exact',False,['-lgmpxx','-lgmp']),('fast',True,[])]:
            options=['-DREPLAY_ARITHMETIC_ONLY'] if macro else []
            output=tmp/('arithmetic' if macro else name)
            subprocess.run(['c++','-O3','-std=c++17',*options,str(here/(name+'.cpp')),*libs,'-o',str(output)],check=True)
        for bs in (COORDINATES[:2] if args.quick else COORDINATES):
            receipt=shards[tuple(bs)]
            expected=sorted(tuple(row) for row in receipt['rows'])
            narrowed=list(filter(accepted,expected))
            key=','.join(map(str,bs))
            settings=[('arithmetic',tmp/'arithmetic',expected),('fast',tmp/'fast',narrowed),
                      ('exact',tmp/'exact',narrowed)]
            if args.baseline and (receipt['native_nodes']<4_000_000 or args.include_slow_baseline):
                settings.append(('legacy',args.baseline,expected))
            out={'prefix':bs,'old_nodes':receipt['native_nodes'],'archived_types':len(expected),
                 'new_types':len(narrowed),'timings':{}}
            for name,exe,answer in settings:
                command=[str(exe),'--shard','30','1',key] if name=='exact' else [str(exe),'30','8','1',key]
                rows,elapsed=run(command)
                if rows!=answer:
                    raise ValueError((name,key,len(rows),len(answer)))
                out['timings'][name]=round(elapsed,5)
            results.append(out)
    data=json.dumps(results,indent=2)+'\n'
    if args.out:
        args.out.write_text(data)
    print(data,end='')

if __name__=='__main__':
    main()
