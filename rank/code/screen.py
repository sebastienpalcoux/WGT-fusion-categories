"""Check the dimension identity, ordered degrees, prime families, and trace rows."""
import argparse
import json
from collections import Counter
from pathlib import Path
from families import dimension_family, factors
from trace_rows import nonunit_trace_row


def screen(degrees, ordinary=True):
    ds=tuple(degrees)
    if not ds or ds!=tuple(sorted(ds)) or ds.count(1)!=1:
        raise ValueError('expected sorted degrees with exactly one unit')
    total=sum(x*x for x in ds)
    if any(total%(x*x) for x in ds):
        raise ValueError('failed half-Frobenius divisibility')
    if dimension_family(total if ordinary else 2*total):
        return 'prime family'
    for pos,x in enumerate(sorted(set(ds[1:]),reverse=True),1):
        # This cutoff can only keep additional candidates.
        if total//(x*x)>120:
            break
        if not nonunit_trace_row(ds,x,ordinary=ordinary):
            return f'trace at degree {pos}'
    return 'open'


def run(rows,ordinary=True):
    counts=Counter()
    survivors=[]
    for d in rows:
        status=screen(d,ordinary)
        counts[status]+=1
        if status=='open':survivors.append(d)
    return {'counts':dict(sorted(counts.items())), 'survivors':survivors}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('--jsonl',action='store_true')
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    if args.jsonl:
        rows=(json.loads(line) for line in args.source.read_text().splitlines())
    else:
        rows=json.loads(args.source.read_text())['types']
    record=run(rows)
    if args.output: args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))
