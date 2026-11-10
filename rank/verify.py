#!/usr/bin/env python3
"""Replay the finite certificates. Use --full to independently regenerate all lists."""
from argparse import ArgumentParser
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent
CODE, DATA = ROOT/'code', ROOT/'data'
sys.path.insert(0,str(CODE))
from families import dimension_family, factors
from screen import run as screen_rows
from check_reduced import check as reduced_verdict
from trace_rows import nonunit_trace_row


def require(ok,message):
    if not ok:
        raise ValueError(message)


def load(name):
    return json.loads((DATA/name).read_text())


def run(cmd):
    return subprocess.run(list(map(str,cmd)),check=True,capture_output=True,text=True)


def compile_tool(name,work):
    exe=work/name
    run(['c++','-O3','-std=c++17',CODE/(name+'.cpp'),'-o',exe])
    return exe


def compare(source,name):
    require(source==(DATA/name).read_bytes(),f'Regenerated output differs: {name}')


def optional_packets(profile):
    # Independent audit uses only the unordered four-prime tests and absolute
    # threshold, matching the earlier of the two enumerations in the paper.
    P=factors(profile)
    E=[e for _,e in P]
    if len(P)<3 or profile<29160000: return True
    if len(P)==3:
        return E[0]<=5 or E[1]<=4 or E[2]<=3 or all(x<=y for x,y in zip(sorted(E),(4,5,7)))
    if len(P)==4:
        d,c,b,a=sorted(E)
        return b<=3 or d==2 and c<=3 and (a//2+b//2<=5 or P[0][0]==2 and E[0]>=4 and a<=9 and b<=5)
    if len(P)==5:
        d,c,b,s,a=sorted(E)
        return a<=3 or a==4 and s<=3 and d==2
    return False


def check_lists():
    require(not load('root18.json')['types'],'unexpected rank <=18 profile')
    ordinary=load('root19.json')['types']
    require(len(ordinary)==41 and all(len(ds)==19 for ds in ordinary),'ordinary profiles')
    verdicts=screen_rows(ordinary)
    require(verdicts==load('ordinary19_verdicts.json') and not verdicts['survivors'],'ordinary trace check')
    counts=Counter()
    for q,M,ds in load('single_packet19.json')['row_feasible']:
        if len(ds)!=18 or optional_packets(q*M*M): continue
        for i,x in enumerate(sorted(set(ds),reverse=True)[:2]):
            if not nonunit_trace_row(sorted([1]+ds),x): counts[i]+=1;break
        else: raise ValueError('independent trace survivor')
    require(counts=={0:595,1:3},'independent trace counts')
    for q,M,ds in load('reduced14.json')['row_feasible']:
        require(1+sum(d*d for d in ds)==q*M*M and 2*q*M*M<29160000,'small reduced profile')
    require(len(load('reduced14.json')['row_feasible'])==8,'small reduced count')
    for rank,expected in [(15,{'BDP':6}),(16,{'BDP':133,'row gcd':104,'four':18,'three':1})]:
        rows=[json.loads(line) for line in (DATA/f'four{rank}.jsonl').read_text().splitlines()]
        require(Counter(map(reduced_verdict,rows))==expected,f'reduced rank {rank}')
    audit=load('odd39_audit.json')['row_feasible']
    require(len(audit)==82 and all(q*M*M==1+sum(x*x for x in ds) and q*M*M//max(ds)**2==7 for q,M,ds in audit),'odd audit')
    require(not load('odd39.json')['row_feasible'],'odd prime screen')


def regenerated_py(name,script,args,work,options=()):
    target=work/name
    run([sys.executable,CODE/script,*args,target,*options])
    compare(target.read_bytes(),name)


def main():
    p=ArgumentParser(description=__doc__); p.add_argument('--full',action='store_true');args=p.parse_args()
    # Reading the checksums here detects damaged or accidentally changed output files.
    for line in (DATA/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1)
        require(hashlib.sha256((DATA/name).read_bytes()).hexdigest()==expected,name)
    check_lists()
    print('All stored ordinary, odd and reduced type lists passed their independent elimination checks.',flush=True)
    with TemporaryDirectory() as tmp:
        work=Path(tmp); formal=compile_tool('formal_sieve',work)
        compare(run([formal]).stdout.encode(),'formal_residual.txt')
        run([sys.executable,CODE/'formal_check.py',DATA/'formal_residual.txt'])
        compare(run([formal,'--odd13']).stdout.encode(),'odd13.txt')
        compare(run([formal,'--odd15']).stdout.encode(),'odd15.txt')
        compare(run([sys.executable,CODE/'odd_formal.py']).stdout.encode(),'odd_formal.txt')
        regenerated_py('odd39.json','odd_rank.py',[39],work)
        frontier=compile_tool('frontier',work)
        log=[]
        for q in (3,10,14,15):
            proc=run([frontier,20,8,q]); final=proc.stderr.splitlines()[-1]
            require(proc.stdout=='' and final.startswith('Complete rank 20,'),'incomplete kernel')
            nodes=int(final.split('nodes ')[1].split(',')[0]);log.append([q,nodes])
        require(log==[[3,814],[10,276],[14,17],[15,17]],'rank 20 easy kernels')
        if args.full:
            for r in (18,19):
                regenerated_py(f'root{r}.json','squares.py',[r],work,['--strong-packets'])
            regenerated_py('single_packet19.json','single_packet.py',[19],work,['--min-codegree','6'])
            regenerated_py('odd39_audit.json','odd_rank.py',[39],work,['--audit-prime-codegrees'])
            regenerated_py('reduced14.json','reduced_types.py',[14],work)
            reduced=compile_tool('reduced_four',work)
            for r in (15,16):
                result=run([reduced,r]);require('Complete rank' in result.stderr,'missing completion')
                compare(result.stdout.encode(),f'four{r}.jsonl')
    replay = ['--replay'] if args.full else []
    for cmd in ([sys.executable,ROOT/'continuation/verify_result.py',DATA/'ordinary/rank20_result.zip',*replay],
                [sys.executable,ROOT/'verify_extended.py',*replay]):
        subprocess.run([str(part) for part in cmd],check=True)
    subprocess.run([sys.executable, str(ROOT/'verify_frontier.py'), *replay], check=True)
    subprocess.run([sys.executable, str(ROOT/'verify_revision.py'), *replay], check=True)
    print('PASS. All archived exclusions and current rank certificates verified.')


if __name__=='__main__': main()
