"""Project the independently returned rank-28 list through the frozen new root tests."""
from pathlib import Path
from zipfile import ZipFile
from math import isqrt
from collections import Counter
import json,sys,hashlib
ROOT=Path(__file__).resolve().parent

def admitted(ds):
    D=sum(x*x for x in ds);f=D//ds[-1]**2
    # q is squarefree and f=q*b1^2. The small f avoids factoring D.
    q=next(f//(b*b) for b in range(isqrt(f),0,-1) if f%(b*b)==0);M=isqrt(D//q)
    assert q*M*M==D
    counts=Counter(M//d for d in ds[1:]);m,j,k=counts[4],counts[5],counts[6]
    if f==9:return 'codegree nine'
    if f==16 and m==14 and j>=2:return 'fourteen and at least two'
    if f==16 and (m,j,k)==(14,1,3):return 'fourteen, one, three'
    if f==16 and (m,j) in ((1,23),(2,21),(3,20),(5,17),(7,14),(10,9),(12,6),(14,3),(15,1)):return 'near-full pair'
    A=6 if not q%3 else (4 if not q%2 else 0)
    if A:
        for p in range(3,2*len(ds)+2,2):
            if D%p or q%p==0 or any(p%i==0 for i in range(2,isqrt(p)+1)):continue
            h=(p-1)//2
            # All occurrences of an equal denominator form one eligible run.
            if sum(m//h for b,m in counts.items() if not b%p)<A:return 'complementary packets'
    return 'kept'

if __name__=='__main__':
    with ZipFile(sys.argv[1]) as ar: raw=ar.read('candidates.jsonl').splitlines()
    counter=Counter(); good=[]
    for w in raw:
        reason=admitted(json.loads(w));counter[reason]+=1
        if reason=='kept':good.append(w+b'\n')
    print(dict(counter),hashlib.sha256(b''.join(good)).hexdigest())
