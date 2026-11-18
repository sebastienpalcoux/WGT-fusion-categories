#!/usr/bin/env python3
"""Independent minor and Pluecker replay of the rank16 symplectic catalogue."""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,math,time
HERE=Path(__file__).resolve().parent

def determinant(a,b,c):
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))
def length(v):return sum(t*t for t in v)
def wedge2(a,b):return tuple(a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3))
def raw(m):
    g=[0]*len(m);k=0
    for inds in combinations(range(len(m)),3):
        z=determinant(*(m[i] for i in inds))**2;k+=z
        for i in inds:g[i]+=z
    return k,tuple(sorted(v for v in g if v))

def replay_unit(targets):
    limit=max(targets);r=math.isqrt(limit-1)
    vectors=[]
    for v in product(range(-r,r+1),repeat=3):
        n=length(v)
        if 0<n<limit and next(x for x in v if x)>0:vectors.append((n,v))
    vectors.sort();out={k:set()for k in targets};count=0;frames=0
    def visit(cols,start,norm,energy):
        nonlocal count,frames
        if norm in out:
            k,p=raw(((1,0,0),(0,1,0),(0,0,1))+cols)
            assert k==norm
            out[k].add(p);frames+=1
        if len(cols)==4:return
        for i in range(start,len(vectors)):
            n,v=vectors[i]
            if energy+n>=limit:break
            count+=1
            # Cauchy--Binet, computed by explicit minors, not det Gram.
            new=norm+n+sum(length(wedge2(v,u))for u in cols)
            new+=sum(determinant(v,u,w)**2 for u,w in combinations(cols,2))
            if new<=limit:visit(cols+(v,),i,new,energy+n)
    visit((),0,1,0)
    return out,{'minor_chart_extensions':count,'minor_frames':frames}

def replay_no_unit(targets):
    triples=list(combinations(range(7),3));where={t:i for i,t in enumerate(triples)}
    out={k:set()for k in targets};inputs=accepted=0;partition_counts={}
    def partitions(n,minimum=2):
        if not n:yield ();return
        for a in range(minimum,math.isqrt(n)+1):
            for tail in partitions(n-a*a,a):yield(a,)+tail
    for k in targets:
        for parts in partitions(k):
            # Put one maximal coefficient on 012; label and global sign
            # changes allow this without loss. All remaining magnitudes in
            # the present range coincide; this is checked rather than assumed.
            a=parts[-1];rest=parts[:-1]
            assert not rest or len(set(rest))==1,(k,parts)
            before=accepted
            for supports in combinations(range(1,35),len(rest)):
                for signs in product((-1,1),repeat=len(rest)):
                    inputs+=1;d=[0]*35;d[0]=a
                    for idx,mag,sgn in zip(supports,rest,signs):d[idx]=mag*sgn
                    m=[(a,0,0),(0,a,0),(0,0,a)]
                    for j in range(3,7):m.append((d[where[(1,2,j)]],-d[where[(0,2,j)]],d[where[(0,1,j)]]))
                    # A form is decomposable iff its affine chart reproduces
                    # every coordinate. The scaled chart has denominator a^2.
                    if any(determinant(*(m[i] for i in t))!=a*a*d[h] for h,t in enumerate(triples)):continue
                    accepted+=1;g=[0]*7
                    for h,t in enumerate(triples):
                        for i in t:g[i]+=d[h]*d[h]
                    out[k].add(tuple(sorted(v for v in g if v)))
            partition_counts[f'{k}:'+','.join(map(str,parts))]=accepted-before
    return out,{'coefficient_arrays':inputs,'decomposable_arrays':accepted,'partitions':partition_counts}

def main():
    started=time.monotonic();p=HERE/'symplectic_profiles.json';data=json.loads(p.read_text());targets=data['targets']
    unit,a=replay_unit(targets);other,b=replay_no_unit(targets)
    for k in targets:
        got=unit[k]|other[k];expected={tuple(x)for x in data['profiles'][str(k)]}
        assert got==expected,(k,got-expected,expected-got)
        assert all(sum(x)==3*k and len(x)<=7 and min(x)>0 for x in got)
    for w in data['unit_witnesses']:
        k,pf=raw(((1,0,0),(0,1,0),(0,0,1))+tuple(map(tuple,w['columns'])))
        assert(k,pf)==(w['k'],tuple(w['profile']))
    report={'status':'PASS','scope':'Complete independent raw-minor replay for every requested norm, and exhaustive no-unit integer coefficient arrays on seven slots.',
            'counts':data['counts'],'minor_audit':a,'no_unit_audit':b,'seconds':time.monotonic()-started,
            'catalogue_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'symplectic_profile_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
