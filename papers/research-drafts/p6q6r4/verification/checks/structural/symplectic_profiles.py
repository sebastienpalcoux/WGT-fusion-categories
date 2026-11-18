#!/usr/bin/env python3
"""Exact decomposable integral 3-form slice profiles on at most seven slots.

Unit-coordinate charts enumerate integral [I_3|A] by positive, sorted columns.
The norm is det(I+AA^t); an independent verifier uses all raw 3x3 minors.
No-unit charts are the finite square partitions explained in the companion proof.
All arithmetic is integral.  No numerical or algebraic-number dependency is used.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGETS = tuple(range(1,16)) + (18,25,27)

def dot(x,y):
    return sum(a*b for a,b in zip(x,y))

def det3(x,y,z):
    return (x[0]*(y[1]*z[2]-y[2]*z[1])
            -x[1]*(y[0]*z[2]-y[2]*z[0])
            +x[2]*(y[0]*z[1]-y[1]*z[0]))

def gram_data(g):
    a,b,c,d,e,f=g
    adj=(d*f-e*e,c*e-b*f,b*e-c*d,a*f-c*c,b*c-a*e,a*d-b*b)
    return a*adj[0]+b*adj[1]+c*adj[2], adj

def qform(v,g):
    x,y,z=v; a,b,c,d,e,f=g
    return a*x*x+d*y*y+f*z*z+2*(b*x*y+c*x*z+e*y*z)

def gram_add(g,v):
    a,b,c,d,e,f=g; x,y,z=v
    return (a+x*x,b+x*y,c+x*z,d+y*y,e+y*z,f+z*z)

def chart_profile(g,cols):
    k,adj=gram_data(g)
    return k, tuple(sorted((adj[0],adj[3],adj[5])+tuple(qform(v,adj) for v in cols)))

def raw_profile(cols,scale=1):
    m=((1,0,0),(0,1,0),(0,0,1))+tuple(cols)
    slices=[0]*len(m); norm=0
    for ijk in itertools.combinations(range(len(m)),3):
        v=scale*det3(*(m[i] for i in ijk)); z=v*v
        norm+=z
        for i in ijk:slices[i]+=z
    return norm,tuple(sorted(x for x in slices if x))

def positive(v):
    return next(x for x in v if x)!=0 and next(x for x in v if x)>0

def unit_charts(targets=TARGETS):
    maximum=max(targets); radius=math.isqrt(maximum-1)
    vectors=[]
    for v in itertools.product(range(-radius,radius+1),repeat=3):
        n=dot(v,v)
        if 0<n<maximum and positive(v):vectors.append((n,v))
    vectors.sort()
    profiles={k:set() for k in targets}; witnesses={}; counters={"vectors":len(vectors),"candidate_extensions":0,"accepted_frames":0}
    def visit(g,cols,start,entry_norm):
        k,p=chart_profile(g,cols)
        if k in profiles:
            profiles[k].add(p)
            witnesses.setdefault((k,p),cols)
            counters["accepted_frames"]+=1
        if len(cols)==4:return
        for i in range(start,len(vectors)):
            n,v=vectors[i]
            if 1+entry_norm+n>maximum:break
            counters["candidate_extensions"]+=1
            ng=gram_add(g,v); nk,_=gram_data(ng)
            if nk>maximum:continue
            visit(ng,cols+(v,),i,entry_norm+n)
    visit((1,0,0,1,0,1),(),0,0)
    return profiles,witnesses,counters

def no_unit_profiles():
    # Every nonzero coefficient has absolute value at least 2.  For the
    # requested norms the following are ALL partitions into such squares.
    # Up to three supported coordinates: decomposability forces the listed
    # common-face configurations. The exceptional 25=9+4+4+4+4 is audited below.
    p={4:{(4,4,4)},8:{(4,4,8,8)},9:{(9,9,9)},
       12:{(8,8,8,12),(4,4,4,12,12)},13:{(4,9,13,13)},
       18:{(9,9,18,18)},25:{(25,25,25),(9,16,25,25)},
       27:{(18,18,18,27),(9,9,9,27,27)}}
    exceptional=set(); arrays=0; integral=0; shapes={}
    # Pivot coefficient 3.  The other four coefficients are ±2. The
    # one-exterior entries form a 3x4 numerator matrix B; all 2x2 minors
    # vanish since 4,8 are not divisible by 3. Hence B has rank one.
    for positions in itertools.combinations(range(12),4):
        for signs in itertools.product((-2,2),repeat=4):
            arrays+=1
            entries=[0]*12
            for i,v in zip(positions,signs):entries[i]=v
            rows=[entries[4*i:4*i+4] for i in range(3)]
            if any(rows[i][a]*rows[j][b]!=rows[i][b]*rows[j][a]
                   for i,j in itertools.combinations(range(3),2)
                   for a,b in itertools.combinations(range(4),2)):continue
            integral+=1
            support_rows=sum(any(row) for row in rows)
            support_cols=sum(any(rows[i][j] for i in range(3)) for j in range(4))
            shape=f"{support_rows}x{support_cols}"
            shapes[shape]=shapes.get(shape,0)+1
            # The coordinate 012 has square 9; each nonzero B_ij has
            # square 4 and replaces pivot index i with outside index 3+j.
            slices=[9,9,9,0,0,0,0]
            for i in range(3):
                for j in range(4):
                    if rows[i][j]:
                        for a in range(3):
                            if a!=i:slices[a]+=4
                        slices[3+j]+=4
            exceptional.add(tuple(sorted(x for x in slices if x)))
    p[25].update(exceptional)
    return p,{"exceptional_arrays":arrays,"exceptional_integral":integral,"exceptional_shapes":shapes}

def build(output):
    p,w,c=unit_charts(); extra,e=no_unit_profiles()
    for k,values in extra.items():p[k].update(values)
    data={"schema":"r36-symplectic-profiles-v1","scope":"All integral decomposable alternating 3-forms of these norms on at most seven active slots, up to slice-profile equality.",
          "targets":list(TARGETS),"profiles":{str(k):[list(x) for x in sorted(p[k])] for k in TARGETS},
          "counts":{str(k):len(p[k]) for k in TARGETS},"unit_enumeration":c,"no_unit_enumeration":e,
          "unit_witnesses":[{"k":k,"profile":list(pr),"columns":[list(v) for v in cols]} for (k,pr),cols in sorted(w.items())],
          "producer_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps({"status":"GENERATED","counts":data["counts"],"unit_enumeration":c,"no_unit_enumeration":e,"path":str(output)},sort_keys=True))

def profiles(k,b=7,path=None):
    data=json.loads(Path(path or HERE/'symplectic_profiles.json').read_text())
    return tuple(tuple([0]*(b-len(p))+p) for p in data['profiles'][str(k)] if len(p)<=b)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'symplectic_profiles.json')
    a=ap.parse_args();build(a.output)
