#!/usr/bin/env python3
"""Independent exact checks for p6q6r4_manuscript.pdf.
This checks arithmetic, not existence or categorifiability, and does not replay
unprovided certificate files. Python 3, NumPy, and SymPy are required.
"""
from pathlib import Path
import argparse
if not __debug__:
    raise RuntimeError("Assertions must be enabled")
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--output",type=Path,required=True)
_args=_parser.parse_args()
_args.output.parent.mkdir(parents=True,exist_ok=True)
from itertools import combinations, combinations_with_replacement, product
from math import isqrt
import json
import numpy as np
import sympy as sp

out={}
limit=5**6*7**6*23**4
out['dimension_arithmetic']={'even':2**6*3**6*5**4,'even_square':5400**2,'odd':limit,'odd_square':22680875**2,'next_least_prime_floor':7**6*11**6*13**4,'5399_factorization':{str(k):int(v) for k,v in sp.factorint(5399).items()}}
triples=[]
primes=list(sp.primerange(3,100))
for p,q,r in combinations(primes,3):
    dim=p**6*q**6*r**4
    if dim>limit: continue
    m=min(r,p*p)//q
    E=m*(m+1)//2
    H=(r-1)//2
    possible=[]
    for h in range(1,E+1):
        k=1
        while k*(k+1)<h: k+=1
        if (H-1)*k<=h*(h-1): possible.append([h,k,h*(h-1)-(H-1)*k])
    triples.append({'p':p,'q':q,'r':r,'dimension':dim,'m':m,'possible_h_k_slack':possible})
# Completeness: p>=7 is already too large; p>=3,q>=5 give r<83.
assert 7**6*11**6*13**4>limit
assert 3**6*5**6*83**4>limit
out['odd_frontier']={'ordered_triples':len(triples),'scalar_excluded':sum(not t['possible_h_k_slack'] for t in triples),'survivors':[t for t in triples if t['possible_h_k_slack']],'all_triples':triples}
A={1,3,4,5,9}
out['difference_set_Z11']={'set':sorted(A),'overlaps':[len(A & {(x+t)%11 for x in A}) for t in range(11)],'order_of_3':int(sp.n_order(3,11)),'minus_one_in_3_subgroup':10 in {pow(3,j,11) for j in range(5)}}
palette=[int(a) for a in sp.divisors(5400) if 10<=a<=900 and len(sp.factorint(a))>=2]
out['canonical_palette']=palette

def gf_counts(maxlength,representation_partitions=False):
    # Bound every intermediate coefficient by a one-variable coefficient,
    # computed with unbounded Python integers, before using int64 arrays.
    bound=[1]+[0]*maxlength
    for _ in palette:
        for m in range(1,(maxlength if representation_partitions else 1)+1):
            for n in range(m,maxlength+1): bound[n]+=bound[n-m]
    assert max(bound)<np.iinfo(np.int64).max
    dp=np.zeros((maxlength+1,5400),dtype=np.int64); dp[0,0]=1
    for a in palette:
        for m in range(1,(maxlength if representation_partitions else 1)+1):
            for length in range(m,maxlength+1):
                cost=a*m
                if cost<5400: dp[length,cost:]+=dp[length-m,:-cost]
    return dp
profiles=gf_counts(13)
assignments=gf_counts(12,True)
out['weighted_profiles']={str(n+1):int(profiles[n,5399]) for n in range(9,14)}
out['representation_multiplicity_assignments']={str(n+1):int(assignments[n,5399]) for n in range(9,13)}
out['one_matrix_placements']={str(r):int(sum(profiles[r-5,5399-2*a] for a in palette if 2*a<=5399)) for r in [13,14,15,16]}
out['two_matrix_placements']={str(r):int(sum(profiles[r-9,5399-2*a-2*b] for a,b in combinations_with_replacement(palette,2) if 2*a+2*b<=5399)) for r in [15,16]}
# Exactly one canonical degree ten: remove its factor, then insert one occurrence.
pal_no10=[a for a in palette if a!=10]
dp=np.zeros((12,5400),dtype=np.int64); dp[0,0]=1
for a in pal_no10:
    for n in range(1,12): dp[n,a:]+=dp[n-1,:-a]
out['unique_degree10_profiles']={str(r):int(dp[r-2,5389]) for r in [12,13]}
# Reconstruct every solution in Theorem 18.4 without its producer or row kernel.
rows=[]; triple_count=0
for x in [-1,0,1]:
    R=5400-3600*x*x
    for u in range(-isqrt(R//400),isqrt(R//400)+1):
        for v in range(-isqrt(R//90),isqrt(R//90)+1):
            for s in range(-isqrt(R//90),isqrt(R//90)+1):
                if s%2: continue
                triple_count+=1
                for d in range(3,52):
                    num=-(d+200*u+45*v+90*s)
                    if num%12: continue
                    w=num//12
                    if 2*d*d+400*u*u+90*v*v+12*w*w+90*s*s==R:
                        rows.append((d,s,x,-x,u,v,w))
rows.sort()
weights=[2,90,1800,1800,400,90,12]
edges=[]
for i,j in combinations(range(len(rows)),2):
    if sum(a*b*c for a,b,c in zip(rows[i],rows[j],weights))==0: edges.append([i+1,j+1])
edge_set={tuple(e) for e in edges}
triangles=[list(c) for c in combinations(range(1,len(rows)+1),3) if all(tuple(e) in edge_set for e in combinations(c,2))]
out['theorem18_4']={'bounded_u_v_s_triples':triple_count,'rows':rows,'edges':edges,'triangles':triangles,'all_row_norms_5400':all(sum(w*x*x for w,x in zip(weights,row))==5400 for row in rows)}
# Exact cubic arithmetic: theta satisfies theta^3-3 theta+1=0.
x=sp.Symbol('x'); poly=sp.Poly(x**3-3*x+1,x)
traces=[3,0,6]
for k in range(3,5): traces.append(3*traces[k-2]-traces[k-3])
gram=sp.Matrix(3,3,lambda i,j:traces[i+j])
a,b,c=sp.symbols('a b c')
assert sp.expand((sp.Matrix([a,b,c]).T*gram*sp.Matrix([a,b,c]))[0]-3*(a+2*c)**2-6*(b-c/sp.Integer(2))**2-sp.Rational(9,2)*c*c)==0
out['cubic_gram']={'gram':gram.tolist(),'determinant':int(gram.det()),'small_norm_elements':[(A,B,C) for A,B,C in product(range(-5,6),repeat=3) if (sp.Matrix([A,B,C]).T*gram*sp.Matrix([A,B,C]))[0]<=4]}
# Small terminal sum-of-squares contradictions printed in the manuscript.
terminal=[]
for degrees,count,target in [([8,25,26,27,33,48],10,2483),([25,26,48,54],10,2483),([14,15,16,17,18,19,20,29,31,32,33,35,48],5,1815)]:
    solutions=[list(t) for t in combinations_with_replacement(degrees,count) if sum(d*d for d in t)==target]
    terminal.append({'dimensions':degrees,'count':count,'target':target,'solutions':solutions})
out['printed_terminal_square_tests']=terminal
out['length13_shapes']=[{'rank':r,'matrix2':a,'matrix3':b,'matrix4':c,'linear':13-2*a-3*b-4*c} for r in range(17,30,2) for a in range(9) for b in range(4) for c in range(2) if 2*a+6*b+12*c==r-13 and 13-2*a-3*b-4*c>0 and (a,b,c)!=(6,0,0)]
out['lemma17_2_counterexample']={'old_rows':[[1,0]],'target':[1,1],'separator':[0,-1],'old_dot_products':[0],'target_dot_product':-1,'added_row':[1,1],'new_dot_product':-1,'target_represented_after_enlargement':True}
# Independent matrix check for the degree-two trace identity.
U=sp.diag(1,-1); V=sp.Matrix([[0,1],[1,0]]); W=sp.Matrix([[0,1],[-1,0]])
H=sp.Matrix(3,3,lambda i,j:2*sp.trace([U,V,W][i]*[U,V,W][j]))
out['degree_two_volume_identity']={'trace_UVW':int(sp.trace(U*V*W)),'omega':int(4*sp.trace(U*V*W)),'minus_det_H':int(-H.det()),'omega_squared':int((4*sp.trace(U*V*W))**2)}
_args.output.write_text(json.dumps(out,indent=2,default=int)+'\n')
for key,value in out.items():
    if key=='odd_frontier': value={k:v for k,v in value.items() if k!='all_triples'}
    print(key,':',value)
