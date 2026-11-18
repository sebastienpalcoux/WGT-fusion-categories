#!/usr/bin/env python3
"""Independent exact induced-representation check of Lemma 3.28.

All arithmetic is in Q[z]/(z^4+1), with Python integer coefficients.
No floating point, eigensolver, or author-supplied certificate is used.
"""
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path
import argparse
if not __debug__:
    raise RuntimeError("Assertions must be enabled")
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--output",type=Path,required=True)
_args=_parser.parse_args()
_args.output.parent.mkdir(parents=True,exist_ok=True)
import json
import numpy as np

N=12
def normalize(a,d=1):
    h=reduce(gcd,(int(x) for x in a.flat),abs(d))
    return a//h,d//h
def zero(): return np.zeros((4,N,N),dtype=object),1
def eye():
    a,d=zero(); a[0]=np.eye(N,dtype=object); return a,d
def add(A,B):
    a,d=A; b,e=B
    return normalize(e*a+d*b,d*e)
def scale(A,n,d=1): return normalize(A[0]*n,A[1]*d)
def mul(A,B):
    a,d=A; b,e=B; c=np.zeros_like(a)
    for i in range(4):
        for j in range(4):
            c[(i+j)%4]+=(-1 if i+j>=4 else 1)*(a[i]@b[j])
    return normalize(c,d*e)
def phase(A,k):
    a,d=A; b=np.zeros_like(a); k%=8
    for i in range(4):
        q,r=divmod(i+k,4); b[r]+=(-1 if q%2 else 1)*a[i]
    return b,d
def iszero(A): return not any(A[0].flat)
def adjoint(A):
    a,d=A; b=np.zeros_like(a)
    for i in range(4):
        q,r=divmod((-i)%8,4)
        b[r]+=(-1 if q%2 else 1)*a[i].T
    return b,d
def same(A,B): return iszero(add(A,scale(B,-1)))
def gm(a,b): return (a@b)%8
def gi(a): return np.array([[a[1,1],-a[0,1]],[-a[1,0],a[0,0]]],dtype=int)%8
cosets=[np.array([[x,-1],[1,0]])%8 for x in range(8)]
cosets += [np.array([[1,0],[2*y,1]]) for y in range(4)]
def induced(g,k):
    a,d=zero()
    for j,rj in enumerate(cosets):
        choices=[]
        for i,ri in enumerate(cosets):
            b=gm(gi(ri),gm(g,rj))
            if b[1,0]==0:
                choices.append((i,(-k*int(b[0,1])*pow(int(b[0,0]),-1,8))%8))
        assert len(choices)==1
        i,t=choices[0]; a[t%4,i,j]=(-1 if t>=4 else 1)
    return a,d

I=eye(); results=[]
for k in range(8):
    U=induced(np.array([[1,1],[0,1]]),k)
    W=induced(np.array([[0,-1],[1,0]])%8,k)
    H=zero()
    for a in [1,3,5,7]: H=add(H,induced(np.array([[a,0],[0,a]]),k))
    powers=[I]
    for j in range(1,8): powers.append(mul(powers[-1],U))
    assert same(mul(powers[-1],U),I)
    assert same(mul(W,W),I) and same(adjoint(W),W)
    Q0=zero()
    for j in range(8): Q0=add(Q0,phase(powers[j],k*j))
    P=scale(mul(H,Q0),1,32)
    assert same(mul(P,P),P) and same(adjoint(P),P)
    A=mul(mul(P,W),P)
    assert same(adjoint(A),A)
    # In the gamma=-1 case the table reads t(t+1)(2t-1).
    if k==0: factors=[A,add(scale(A,2),I),add(A,scale(I,-1))]
    elif k==4: factors=[A,add(A,I),add(scale(A,2),scale(I,-1))]
    else: factors=[A,add(scale(A,2),I),add(scale(A,2),scale(I,-1))]
    assert iszero(mul(mul(factors[0],factors[1]),factors[2]))
    for t in range(1,8):
        Q=zero()
        for j in range(8): Q=add(Q,phase(powers[j],(k-t)*j))
        Q=scale(Q,1,8)
        assert same(mul(Q,Q),Q) and same(adjoint(Q),Q)
        C=mul(mul(mul(mul(P,W),Q),W),P)
        assert same(adjoint(C),C)
        if k in [0,4] and t==4:
            assert same(mul(C,C),scale(C,3,4))
            assert same(mul(A,C),scale(C,(-1 if k==0 else 1),2))
        else:
            assert iszero(mul(mul(C,add(scale(C,4),scale(I,-1))),add(scale(C,2),scale(I,-1))))
    results.append({'gamma_exponent':k,'all_seven_twists_pass':True,
                    'B_projector_rank':int(sum(P[0][0,i,i] for i in range(N))//P[1])})
out={'coefficient_ring':'Q[z]/(z^4+1)','gamma_values':8,'nonidentity_theta_values':7,
     'annihilating_and_equality_identities':'all passed','results':results}
_args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
