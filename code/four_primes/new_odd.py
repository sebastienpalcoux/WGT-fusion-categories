"""Simultaneous vacuum congruences in the odd P^h Q^2 u v window.

After reserving both mandatory axes, test whether one negative common
term fits both budgets. The octic table uses h=4 and P=3,5,7,11,13,17,19.

Standard library only. Inputs are a least prime P and root height h.
Every other prime position ranges over its entire proven window.
"""
import argparse
from array import array
from functools import lru_cache
from math import isqrt
from time import monotonic

@lru_cache(maxsize=100000)
def prime(n):
    if n <= 1 or n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3,isqrt(n)+1,2))



@lru_cache(maxsize=128)
def canonical_roots(p):
    """Store one of the two roots for every nonzero square modulo p."""
    roots = array("L",[0])*p
    for x in range(1, (p+1)//2):
        roots[x*x % p] = x
    return roots


def square_roots(p,residue):
    x=canonical_roots(p)[residue]
    return (x,p-x) if x else ()


@lru_cache(maxsize=None)
def deficits(p):
    result = []
    for d in range(1, isqrt(p-1)+1):
        if (p-1) % d:
            continue
        for j in {d,(p-1)//d}:
            if j%2 == 0 and ((p-1)//j)%2 == 1 and 3*j <= p+1:
                result.append(j)
    return tuple(sorted(result))


def residue_cases(P,h):
    power = P**h
    bound = ((9*P**(h-1)+2)//4 if P==3 else
             P**(h-1)*(2*(P-1)**2+P**2)//(2*(P-1)**2))
    limit = P*bound//(P-1)
    small_primes = [p for p in range(P+1,limit+1)
                    if p not in (5,7,13) and prime(p) and deficits(p)]
    for u_index,u in enumerate(small_primes):
        if u > bound:
            break
        for v_index in range(u_index+1,len(small_primes)):
            v = small_primes[v_index]
            if v*(P-1)>P*u:
                break
            inv_u = pow(u,-1,v)
            for ju in deficits(u):
                quota = max(2, (u+1)//ju-1)
                if v > power*quota+1:
                    continue
                q_limit = v*quota+1
                mod_u = pow(power*v*ju,-1,u)
                for sign in (-1,1):
                    left_roots = square_roots(u,-sign*mod_u % u)
                    if not left_roots:
                        continue
                    for jv in deficits(v):
                        mod_v = (-sign*pow(power*u*jv,-1,v))%v
                        for b in square_roots(v,mod_v):
                            for a in left_roots:
                                # Chinese remainder theorem, 0 <= q0 < uv.
                                q0 = a+u*((b-a)*inv_u % v)
                                for q in range(q0,q_limit+1,u*v):
                                    if q<=P or not prime(q):
                                        continue
                                    if u*(P-1) > P*min(power,q*q):
                                        continue
                                    core=power*q*q
                                    n=core*u*v
                                    if 0 < n-sign-core*(v*ju+u*jv) < n:
                                        yield q,u,v,ju,jv,sign


def coordinate(t,j,eps,d):
    residue=(d+t//2)%t-t//2
    signed=(-eps*j*d+t//2)%t-t//2
    size=abs(residue)
    allowed=j*size*(t-size)<t*t
    unramified=(signed==-eps*j*residue and j*(t-1)*residue*residue<t*t)
    return signed,allowed,unramified


def capacity(P,Q,h,t,w,j,eps):
    """Relax the two compulsory axis costs using odd duality multiplicity two."""
    budget=j*(t+1-j)
    for base, power in ((P,h),(Q,2)):
        costs=[]
        for exponent in range(1,power+1):
            z,legal,_=coordinate(t,j,eps,base**exponent*w)
            if legal:
                costs.append(2*z*z)
        if not costs:
            return 0
        budget -= min(costs)
    return budget


def obstruction_stage(P,h,case):
    """0 = insufficient axis budget, 1 = no negative common term, 2 = open."""
    q,u,v,ju,jv,eps=case
    cu=capacity(P,q,h,u,v,ju,eps)
    cv=capacity(P,q,h,v,u,jv,eps)
    if min(cu,cv)<=0:
        return 0
    for pi in (P**i for i in range(1,h+1)):
        for qj in (q,q*q):
            degree=pi*qj
            g=min(pi,qj)
            a,ok_u,free_u=coordinate(u,ju,eps,degree)
            b,ok_v,free_v=coordinate(v,jv,eps,degree)
            if not(ok_u and ok_v) or a*b>=0 or 2*a*a>cu or 2*b*b>cv:
                continue
            if (free_u and (P-1)*u<=P*g) or (free_v and (P-1)*v<=P*g):
                return 2
    return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prime",type=int)
    parser.add_argument("height",type=int,nargs="?",default=4)
    args=parser.parse_args()
    p,h=args.prime,args.height
    if p<3 or not prime(p) or h<3:
        parser.error("P must be an odd prime and h at least 3")
    start=monotonic()
    counts=[0,0,0]
    for case in residue_cases(p,h):
        stage=obstruction_stage(p,h,case)
        counts[stage]+=1
        if stage==2 and counts[stage]==1:
            print("First inconclusive residue class",case)
    print(f"P={p} h={h}: {sum(counts)} vacuum classes, "
          f"{counts[1]+counts[2]} with positive residual axis budgets, "
          f"{counts[2]} allowing a negative common term; {monotonic()-start:.2f}s")
    return 2 if counts[2] else 0


if __name__=="__main__":
    raise SystemExit(main())
