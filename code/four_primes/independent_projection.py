"""Independent audit of the h=3,4 dyadic projections. Standard library only.

Vacuum classes are generated from j_u,j_v and CRT. The modular multiplicative
order bounds the exponent signatures. Exact cost sets and Pareto minima
implement the two axis congruences and the shared Gram constraint.
"""
import argparse
from functools import lru_cache
from math import isqrt


@lru_cache(None)
def primes(top):
    return [n for n in range(2,top+1) if all(n%d for d in range(2,isqrt(n)+1))]


def factors(n):
    result = []
    d = 2
    while d*d<=n:
        if n%d==0:
            result.append(d)
            while n%d==0:
                n//=d
        d+=1
    if n>1:
        result.append(n)
    return result


@lru_cache(maxsize=256)
def signatures(q,u,v,height):
    # Once q^a >= 2^h its size can be replaced by 2^h. Before then the
    # exact size matters for the inertia floor. q is coprime to uv.
    modulus=u*v
    period=(u-1)*(v-1)
    for prime in factors(period):
        while period%prime==0 and pow(q,period//prime,modulus)==1:
            period//=prime
    cap = 2**height
    first_full=0
    while q**first_full<cap:
        first_full+=1
    return tuple((pow(q,i,modulus),min(q**min(i,first_full),cap))
                 for i in range(1,first_full+period))


def centered(x,t):
    return (x+t//2)%t-t//2


def coord(d,t,j,sign):
    x=centered(d,t)
    a=centered(-sign*j*d,t)
    if j*abs(x)*(t-abs(x))>=t*t:
        return None
    possible=(t-1)*j*x*x<t*t and a==-sign*j*x
    if j==1 and not possible:
        return None
    return a,possible


@lru_cache(None)
def inertia(t,g):
    if 2*g>=t:
        return 1
    return min(i for i in range(1,t) if (t-1)%i==0 and (t-1)//i<=g)


def axis_cost(t,w,j,sign,omitted,ds,capacity):
    # ds contains (degree modulo t, degree modulo omitted).
    types={}
    for dt,dr in ds:
        item=coord(dt,t,j,sign)
        if item is None:
            continue
        a,free=item
        signature=(a*a,(dr*a)%omitted)
        cost=signature[0]
        n0=max(inertia(t,w),inertia(omitted,w),1 if free else 2)
        if cost*n0<=capacity:
            types[signature]=min(types.get(signature,n0),n0)
    reachable=set()
    for (cost,step),floor in types.items():
        additions=set()
        for paid,residue in reachable|{(0,0)}:
            for n in range(floor,(capacity-paid)//cost+1):
                additions.add((paid+n*cost,(residue+n*step)%omitted))
        reachable.update(additions)
    return min((cost for cost,residue in reachable if residue==sign*j%omitted),default=capacity+1)


def has_gram(q,u,v,ju,jv,sign,height):
    powers=signatures(q,u,v,height)
    dyadic=[2**i for i in range(1,height+1)]
    capacities=[]
    for t,w,j in ((u,v,ju),(v,u,jv)):
        capacity=j*(t+1-j)-1
        twos={(d*w%t,d*w%q) for d in dyadic}
        qs={(d*w%t,1) for d,_ in powers}
        left=capacity-axis_cost(t,w,j,sign,q,twos,capacity)
        if left<0:
            return False,False
        left-=axis_cost(t,w,j,sign,2,qs,left)
        if left<0:
            return False,False
        capacities.append(left)
    cu,cv=capacities
    types={}
    for d in dyadic:
        for residue,cap in powers:
            g=min(d,cap)
            co=[coord(d*residue,t,j,sign) for t,j in ((u,ju),(v,jv))]
            if None in co:
                continue
            (a,fu),(b,fv)=co
            f=int(fu and u<=2*g)+int(fv and v<=2*g)
            if f==0:
                continue
            lower=max(3-f,inertia(u,g),inertia(v,g))
            if lower*a*a<=cu and lower*b*b<=cv:
                types[a,b]=min(types.get((a,b),lower),lower)
    # State for a fixed (u-cost, Gram sum) needs only smallest v-cost.
    states={(0,0):0}
    target=-ju*jv
    for (a,b),lower in types.items():
        delta={}
        for (paid,mixed),other in states.items():
            upper=min((cu-paid)//(a*a),(cv-other)//(b*b))
            for n in range(lower,upper+1):
                new_mixed=mixed+n*a*b
                if new_mixed==target:
                    return True,True
                key=(paid+n*a*a,new_mixed)
                delta[key]=min(delta.get(key,cv+1),other+n*b*b)
        for key,cost in delta.items():
            states[key]=min(states.get(key,cv+1),cost)
    return True,False


def audit(height):
    power=2**height
    totals=[0,0,0]
    odd_primes=set()
    for u in primes(2*power):
        if u<5:
            continue
        for ju in range(1,u):
            if (u-1)%ju and (u+1)%ju:
                continue
            quota=max(2,(u+1)//ju-1)
            for v in primes(power*quota+1):
                if v<=u:
                    continue
                upper=min(2*u,v*quota+1) if v>2*u else v*quota+1
                js=[j for j in range(1,v) if (not (v-1)%j or not (v+1)%j)
                    and ju*v+j*u<u*v]
                if not js:
                    continue
                for q in primes(upper):
                    if q in (2,u,v):
                        continue
                    powers=signatures(q,u,v,height)
                    available={residue for residue,_ in powers}
                    for jv in js:
                        for sign in (-1,1):
                            # q^k*(2^h v)*ju = -sign (mod u), with a
                            # second equation at v. CRT gives just one class.
                            a=-sign*pow(power*v*ju,-1,u)%u
                            b=-sign*pow(power*u*jv,-1,v)%v
                            residue=a+u*((b-a)*pow(u,-1,v)%v)
                            if residue not in available:
                                continue
                            totals[0]+=1
                            axis,gram=has_gram(q,u,v,ju,jv,sign,height)
                            totals[1]+=axis
                            totals[2]+=gram
                            if gram:
                                odd_primes.add(q)
    return totals,sorted(odd_primes)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("heights",type=int,nargs="+")
    args=parser.parse_args()
    for height in args.heights:
        if height<2:
            parser.error("root height must be at least two")
        totals, surviving=audit(height)
        print(f"h={height}: {totals}; heavy-prime projection {surviving}",flush=True)
