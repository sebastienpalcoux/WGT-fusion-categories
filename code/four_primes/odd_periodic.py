"""Necessary simultaneous two-axis tests at an odd least prime P, all Q-heights.

Only the fixed P-axis degrees P^i w enter the first trace. Computing this
first removes most Q before any power signatures or shared Gram steps.
"""
from argparse import ArgumentParser
from functools import lru_cache
from math import isqrt
from independent_projection import centered, inertia, primes


@lru_cache(maxsize=10000)
def powers(q,mod,cap):
    output=set()
    r,s=q%mod,min(q,cap)
    while (r,s) not in output:
        output.add((r,s))
        r,s=(r*q)%mod,min(cap,s*q)
    return tuple(output)


def deficits(t):
    return [j for j in range(2,(t+1)//3+1,2)
            if (t-1)%j==0 and (t-1)//j%2]


def coordinate(d,t,j,sign):
    r=centered(d,t)
    a=centered(-sign*j*d,t)
    if j*abs(r)*(t-abs(r))>=t*t:
        return None
    free=a==-sign*j*r and j*(t-1)*r*r<t*t
    if a%2==0 and not free:
        return None  # Every proper row with even signed coordinate is unramified.
    return a,free


def axis(t,w,j,sign,other,degree_pairs,capacity,allow_ram=True):
    options={}
    for dt,do in degree_pairs:
        result=coordinate(dt,t,j,sign)
        if result is None:
            continue
        a,free=result
        floor=max(2,inertia(t,w),inertia(other,w))
        floor+=floor%2
        if floor*a*a>capacity:
            continue
        if a%2 and (not allow_ram or floor>j):
            continue
        key=(a*a,do*a%other,bool(a%2))
        options[key]=min(options.get(key,floor),floor)
    minima={}
    for (unit_cost,step,ramified),floor in options.items():
        updated=dict(minima)
        for r,cost in [(0,0),*minima.items()]:
            upper=(capacity-cost)//unit_cost
            counts=([j] if j<=upper else []) if ramified else range(floor,upper+1,2)
            for m in counts:
                residue=(r+m*step)%other
                price=cost+m*unit_cost
                updated[residue]=min(updated.get(residue,capacity+1),price)
        minima=updated
    return minima.get(sign*j%other,capacity+1)


def test(P,H,q,u,v,ju,jv,sign,t_ram=None):
    power=P**H
    pf=[P**i for i in range(1,H+1)]
    ends=((u,v,ju),(v,u,jv))
    budgets=[]
    for t,w,j in ends:
        budget=j*(t+1-j)-j  # subtract the guaranteed inactive full orbit
        p_options={(v0*w%t,v0*w%q) for v0 in pf}
        budget-=axis(t,w,j,sign,q,p_options,budget,t!=t_ram)
        if budget<=0:
            return 0
        budgets.append(budget)
    # Vacuum equations uniquely determine Q^k modulo u*v. Keep all k.
    x=-sign*pow(power*v*ju,-1,u)%u
    y=-sign*pow(power*u*jv,-1,v)%v
    m=x+u*((y-x)*pow(u,-1,v)%v)
    if not any(r==m for r,_ in powers(q,u*v,power)):
        return 0
    for i,(t,w,j) in enumerate(ends):
        q_options={(r*w%t,r*w%P) for r,_ in powers(q,P*t,power)}
        budgets[i]-=axis(t,w,j,sign,P,q_options,budgets[i],t!=t_ram)
        if budgets[i]<0:
            return 1
    cu,cv=budgets
    choices={}
    for r,cap in powers(q,u*v,power):
        for p_part in pf:
            g=min(p_part,cap)
            du=coordinate(r*p_part,u,ju,sign)
            dv=coordinate(r*p_part,v,jv,sign)
            if du is None or dv is None:
                continue
            a,fu=du
            b,fv=dv
            free_ends=int(fu and (P-1)*u<=P*g)+int(fv and (P-1)*v<=P*g)
            if not free_ends:
                continue
            floor=max(2,3-free_ends,inertia(u,g),inertia(v,g))
            floor+=floor%2
            orbit=ju if a%2 else jv if b%2 else 0
            if orbit and floor>orbit:
                continue
            if cu>=floor*a*a and cv>=floor*b*b:
                choices[a,b]=min(choices.get((a,b),floor),floor)
    target=-ju*jv
    states={(0,0,False):0}
    for (a,b),floor in choices.items():
        updated=dict(states)
        for (cost,mixed,hit),price in states.items():
            orbit=ju if a%2 else jv if b%2 else 0
            upper=min((cu-cost)//(a*a),(cv-price)//(b*b))
            counts=([orbit] if orbit<=upper else []) if orbit else range(floor,upper+1,2)
            for m in counts:
                total=mixed+m*a*b
                use=bool((t_ram==u and a%2) or (t_ram==v and b%2))
                if hit and use:
                    continue  # only one proper ramified orbit at the chosen endpoint
                seen=hit or use
                if total==target and (t_ram is None or seen):
                    return 3
                key=(cost+m*a*a,total,seen)
                updated[key]=min(updated.get(key,cv+1),price+m*b*b)
        states=updated
    return 2


def scan(P,H,on_survivor=None):
    if P==2 or P<3 or any(P%d==0 for d in range(2,isqrt(P)+1)) or H<2:
        raise ValueError("P must be an odd prime and H an integer at least two")
    stats=[0,0,0,0]
    remaining=set()
    for u in primes(P**(H+1)//(P-1)):
        if u<=P or u in (5,7,13):
            continue
        for v in primes(P*u//(P-1)):
            if v<=u or v in (5,7,13):
                continue
            for ju in deficits(u):
                bound=v*max(2,(u+1)//ju-1)+1
                for q in primes(bound):
                    if q<=P or q==u or q==v:
                        continue
                    for jv in deficits(v):
                        if ju*v+jv*u>=u*v:
                            continue
                        for sign in (-1,1):
                            level=test(P,H,q,u,v,ju,jv,sign)
                            if level==3:
                                if all(test(P,H,q,u,v,ju,jv,sign,t)<3 for t in (u,v)):
                                    level=2
                            for i in range(level+1):
                                stats[i]+=1
                            if level==3:
                                remaining.add(q)
                                if on_survivor is not None:
                                    on_survivor((q,u,v,ju,jv,sign))
    return stats,sorted(remaining)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("P",type=int)
    parser.add_argument("heights",nargs="+",type=int)
    args=parser.parse_args()
    for h in args.heights:
        counts,qs=scan(args.P,h)
        print(f"P={args.P} h={h} classes {counts}; Q projection {qs}",flush=True)
