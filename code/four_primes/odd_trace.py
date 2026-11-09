"""Final singleton trace after reserving one full orbit at an odd endpoint.

Test candidate ordinary roots from odd_periodic.py. At one endpoint a common
degree ramifies, and uses the only proper ramified orbit. Measure cost per
dual pair and compute the exact remaining norm and degree trace modulo successively larger powers of Q.
"""
from argparse import ArgumentParser
from math import lcm
from odd_periodic import coordinate, powers, scan


def endpoint_possible(P,h,q,t,w,j,jw,sign,modulus=None):
    modulus=modulus or q*q
    budget=j*(t-j)//2
    target=(sign*j//2)%modulus
    # The exponent zero is separate, because being on an axis matters.
    q_powers=[(1,False,1)]
    q_powers.extend((r,True,cap) for r,cap in powers(q,lcm(t*w,modulus),P**h))
    types=set()
    for i in range(h+1):
        for r,has_q,cap in q_powers:
            for has_w in (False,True):
                # Every nonunit degree contains at least two distinct primes.
                if int(i>0)+int(has_q)+int(has_w)<2:
                    continue
                d=P**i*r*(w if has_w else 1)
                co=coordinate(d,t,j,sign)
                if co is None:
                    continue
                a,free=co
                common=not has_w
                # At the chosen endpoint every proper ramified label is common.
                if a%2 and not common:
                    continue
                if common:
                    other=coordinate(d,w,jw,sign)
                    if other is None:
                        continue
                    at_w,free_w=other
                    g=min(P**i,cap)
                    usable_t=free and (P-1)*t<=P*g
                    usable_w=free_w and (P-1)*w<=P*g
                    if not (usable_t or usable_w):
                        continue
                    # The mixed Gram parity permits selecting the ramified orbit
                    # with an odd centered degree at the opposite endpoint.
                    if a%2 and (not usable_w or (at_w//jw)%2==0):
                        continue
                kind=1 if not has_q else 2 if not i else 4 if common and a%2 else 0
                if a*a<=budget:
                    types.add((a*a,d*a%modulus,kind))
    # Only one ramified orbit can occur. Its j labels form j/2 dual pairs.
    states={(0,0):1}
    mask=(1<<modulus)-1
    for cost,step,kind in sorted(types,reverse=True):
        updated=dict(states)
        for (used,seen),residues in states.items():
            upper=(budget-used)//cost
            if kind==4:
                counts=(j//2,) if not seen&4 and upper>=j//2 else ()
            else:
                counts=range(1,upper+1)
            for m in counts:
                shift=m*step%modulus
                reached=((residues<<shift)&mask)|(residues>>(modulus-shift))
                key=(used+m*cost,seen|kind)
                updated[key]=updated.get(key,0)|reached
        states=updated
    return bool(states.get((budget,7),0)&(1<<target))


def audit(P,h):
    rows=[]
    scan(P,h,on_survivor=rows.append)
    stages=[rows]
    for power in (2,3,4):
        survivors=[]
        for q,u,v,ju,jv,sign in stages[-1]:
            modulus=q**power
            if (endpoint_possible(P,h,q,u,v,ju,jv,sign,modulus) or
                    endpoint_possible(P,h,q,v,u,jv,ju,sign,modulus)):
                survivors.append((q,u,v,ju,jv,sign))
        stages.append(survivors)
    return stages


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("P",type=int)
    parser.add_argument("heights",type=int,nargs="+")
    args=parser.parse_args()
    for h in args.heights:
        rows,*levels=audit(args.P,h)
        print(f"P={args.P} h={h} projection rows {rows}",flush=True)
        for exponent,survivors in zip((2,3,4),levels):
            print(f"  Q^{exponent} trace survivors {survivors}",flush=True)
