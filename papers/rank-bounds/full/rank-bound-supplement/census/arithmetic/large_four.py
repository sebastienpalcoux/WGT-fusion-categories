"""Upper bounds on reciprocal-square mass by prime omissions through rank 40."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import ceil, prod, lcm

from row_quota import row_quota

PURE, MIXED, PRIME = range(3)
PRIMES = (3,5,7,11,13,17,19)
EVEN = tuple(mask for mask in range(8) if bin(mask).count("1")<=2)
ODD = (0,1,2,4)


@lru_cache(None)
def mask_lists(size, masks):
    return [(word,tuple(sum(bool(m&bit) for m in word) for bit in (1,2,4)))
            for word in combinations_with_replacement(masks,size)]


@lru_cache(None)
def quotas(p,rows):
    return row_quota(p,rows,False), row_quota(p,rows,True)


def profiles(rank, rows, five, support, modes):
    """Yield possible omission multisets and their upper mass."""
    if five and 5 not in support:
        return
    classes = [i for i,mode in enumerate(modes) if mode==PRIME]
    seven = 7 in (support[i] for i in classes)
    if (five and rows==13) or (seven and rows==13):
        return  # Positive mass, or T7=0 and odd T4 with too little residual.
    powers = [p**3 if p>=5 and mode==PURE else p*p if mode==PURE else p for p,mode in zip(support,modes)]
    min_b = 8 if rows==12 and five else 6
    active = []
    for p,mode in zip(support,modes):
        count,b = ((p,p*p),(4*(p-1),2*p),(p,p))[mode]
        if mode==MIXED and b<min_b:
            b *= (min_b+b-1)//b
        active.append((count,b))
    room = rank-rows-1-sum(n for n,b in active)
    base = F(rows,16)+F(1,(8*prod(powers))**2)+sum(F(n,b*b) for n,b in active)
    if (room<3 or base+F(room,min_b*min_b)<1 or
        F(rows,16)+sum(F(n,b*b) for (n,b),mode in zip(active,modes) if mode==PRIME)>=1):
        return
    weighted,plain = zip(*(quotas(p,rows) for p in support))
    # Omitted-prime information outside the hit set only enters through
    # the gcd of the proper free degrees at a prime denominator.
    class_mask = sum(1<<i for i in classes)
    other_masks = (0,)+tuple(mask for mask in EVEN if mask&class_mask)
    j5 = support.index(5) if five else None
    bit5 = (1<<j5) if five else 0
    j7 = support.index(7) if seven else None
    bit7 = (1<<j7) if seven else 0
    weights = {}
    for t in (1,4,8):
        for mask in EVEN:
            k = t*prod(pe for i,pe in enumerate(powers) if mask & (1<<i))
            b = ((min_b-1)//k+1)*k
            while (any(b%p==0 and b%(p*p)!=0 for p,mode in zip(support,modes) if p>=7 and mode==PURE)
                   or (b in support and modes[support.index(b)]!=PRIME)
                   or (mask&bit7 and rows==12 and (b<=7 or b//7%7 not in (1,6)))
                   or (mask&bit5 and (b<=5 or b//5%5 not in (1,4)
                       or (rows==12 and (b<80 or b%20))))):
                b += k
            weights[t,mask] = F(1,b*b)
    # Integer weights avoid repeated reduction in the omission sums.
    scale=lcm(base.denominator,min_b**2,64,81,*(value.denominator for value in weights.values()))
    base=int(base*scale)
    weights={key:int(value*scale) for key,value in weights.items()}
    ordinary=scale//(min_b*min_b)
    five_hit_min=15//(16-rows) if five else 0
    for o in range(3,room+1,4):
        for odd,oc in mask_lists(o,ODD):
            if five and oc[j5]!=1:
                continue  # Three odd proper labels leave no room to omit both other odd primes.
            ov = sum(weights[8,m] for m in odd)
            if base+ov+(room-o)*ordinary<scale:
                continue
            for e in range(room-o+1):
                for even,ec in mask_lists(e,EVEN):
                    if any(4*a+b<w or a+b<u for a,b,w,u in zip(oc,ec,weighted,plain)):
                        continue
                    if five and oc[j5]+ec[j5]<five_hit_min:
                        continue
                    if five and rows==12 and oc[j5]+ec[j5]!=4:
                        continue
                    total = base+ov+sum(weights[4,m] for m in even)
                    if total+(room-o-e)*ordinary<scale:
                        continue
                    for other,nc in mask_lists(room-o-e,other_masks):
                        if five and oc[j5]+ec[j5]+nc[j5]!=4:
                            continue
                        if seven and rows==12 and oc[j7]+ec[j7]+nc[j7]!=6:
                            continue
                        if any((1<<i) not in odd or any((1<<i | 1<<j) not in even+other
                               for j in range(3) if j!=i) for i in classes):
                            continue
                        mass = total+sum(weights[1,m] for m in other)
                        if five and rows==12:
                            eights=sum(weights[t,m]==scale//64 for word,t in
                                       ((odd,8),(even,4),(other,1)) for m in word)
                            mass-=max(0,eights-2)*(scale//64-scale//81)
                        if mass>=scale:
                            yield odd,even,other,F(mass,scale)


def preliminary(rank):
    for rows in range(9,14):
        for five in (False,True):
            for ps in combinations(PRIMES,3):
                if sum(p-1 for p in ps)>rank-1-rows-3:
                    continue
                choices=[((PRIME,) if five else (PURE,MIXED)) if p==5
                         else (PURE,MIXED) if p==3 else (PURE,PRIME) for p in ps]
                for modes in product(*choices):
                    yield rows,five,ps,modes

def main():
    from phase_deficits import pure_check

    mixed=set()
    for rows,five,support,modes in preliminary(40):
        candidates=profiles(40,rows,five,support,modes)
        first=next(candidates,None)
        if first is None:
            continue
        if MIXED not in modes:
            neutral=6 if five and modes.count(PRIME)>1 else 4 if five else 3 if PRIME in modes else 2
            roots=pure_check(40,rows,support,modes,neutral)
            exception=rows==9 and five and support==(3,5,7) and modes[2]==PRIME
            assert set(roots)<=({15482880} if exception else set())
            continue
        assert modes[0]==MIXED and MIXED not in modes[1:]
        assert rows<=12 and support==(3,5,7)
        if rows<=11 and not five:
            assert modes[2]==PRIME
        if rows==12 or rows==11 and five:
            assert modes[2]!=PRIME
        mixed.add((rows,five,support,modes))
    assert len(mixed)==10
    print("Mixed prime and four-row profiles through rank 40:",len(mixed))

if __name__=='__main__':
    main()
