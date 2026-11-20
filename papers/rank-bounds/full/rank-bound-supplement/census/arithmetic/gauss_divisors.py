"""Complete the two possible primary dimension roots with nine denominator-four rows."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement as multisets
from math import ceil, gcd, isqrt

PRIMES = (2, 3, 5, 7)


def divisors(root):
    small = [d for d in range(1,isqrt(root)+1) if root % d == 0]
    return sorted(set(small) | {root//d for d in small})


def free_classes(root, ds, p):
    """The proper five- or seven-free degrees and their signed Weyl coefficients."""
    n = root//p
    epsilon = -1 if n % p == 1 else 1
    assert (n+epsilon) % p == 0
    options = [(d,a) for d in ds if d < n and n % d == 0 and sum(d % q == 0 for q in PRIMES)>=2
               for a in (-2,-1,1,2) if (n//d)*a % p == 1]
    singles = [(d,a) for d,a in options if abs(a)==1]
    doubles = [((d,a),) for d,a in options if abs(a)==2]
    result=set()
    for count, has_double in ([(4,False)] if p==5 else [(6,False),(2,True)]):
        lefts=defaultdict(list)
        for left in multisets(singles,count//2):
            lefts[sum(d*a for d,a in left)].append(left)
        for double in doubles if has_double else [()]:
            for right in multisets(singles,count-count//2):
                tail=right+double
                for left in lefts.get(epsilon-sum(d*a for d,a in tail),()):
                    degrees=tuple(sorted(d for d,a in left+tail))
                    if gcd(*degrees)==1 and len(set(degrees))>=3 and all(
                        gcd(d,e)>1 for d,e in multisets(degrees,2)):
                        result.add(degrees)
    return result


def dyadic_degree_denominators(value, slots, top):
    """Split binary digits to enumerate all allowed dyadic deficit multisets."""
    if value<0 or value>F(slots,32):
        return
    if not value:
        yield ()
        return
    exponent=value.denominator.bit_length()-1
    if 2**exponent!=value.denominator:
        return
    first=tuple(sorted(exponent-j for j in range(value.numerator.bit_length()) if value.numerator & (1<<j)))
    pending=[first]; seen=set()
    while pending:
        parts=pending.pop()
        if parts in seen or len(parts)>slots or max(parts)>top:
            continue
        seen.add(parts)
        if min(parts)>=5:
            yield tuple(2**((h+1)//2) for h in parts)
        if len(parts)<slots:
            for h in set(parts):
                if h<top:
                    next_parts=list(parts)
                    next_parts.remove(h)
                    pending.append(tuple(sorted(next_parts+[h+1,h+1])))


def deficit_profiles(root, budget, fours=9):
    powers=[]
    b=9
    while root % b == 0:
        powers.append(b)
        b*=3
    top=2*((root & -root).bit_length()-1)
    result=set()
    for count in range(1,budget//2+1):
        for terms in multisets(powers,count):
            ternary=sum(F(3,b*b) for b in terms)
            max_trace=int(16*(1-F(fours,16)-F(1,5)-F(1,7)+F(1,root)))
            for trace in range(-max_trace,max_trace+1):
                base=1-F(fours-trace,16)-F(1,5)-F(1,7)-ternary
                for sign in (-1,1):
                    for dyadic in dyadic_degree_denominators(base-F(sign,root),budget-2*count,top):
                        result.add((ternary,dyadic,count))
    return result


def verify_root(root, rank=42, fours=9, sixes=0):
    ds=divisors(root)
    five,seven=(free_classes(root,ds,p) for p in (5,7))
    unions={tuple(sorted((Counter(x)|Counter(y)).elements())) for x in five for y in seven
            if sorted(d for d in x if d%7)==sorted(d for d in y if d%5)}
    fixed=(1,)+fours*(root//4,)+5*(root//5,)+7*(root//7,)+sixes*(root//6,)
    profiles=deficit_profiles(root,rank-len(fixed)-min(map(len,unions)),fours)
    targets=set()
    for ternary,dyadic,count in profiles:
        for proper in unions:
            known=fixed+proper+tuple(root//b for b in dyadic)
            room=rank-len(known)
            if room<2*count:
                continue
            # The active ternary mass is >=2g_3/3. Each slot not at nine
            # contributes at most 1/27^2, including a nonternary slot.
            nines=max(0,ceil((F(2,3)*ternary-F(room,729))/(F(1,81)-F(1,729))))
            left=room-nines
            energy=root**2-sum(d*d for d in known)-nines*(root//9)**2
            if left>=0 and energy>=0 and energy%1225==0:
                targets.add((energy//1225,left))
    squares=[(d//35)**2 for d in ds if d%35==0 and root//d>4 and root//d not in (5,7)]
    for energy,left in targets:
        assert left<=5
        # Zero permits fewer labels and therefore covers every rank through rank.
        possible=[0]+[d2 for d2 in squares if d2<=energy]
        first={sum(part) for part in multisets(possible,left//2)}
        assert not any(energy-sum(part) in first for part in multisets(possible,left-left//2))
    return len(unions),len(profiles),len(targets)


if __name__=='__main__':
    for root,fours,expected in [(136080,9,(13,0,0)),
                                (15482880,9,(39,101,2698))]:
        result=verify_root(root,rank=43,fours=fours)
        assert result==expected
        print('43',root,fours,'fours; proper unions, deficit profiles, square-sum targets:',result)
