"""Exact mass and Gauss certificates for separated odd primes through rank44."""
from fractions import Fraction
from itertools import combinations, product
from math import prod, isqrt
import json

from phase_deficits import sparse_budget


# Each odd dimension prime needs at least p-1 active positions in addition to the vacuum.
PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def alternatives(p):
    """Reserved labels and denominator floor, then proper labels and their floor."""
    yield p, p*p, 0, 1
    if p in (3, 5):
        yield 4*(p - 1), 2*p, 0, 1
    if p == 5:
        yield 5, 5, 4, 20
    elif p > 5:
        yield p, p, 3, 2*p


def cases(rank):
    for primes in combinations(PRIMES, 4):
        for choices in product(*(alternatives(p) for p in primes)):
            room = rank - 1 - sum(count for count, _, _, _ in choices)
            if room < 0:
                continue
            fixed = Fraction(1, 2310**2) + sum(Fraction(n, b*b) for n, b, _, _ in choices)
            for fours in range(min(13, room) + 1):
                witnesses = [(n, b) for _, _, n, b in choices if n]
                if fours >= 9:
                    witnesses.append((3, 8))
                # The witness sets may overlap. Sort their union by decreasing denominator.
                size = max((n for n, _ in witnesses), default=0)
                floors = [max(b for n, b in witnesses if n >= j) for j in range(1, size + 1)]
                remaining = room - fours - size
                if remaining < 0:
                    continue
                mass = fixed + Fraction(fours, 16) + sum(Fraction(1, b*b) for b in floors)
                if mass + Fraction(remaining, 36) >= 1:
                    yield fours, primes, choices


def gauss_errors(rank, profiles):
    for fours, primes, choices in profiles:
        assert all(n == p or b == p for p, (n, b, _, _) in zip(primes, choices))
        present = tuple(p for p, (_, b, _, _) in zip(primes, choices) if p == b)
        assert present
        absent = tuple(p for p in primes if p not in present)
        required = 4 * prod(p if p in present else p*p for p in primes)
        proper = max(n for _, _, n, _ in choices)
        packets = [(1,) if p == 5 else range(1, p) for p in present]
        for counts in product(*packets):
            classes = [p if j == 1 else (p - 1)*j for p, j in zip(present, counts)]
            budget = rank - 1 - fours - sum(classes) - proper
            if budget < sum(p - 1 for p in absent):
                continue
            available_mass = 1 - Fraction(fours, 16) - sum(
                Fraction(n, p*p) for n, p in zip(classes, present))
            for trace in range(-fours, fours + 1):
                residual = 1 - Fraction(fours - trace, 16) - sum(
                    Fraction(j, p) for j, p in zip(counts, present))
                if -Fraction(1, required) < residual <= 2*available_mass + Fraction(1, required):
                    yield from sparse_budget(residual, required, absent, budget)


F = Fraction

def square_partitions(n,least=1):
    if n==0:
        yield ()
        return
    for a in range(least,isqrt(n)+1):
        for tail in square_partitions(n-a*a,a):
            yield (a,)+tail


def verify():
    # The three-packet union lemma: possible coefficient supports within
    # five coordinates. A support needs at least three distinct degrees.
    shapes={p:[x for x in square_partitions(p-1) if 3<=len(x)<=5]
            for p in (5,7,11)}
    assert shapes=={5:[(1,1,1,1)],7:[(1,1,2)],11:[(1,1,2,2)]}
    # The underlying two-prime intersection statement is exhaustive.
    sets=[1,2,3]
    assert all((a&b&c)!=0 for a,b,c in product(sets,repeat=3)
               if a&b and a&c and b&c)

    rows=list(cases(44))
    pure=[r for r in rows if all(n==p or b==p
          for p,(n,b,_,_) in zip(r[1],r[2]))]
    mixed=[r for r in rows if r not in pure]
    assert (len(rows),len(pure),len(mixed))==(49,38,11)
    assert not list(gauss_errors(44,pure))

    survivors=[]
    for m,primes,choices in mixed:
        assert choices[0][:2]==(8,6)
        room=44-1-m-sum(n for n,_,_,_ in choices)
        present=tuple(p for p,(_,b,_,_) in zip(primes,choices) if b==p)
        if room==3:
            assert 11 in present and 5 not in present
            assert not any(sum(a*a for a in x)==10
                           for x in product(range(1,4),repeat=3))
            continue
        if room==4:
            assert 5 in present and any(p>5 for p in present)
            continue # Four signs at5 and even norm at7 require one more label.
        assert primes==(3,5,7,11)
        assert choices[1][1]==5 and choices[2][1]==7
        survivors.append((m,present,room))
    assert survivors==[(6,(5,7),6),(7,(5,7),5),
                      (4,(5,7,11),8),(5,(5,7,11),7),
                      (6,(5,7,11),6),(7,(5,7,11),5)]

    # Minimum packets: another packet uses at least p-1 additional labels,
    # leaving fewer than the four mandatory proper five-free labels.
    for m,present,room in survivors:
        assert all(room-(p-1)<4 for p in present if p>5)

    P=F(1,5)+F(1,7)+F(1,11)
    W5=2*F(1,400)+F(1,900)+F(1,2025)+F(1,441)
    W6=W5+F(1,441)
    vac=F(1,2310**2)
    upper={
      'absent11_m6':F(6,16)+F(8,36)+F(1,5)+F(1,7)+F(11,121**2)+W5+F(1,36)+vac,
      'absent11_m7_gap':F(7,16)+F(13,72)+F(1,5)+F(1,7)+F(11,121**2)+W5+vac,
      'present_m4':F(4,16)+F(8,36)+P+W6+F(2,36)+vac,
      'present_m5_at_most8sixes':F(5,16)+F(8,36)+P+W6+F(1,64)+vac,
      'present_m6_gap':F(6,16)+F(14,81)+P+W6+vac,
    }
    assert all(x<1 for x in upper.values())
    lower={
      'absent11_m7_eight_sixes':F(7,16)+F(8,36)+F(1,5)+F(1,7)+F(11,121**2),
      'present_m5_ten_sixes':F(5,16)+F(10,36)+P,
      'present_m6_eight_sixes':F(6,16)+F(8,36)+P,
    }
    assert all(x>1 for x in lower.values())
    gauss=144*(1-F(5,16)-F(9,36)-P+F(1,2310))
    assert gauss==F(3,5)<1
    assert 1+7+5+7+11+8+6==45 # present_m7 has too many labels.

    # With at least five odd primes, the least isolated reservation is
    # 3+5+7+11+13=39. Through47 only3 can be mixed with two.
    # The sole rank-44 reservation with a dyadic interaction at p>=7.
    assert F(24,196)+F(1,27)+F(1,5)+F(1,11)+vac<1
    upper_six={
      'all_isolated':F(1,27)+sum(F(1,p) for p in (5,7,11,13))+F(7,16)+vac,
      'mixed_three':F(8,36)+sum(F(1,p) for p in (5,7,11,13))+F(2,16)+vac,
    }
    assert all(x<1 for x in upper_six.values())
    return dict(rank=44,initial_profiles=len(rows),all_primary=len(pure),
                mixed=len(mixed),remaining_mixed_profiles=len(survivors),survivors=0,
                upper={k:str(v) for k,v in upper.items()},
                lower={k:str(v) for k,v in lower.items()},gauss_house=str(gauss),
                at_least_six_primes_through47={k:str(v) for k,v in upper_six.items()})


if __name__=='__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
