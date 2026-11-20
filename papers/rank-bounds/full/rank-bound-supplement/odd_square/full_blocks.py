"""Residue and omission bounds for prime-free neutral labels."""
from fractions import Fraction
from functools import lru_cache
from math import ceil, lcm


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def seven_free(nines,marked):
    """Use the mass of rows of other prime-power codegrees."""
    other=Fraction(nines,81)+sum(Fraction(v[1],p*p) for p,v in marked.items() if p!=7)
    return 8 if other>Fraction(3,7) else 6


@lru_cache(None)
def count(q, nines):
    """Lower bound from the joint residue spectrum and the dimension modulo q."""
    half = (q - 1) // 2
    floors = []
    for residue in range(1, half + 1):
        upper = next(n for n in range(3, 3*q, 3) if (n + residue) % q == 0)
        lower = upper - 3*q
        costs = []
        for count_upper in range(0, nines + 1, 2):
            trace = (nines - count_upper)*lower + count_upper*upper
            squares = (nines - count_upper)*lower**2 + count_upper*upper**2
            mean = Fraction(-trace, 81 - nines)
            left = -residue + q*((mean + residue)//q)
            number = squares - (2*left + q)*trace - (81 - nines)*left*(left + q)
            costs.append(Fraction(number, 162))
        floors.append(ceil(min(costs)))
    for pairs in range(seven_free(nines,{})//2 if q==7 else 1, 32):
        for multiplicities in compositions(pairs, half):
            if (1 + 2*sum(m*t*t for t,m in enumerate(multiplicities, 1))) % q:
                continue
            if all(sum(m*min(a*t % q, -a*t % q)**2
                       for t,m in enumerate(multiplicities, 1)) >= lower
                   for a,lower in enumerate(floors, 1)):
                return 2*pairs
    return 64  # All smaller numbers of pairs have been excluded.

def neutral_pair_floors(support, marked, parts, nines, nine=False):
    """Denominator floors of disjoint neutral pairs on four or five primes.

    parts contains certified prime-power divisors at a search node, or
    the exact primary parts when the dimension root is known.
    """
    isolated=[p for p in support[2:] if p not in marked]
    if len(support)==5:
        if 7 not in isolated:
            raise ValueError("unhandled five-prime neutral support", support, marked)
        # Least-residue bound for seven-free labels uses all coprime-codegree
        # rows; the recorded nine and prime rows exceed mass 3/7 here.
        divisible=[parts[7]]*(seven_free(nines,marked)//2)
        generic=0
    elif len(marked)==2:
        p,q=sorted(marked)
        divisible=[parts[3]*p,parts[5]*p,parts[3]*q,parts[5]*q,p*q]
        generic=0
    elif marked:
        p,=marked
        q,=isolated
        divisible=[parts[3]*p,parts[5]*p,parts[q]*p,parts[q]]
        if q==7:
            divisible.extend([parts[q]]*(seven_free(nines,marked)//2-2))
        generic=0
    else:
        if 7 in isolated:
            divisible=[parts[7]]*(seven_free(nines,marked)//2)
        else:
            divisible=[min(parts[p] for p in isolated)]*2
        generic=0
    if nine:
        return [lcm(9,d) for d in divisible]+[27]*generic
    return divisible+[15]*generic


@lru_cache(None)
def omission_count(p, nines):
    """Deletion supplies two nonunit dual pairs for every larger prime."""
    return count(p, nines) if p <= 13 else 4
