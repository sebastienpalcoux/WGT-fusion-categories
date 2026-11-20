"""Prime-free census from the integral spectrum on all character rows."""
from fractions import Fraction
from functools import lru_cache
from math import ceil


def norm_floor(p, a, m):
    """m denominator-four rows have even spectrum; all rows are -a modulo p."""
    f = 16
    upper = next(z for z in range(1, 2*p) if z % 2 == 0 and (z+a) % p == 0)
    lower = upper - 2*p
    result = []
    for j in range(m+1):
        trace = (m-j)*lower+j*upper
        squares = (m-j)*lower**2+j*upper**2
        mean = Fraction(-trace, f-m)
        left = -a + p*((mean+a)//p)
        right = left+p
        result.append(Fraction(squares-(left+right)*trace-(f-m)*left*right, f))
    return ceil(min(result))


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for k in range(total+1):
        for tail in compositions(total-k,parts-1):
            yield (k,)+tail


@lru_cache(None)
def census(p, m):
    half = (p-1)//2
    floors = [norm_floor(p,a,m) for a in range(1, half+1)]
    for total in range(1, 100):
        for counts in compositions(total, half):
            if (1+sum(n*j*j for j,n in enumerate(counts,1))) % p:
                continue
            if all(sum(n*min((a*j)%p,(-a*j)%p)**2 for j,n in enumerate(counts,1)) >= need
                   for a,need in enumerate(floors,1)):
                return total
    raise ArithmeticError('prime-free census not reached')
