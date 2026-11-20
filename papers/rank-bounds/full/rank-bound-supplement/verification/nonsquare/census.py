#!/usr/bin/env python3
"""Exact four-prime ternary census through rank 211. Standard library only."""
from bisect import bisect_left
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil, isqrt, lcm, prod
import json

LIMIT = 211
PRIMES = (5, 7, 11, 13, 17, 19, 23)


def upper(slots, vacuum, demands):
    remaining = slots - sum(count for count, _ in demands)
    if remaining < 0:
        return -1
    return vacuum + Q(remaining, 225) + sum(
        Q(count, max(15, floor)**2) for count, floor in demands)


@lru_cache(None)
def signatures(m):
    """Real cube discrepancy and oriented cube count of a ternary degree class."""
    return sorted({(2*neutral-primitive, 2*positive-primitive)
                   for higher in range(0, m+1, 6)
                   for neutral in range(0, m-higher+1, 2)
                   for primitive in [m-higher-neutral]
                   for positive in range(0, primitive+1, 2)})


@lru_cache(None)
def signed_blocks(p, m):
    if m == 6*(p-1):
        return (-Q(4, p*p), Q(4, p*p))
    if p in (7, 11) and m == 6*p:
        return (-Q(8,p*p), 0, Q(8,p*p))
    if m < 6*(p-1):
        orbit = (p-1) & -(p-1)
        return tuple(Q(2*i, p*p) for i in range(-m//3, m//3+1, 2*orbit))
    limit = (m//3)*(p-1)
    return tuple(Q(2*i,(p-1)*p*p) for i in range(-limit,limit+1))


def census(ps, rank=LIMIT):
    ds = sorted({3, 9, *ps})
    base = 3*prod(ps)
    scale = lcm(base, *ds, 15)**2
    weights = {b: scale//(b*b) for b in {15, base, *ds}}
    vacuum = Q(1, base*base)
    quota = {p: 6*(p-1) for p in ps}
    assert all(30*(p-1) + sum(quota[q] for q in ps if q != p) > rank for p in ps)
    counts = [0, 0, 0]
    survivors = []

    @lru_cache(None)
    def possible_full(p, m):
        return any(
            j*(p-1) <= m//3 and
            (j != 2 or (p >= 11 and p % 4 == 3 and m//3 == 2*p)) and
            p*abs(p*j-m//3) <= j*(p-j)+Q((p+1)*j*p, base)
            for j in range(2, p, 2))

    def check_height(m, slots, energy, low, activity, g, signs):
        if g == 3 and m[9]:
            return
        balanced = g == 27 or (g == 9 and (m[9] <= 2 or m[3] > 3))
        if balanced and m[3] % 6:
            return
        if g == 9 and energy < sum(Q(c, b*b) for c,b in activity if b//3 in low):
            return
        heights = {p: p if p in low or not low else p*p for p in ps}
        heights[3] = g
        lower_m = prod(heights.values())
        vacuum = Q(1, lower_m**2)
        witnesses = []
        free = sum(m[b] for b in ds if b % g == 0)
        for a, b in combinations((3,)+ps, 2):
            required = (a in low or b in low or (a == 3 and
                        (g == 3 or (m[3] and balanced))))
            if not required:
                continue
            count = 6 if a != 3 and g == 3 else 2
            witnesses.append((count, heights[a]*heights[b], a != 3))
            if a == 3:
                free += count
        if g == 9 and m[3]:
            choices = (9,15,21) if m[3] > 3 else (9,15,21,27)
            full_free = next((n for n in choices if n >= 1+free), None)
            if full_free is None:
                return
            witnesses.append((full_free-1-free,45,False))
        else:
            witnesses.append((max(0,8-free), g, False))

        spare = slots - sum(c for c, _ in activity) - sum(c for c,_,_ in witnesses)
        reserve = [(c, 3*b if allow_ternary and g > 3 and spare < 4 else b)
                   for c,b,allow_ternary in witnesses]
        if energy > upper(slots, vacuum, activity+reserve):
            return
        counts[2] += 1
        neutral = upper(slots-sum(c for c, _ in activity), vacuum, reserve)
        weighted = sum(Q(c, (b//3-1)*b*b) for c,b in activity)
        for d,u in ([(0,0)] if balanced else signatures(m[3])):
            for e,v in signatures(m[9]):
                z = Q(9*d+e,81)
                if not -2*neutral <= z <= neutral:
                    continue
                error = Q(2,lower_m) + min(energy, (2*(neutral+weighted)+z)/3)
                target = -Q(9*u+v,81)
                i = bisect_left(signs, target-error)
                if i < len(signs) and signs[i] <= target+error:
                    survivors.append((tuple(m[b] for b in ds), g))
                    return

    def accept(ms, slots, energy):
        counts[0] += 1
        m = dict(zip(ds, ms))
        low = {p for p in ps if m[p]}
        if len(low) >= 2:
            return
        activity = []
        for p in ps:
            if m[p] < quota[p]:
                activity.append((max(quota[p], ceil(Q(9*(p-1)*m[p],2*p))), 3*p))
            elif not possible_full(p, m[p]):
                activity.append((quota[p], 3*p))
        if energy > upper(slots, vacuum, activity):
            return
        counts[1] += 1
        if (m[3]>=8 and m.get(5,0)>=36) or (m[3]>=6 and m.get(5,0)>=48):
            return
        assert m[3] or m[9]
        for p in ps:
            if m[p] == quota[p]:
                assert (1+m[3]+m[9]+m[p]+ceil(Q(9*(p-1)*m[p],2*p))
                        +sum(quota[q] for q in ps if q!=p) > rank)
        signs = sorted({sum(v) for v in product(*(signed_blocks(p,m[p]) for p in ps))})
        for g in (3,9,27):
            check_height(m, slots, energy, low, activity, g, signs)

    def visit(i, ms, slots, residual):
        if residual <= 0:
            return
        floor = min(15, ds[i]) if i < len(ds) else 15
        if residual > weights[base] + slots*weights[floor]:
            return
        if i == len(ds):
            accept(ms, slots, Q(residual, scale))
            return
        b = ds[i]
        step = 2 if b % 3 == 0 else 3*((b-1) & -(b-1))
        maximum = min(slots, 3*b*b-1, {3:16,5:48}.get(b, slots))
        for m in range(0, maximum+1, step):
            if b in ps and m == 3*b*(b-1):
                continue
            visit(i+1, ms+(m,), slots-m, residual-m*weights[b])
    visit(0, (), rank-1, 3*scale)
    return counts, survivors


def sector_denominators(m9):
    """Enclose every B_p above 500 by 0 <= 1/B_p^2 <= 1/501^2."""
    M0 = 27*5*49*121
    energy = 3-Q(12,9)-Q(24,25)-Q(m9,81)
    slots = 54-m9
    reserves = [(2,b) for b in (135,1323,3267,27,245,605)]
    neutral = upper(slots,Q(1,M0**2),reserves)
    primes = (5,7,11)
    options = []
    for p in primes:
        values = sorted({3**a*p**j for a in range(6) for j in range(1,5)
                         if (p != 5 or j == 1) and 3*p <= 3**a*p**j <= 500})
        options.append([(b,Q(1,b*b),Q(1,b*b)) for b in values] + [(None,0,Q(1,501**2))])
    possible = []
    for triple in product(*options):
        least = sum(6*(p-1)*lo for p,(_,lo,_) in zip(primes,triple))
        most = sum(6*(p-1)*hi for p,(_,_,hi) in zip(primes,triple))
        if least < energy <= most+neutral:
            possible.append(tuple(b for b,_,_ in triple))
    return possible


def terminal_tail():
    """Enumerate nine unordered dual pairs, preserving every possible witness matching."""
    primes = (3,5,7,11)
    minimum_m = 27*5*49*121
    floors = (135,1323,3267,27,245,605)
    nodes = [0]*10
    solutions = []

    def denominators(lo, hi):
        def powers(i, value):
            if i == len(primes):
                if value>=lo:
                    yield value
                return
            p = primes[i]
            while value<=hi:
                yield from powers(i+1,value)
                value *= p
                if p == 5 and value%25 == 0:
                    break
        return sorted(powers(0,1))

    def search(residual, lo, remaining, matchings, fixed):
        nodes[remaining] += 1
        if remaining == 0:
            if residual.numerator == 1 and 63 in matchings:
                root = isqrt(residual.denominator)
                if root*root == residual.denominator and root % lcm(minimum_m,*fixed) == 0:
                    solutions.append((fixed,root))
            return
        hi = isqrt(((2*remaining+1)*residual.denominator)//residual.numerator)
        lo = max(lo, 1+isqrt(2*residual.denominator//residual.numerator))
        for b in denominators(lo,hi):
            next_residual = residual-Q(2,b*b)
            candidates = set(matchings)
            for state in matchings:
                for i,floor in enumerate(floors):
                    if not (state>>i)&1 and b%floor == 0:
                        candidates.add(state | (1<<i))
            keep = []
            for state in candidates:
                missing = [val for i,val in enumerate(floors) if not (state>>i)&1]
                slack = remaining-1-len(missing)
                if slack < 0:
                    continue
                cap = (Q(1,max(b,minimum_m)**2) + Q(2*slack,b*b) +
                       sum(Q(2,max(b,val)**2) for val in missing))
                if cap >= next_residual:
                    keep.append(state)
            if keep:
                search(next_residual,b,remaining-1,keep,fixed+(b,))
    residual = 3-Q(12,9)-Q(36,81)-Q(24,25)-Q(24,225)-Q(36,441)-Q(60,1089)
    assert residual == Q(5023,266805)
    search(residual,15,9,[0],())
    return nodes,solutions


def reductions():
    # r<=218, q squarefree and 9*q<r.
    assert 135+120*(1-Q(9,49)) > 218
    assert 189+112*(1-Q(9,25)) > 218
    assert 171+152*(1-Q(9,25)) > 218
    assert 207+184*(1-Q(9,25)) > 218
    assert min(7*(1+2*(2+4+10)),10*(1+2*(2+4+6)))>218
    assert 12*(1+Q(1600,89)) > 218 and 6*(1+Q(7056,193)) > 218
    # q=5: possible rank at most 213, with at most eight more slots than at 205.
    assert Q(62,25)+Q(151,81)<5
    assert Q(112,25)+Q(60,441)+Q(20,81)+Q(12,625)+Q(4,2025)+Q(4,60025)+Q(1,11025**2)<5
    assert Q(20,9)+Q(20,81)+Q(60,441)+Q(100,121)+Q(12,25)+Q(1,1155**2)<5
    assert 1+12*4+6*(6+10+12) == 217
    assert Q(18,9)+Q(24,225)+Q(36,441)+Q(60,121)+Q(72,1521)+Q(1,15015**2)<3


def main():
    reductions()
    supports = [ps for ps in combinations(PRIMES,3) if 3+sum(6*(p-1) for p in ps)<=LIMIT]
    report = {}
    terminal = []
    for ps in supports:
        counts, rows = census(ps)
        report[','.join(map(str,ps))] = [*counts,len(rows)]
        terminal.extend((ps,values,g) for values,g in rows)
    expected = [((5,7,11),(12,24,0,m,0),27) for m in (36,38,42)]
    assert sorted(terminal) == expected, terminal
    triples = {m:sector_denominators(m) for m in (36,38,42)}
    assert triples == {36:[(15,21,33)],38:[],42:[]}, triples
    nodes, solutions = terminal_tail()
    assert not solutions
    print(json.dumps({'rank_cap':LIMIT,'counts':report,
                      'height_survivors':terminal,'terminal_denominators':triples,
                      'tail_nodes_by_remaining_pairs':nodes,'solutions':solutions},indent=2))


if __name__ == '__main__':
    main()
