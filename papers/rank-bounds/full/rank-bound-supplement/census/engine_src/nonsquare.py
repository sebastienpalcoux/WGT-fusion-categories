"""Complete nonsquare starts and unbounded continuations at ranks 1--44.

Native build: c++ -O3 -std=c++17 engine_src/nonsquare.cpp -o BUILD/nonsquare
Native command: BUILD/nonsquare R Q B M [--prefix comma,separated,denominators]
Here M in the command is the top multiplicity, not the dimension root.
Standard output is newline-delimited degree lists or objects with ``exact``
and ``primes`` keys.  Every exact object covers its entire prefix subtree;
the Python continuation has no search cap.  Successful standard error ends
with: rank R, q Q, B B, multiplicity M, nodes N, candidates C, exact E.

The shared engine directory supplies arithmetic helpers only.  In particular
we do not invoke the obsolete rank<=34 complete() functions in nonsquare_exact.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt

from defect import defect_tiles
from denominators import smooth_interval, smooth_over, terminal
from nonsquare_exact import (joint_packets, possible_packets,
                             remaining_mass_possible, seed_mass_possible)
from primes import primes_up_to
from screen import factor, fiber_floor
from starts import split_prefixes


RANKS = range(1, 45)
OBSTRUCTIONS = {'frontier', 'four_prime', 'effective_prime_height'}


def _check_rank(rank):
    if type(rank) is not int or rank not in RANKS:
        raise ValueError('the local runner supports exact ranks1 through44 only')


@lru_cache(None)
def start_pairs(rank):
    """All nonsquare, nonsquarefree maximal codegrees 8<=f<=rank-2.

    Codegree eight requires rank>=46 by the new theoretical exclusion.
    Squarefree nonunit codegrees are impossible.  Square dimension is owned
    by the separate square generator.
    """
    _check_rank(rank)
    pairs = []
    for f in range(8, rank-1):
        b = q = 1
        for p, exponent in factor(f):
            b *= p**(exponent//2)
            if exponent % 2:
                q *= p
        if b > 1 and q > 1 and f != 8:
            pairs.append((b, q))
    return tuple(sorted(pairs))


@lru_cache(None)
def supports(b, q, rank):
    """Every support permitted by the ramified fiber and joint-packet bounds."""
    if (b, q) not in start_pairs(rank):
        raise ValueError('inadmissible nonsquare starting pair')
    required = {p for p, _ in factor(2*b*q)}
    fiber = fiber_floor(q)
    inactive = max(4, fiber) if q == 3 else fiber
    budget = ((rank-inactive)//fiber
              - sum(p-1 for p in required if p != 2 and q % p))
    if budget < 0:
        return ()
    optional = [p for p in primes_up_to(rank)
                if p not in required and p-1 <= budget]
    minimum = 4 if q*b*b == 18 else 3
    result = []
    for size in range(len(optional)+1):
        for group in combinations(optional, size):
            if (len(required)+size >= minimum
                    and sum(p-1 for p in group) <= budget):
                result.append(tuple(sorted(required.union(group))))
    return tuple(result)


@lru_cache(None)
def multiplicities(rank, b, q):
    """Rank-independent top-class bounds followed by an optimistic mass test."""
    if (b, q) not in start_pairs(rank):
        raise ValueError('inadmissible nonsquare starting pair')
    f = q*b*b
    upper = min([f-1] + [f-1-(p-1)//2 for p, exponent in factor(f)
                         if p > 2 and exponent == 1])
    if f == 18:
        upper = min(upper, 10)
    step = 1
    if b == 2 and q % 2 and factor(q) == [(q, 1)]:
        step = (q-1)//2
        upper = min(upper, 7*step)
    fiber = fiber_floor(q)
    result = []
    for m in range(step, upper+1, step):
        if q % 2 == 0 and b % 2 and m % 2:
            continue
        if any(seed_mass_possible(q, b, m, rank, ps, fiber)
               for ps in supports(b, q, rank)):
            result.append(m)
    return tuple(result)


def tasks(rank):
    """A complete disjoint breadth partition, with no historical input files."""
    _check_rank(rank)
    result = []
    for b, q in start_pairs(rank):
        for support in supports(b, q, rank):
            if (len(support) > 5
                    or sum(p > 2 and q % p != 0 for p in support) > 4):
                raise ValueError('native fixed-size support counters would overflow')
        for m in multiplicities(rank, b, q):
            task = dict(q=q, B=b, m=m)
            if q == 3 and b == 2 and rank >= 36 and m >= 4:
                for prefix in _deep_split_prefixes(rank, b, q, m):
                    result.append(dict(task, prefix=list(prefix)))
            elif m >= 8 or q == 2:
                for prefix in split_prefixes(rank, b, q, m, goal=64):
                    result.append(dict(task, prefix=list(prefix)))
            else:
                result.append(task)
    groups = {}
    for task in result:
        groups.setdefault((task['B'], task['q'], task['m']), []).append(
            task.get('prefix', [task['B']]))
    for prefixes in groups.values():
        ordered = sorted(prefixes)
        if any(y[:len(x)] == x for x, y in zip(ordered, ordered[1:])):
            raise ValueError('overlapping nonsquare task prefixes')
    return result


def _deep_split_prefixes(rank, b, q, m, goal=64):
    """Disjoint elementary frontier without the historical nine-level limit.

    This is the same all-children-or-parent breadth split as starts.py.
    It does not consult any stored search output or new obstruction.  The
    budget controls task granularity only; an unsplit parent remains a task.
    """
    queue = deque([((b,)*m, b, q*b*b-m, rank-1-m)])
    settled = []
    while queue:
        prefix, L, E, slots = queue.popleft()
        if slots <= 2:
            settled.append(prefix)
            continue
        lo = max(prefix[-1], isqrt(L*L//E)+1, b+1)
        hi = isqrt((36*slots+1)*L*L//(36*E))
        if hi < lo:
            continue
        if len(settled)+len(queue)+hi-lo+1 > goal:
            settled.append(prefix)
            continue
        for z in range(lo, hi+1):
            g = gcd(L, z)
            k = z//g
            queue.append((prefix+(z,), L*k, E*k*k-(L//g)**2, slots-1))
    return sorted(settled)


def _effective_prime_heights(M, counts, q, primes, suffix=()):
    """The twice-square prime-denominator theorem, valid through rank59."""
    return q != 2 or not any(
        p >= 5 and M % (p*p) == 0 and (counts[p] or p in suffix)
        for p in primes)


def _three_prime_block_capacity(slots, prefix, counts, ps, q):
    """Necessary partner capacity, including an optimistic open-class tail.

    At fixed final support {2,3,p}, the prime-free trace theorem forbids
    equal ternary blocks at denominator 2 or p. A denominator-p_b label
    must be the larger member of (n/2,n/2,n), using two denominator-2p_b
    partners. For p_b=2 the other orientation would need denominator 1,
    which is excluded by the squarefree nonunit-codegree obstruction.
    Until 2p_b is passed, all remaining positions are allowed as partners.
    """
    if q != 3 or len(ps) != 3 or ps[:2] != (2, 3):
        return True
    return all(2*counts[p] <= counts[2*p]
               + (slots if prefix[-1] <= 2*p else 0) for p in (2, ps[2]))


def exact_rows(item, task, rank):
    """Continue one native width handoff with unlimited-size integer arithmetic."""
    _check_rank(rank)
    b, q, m = task['B'], task['q'], task['m']
    if (b, q) not in start_pairs(rank) or m not in multiplicities(rank, b, q):
        raise ValueError('inadmissible nonsquare task')
    key = task.get('prefix', [b])
    raw = item.get('exact')
    if (not isinstance(raw, list) or raw[:len(key)] != key
            or not 1 <= len(raw) <= rank-3
            or any(type(v) is not int or v <= 0 for v in raw)
            or raw != sorted(raw)):
        raise ValueError('invalid arbitrary-precision prefix')
    prefix = list(raw)
    ps = tuple(item.get('primes', ()))
    if ps not in supports(b, q, rank) or not all(smooth_over(v, ps) for v in prefix):
        raise ValueError('invalid nonsquare handoff support')
    L, E = 1, q
    for z in prefix:
        g = gcd(L, z)
        scale = z//g
        L, E = L*scale, E*scale**2-(L//g)**2
        if E <= 0:
            return
    if prefix.count(b) > m or (prefix[-1] > b and prefix.count(b) != m):
        return
    fiber = fiber_floor(q)
    tile_primes = tuple(p for p in (2, 3, 5) if q % p == 0)

    def visit(L, E, slots):
        counts = Counter(prefix)
        last = prefix[-1]
        if not _three_prime_block_capacity(slots, prefix, counts, ps, q):
            return
        if not _effective_prime_heights(L, counts, q, ps):
            return
        if q == 2 and counts[3] and (counts[3] > 10 or L % 9 == 0):
            return
        if not possible_packets(prefix, slots, ps, q, fiber):
            return
        if not all(defect_tiles(counts, p, last=last, slots=slots)
                   for p in tile_primes):
            return
        lo = max(last, isqrt(L*L//E)+1)
        hi = isqrt((36*slots+1)*L*L//(36*E))
        if last == b:
            if len(prefix) < m:
                hi = min(hi, b)
            else:
                lo = max(lo, b+1)
        if lo > hi or not remaining_mass_possible(
                L, E, lo, slots, prefix, ps, q, fiber):
            return
        if slots == 2:
            for z, k, d, e in terminal(L, E, lo, hi, ps, fixed_support=True):
                M = k*L
                suffix = (z, M//d)
                if not _effective_prime_heights(M, counts, q, ps, suffix):
                    continue
                full_prefix = prefix+list(suffix)
                if not _three_prime_block_capacity(0, full_prefix,
                        Counter(full_prefix), ps, q):
                    continue
                row = [1, d, e] + [M//w for w in reversed(prefix)]
                if row.count(row[-1]) != m:
                    continue
                if any(sum(x % p == 0 for p in ps) < 2 for x in row[1:]):
                    continue
                if any(sum(x % p != 0 for x in row[1:])
                       < (3 if p in (2, 5) else 2) for p in ps):
                    continue
                if joint_packets(full_prefix, ps, q, fiber):
                    yield row
            return
        for z in smooth_interval(lo, hi, ps):
            g = gcd(L, z)
            scale = z//g
            prefix.append(z)
            yield from visit(L*scale, E*scale**2-(L//g)**2, slots-1)
            prefix.pop()

    yield from visit(L, E, rank-1-len(prefix))


def row_failure(row, task, rank):
    """Extra nonsquare rules; terminal.py owns validation and shared screens."""
    _check_rank(rank)
    b, q = task['B'], task['q']
    M = b*row[-1]
    counts = Counter(M//d for d in row[1:])
    if not _effective_prime_heights(M, counts, q, (p for p, _ in factor(M))):
        return 'effective_prime_height'
    if b == 2 and q % 2 and factor(q) == [(q, 1)]:
        half = (q-1)//2
        if task['m'] % half or task['m'] > 7*half:
            return 'four_prime'
    if q == 2 and M % 9 == 0 and counts[3]:
        return 'frontier'
    if q == 3 and counts[2] and counts[4] < 2*counts[2]:
        return 'frontier'
    if q == 5 and M % 9 == 0 and counts[3] > counts[6]:
        return 'frontier'
    return None
