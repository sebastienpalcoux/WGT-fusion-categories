"""Provable starting cases and disjoint worker tasks for ranks nine to thirty-four."""
from collections import deque
from math import gcd, isqrt


def start_pairs(rank):
    """(B,q), where q is square-free and B=M/d_max, as in the five-start reduction."""
    if type(rank) is not int or not 9 <= rank <= 34:
        raise ValueError('supported ranks are 9 through 34')
    pairs = ((4, 1), (5, 1), (2, 3), (2, 5), (3, 2))
    return tuple((b,q) for b,q in pairs if q*b*b <= rank-2)


def top_multiplicities(rank, b, q):
    """The multiplicity restrictions before testing any individual prime support."""
    if (b,q) not in start_pairs(rank) or q == 1:
        raise ValueError('expected an admissible nonsquare start')
    f = q*b*b
    # The mass at trivial and order-p phases for each odd p of valuation one
    # gives m <= f - 1 - (p-1)/2. At 2, odd b comes in pairs.
    bound = min([f-1] + [f-1-(p-1)//2 for p in (3,5) if f%p == 0 and f%(p*p)])
    return range(2,bound+1,2) if q == 2 and b%2 else range(1,bound+1)


def split_prefixes(rank, b, q, m, goal=48):
    """Group a breadth-first exact denominator frontier into at most goal siblings/leaves.

    Begin with the fixed top run. Expanding a node replaces it by *all*
    children permitted by the elementary residual interval. If that would
    overfill the budget, retain the whole parent. Empty intervals prove that
    a node has no completions. This uses no native pruning assumptions.
    """
    initial = (b,)*m
    slots = rank-1-m
    queue = deque([(initial,b,q*b*b-m,slots)])
    settled = []
    while queue:
        prefix, L, E, slots = queue.popleft()
        if slots <= 2 or len(prefix) >= m+9:
            settled.append(prefix)
            continue
        lo = max(prefix[-1], isqrt(L*L//E)+1, b+1)
        hi = isqrt((36*slots+1)*L*L//(36*E))
        if hi < lo:
            continue
        if len(settled)+len(queue)+hi-lo+1 > goal:
            settled.append(prefix)
            continue
        for z in range(lo,hi+1):
            g = gcd(L,z)
            scale = z//g
            queue.append((prefix+(z,), scale*L, scale*scale*E-(L//g)**2, slots-1))
    return sorted(settled)


def nonsquare_tasks(rank):
    """Keep the same multiplicity coordinates at every rank; split expensive higher-rank jobs."""
    tasks = []
    for b,q in start_pairs(rank):
        if q == 1:
            continue
        for m in top_multiplicities(rank,b,q):
            task = dict(B=b,q=q,m=m)
            # The largest ternary branches dominate the nonsquare search at
            # the next ranks. Split the expensive, high-multiplicity branches
            # without changing the arithmetic search or earlier partitions.
            if rank > 32 and ((q == 3 and m >= 8) or (q == 2 and m >= 12)
                              or (q == 5 and m >= 13)):
                for prefix in split_prefixes(rank,b,q,m):
                    tasks.append(dict(task,prefix=list(prefix)))
            else:
                tasks.append(task)
    return tasks
