"""Conservative unbounded square continuation for ranks below 44.

This is the square branch of the former frontier continuation, expressed
without importing any historical driver, archive, receipt or partition.
Every native width handoff is exhausted with Python integers. The final
row classifier owns all further necessary-condition tests.
"""
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt

from denominators import smooth_interval, smooth_over, terminal
from primes import primes_up_to
from square_exact import packet_room
from prime_cut_prefix import prime_cut_prefix_failure


def complete(item, task, rank):
    """Yield every completion of one native square-prefix handoff."""
    if type(rank) is not int or not 1 <= rank <= 43 or task.get('q') != 1:
        raise ValueError('this square continuation supports ranks 1 through 43')
    key = task.get('prefix', [task['B']])
    raw = item.get('exact')
    if (not isinstance(raw, list) or raw[:len(key)] != key
            or not 1 <= len(raw) <= rank - 3
            or any(type(v) is not int or v <= 0 for v in raw)
            or raw != sorted(raw)):
        raise ValueError('invalid arbitrary-precision square prefix')
    prefix = list(raw)
    all_primes = set(primes_up_to(rank - 1))
    if item.get('primes') or not all(smooth_over(v, all_primes) for v in prefix):
        raise ValueError('invalid square handoff support')
    L, E = 1, 1
    for z in prefix:
        g = gcd(L, z)
        scale = z // g
        L, E = L * scale, E * scale**2 - (L // g)**2
        if E <= 0:
            return

    def visit(L, E, slots):
        counts = Counter(prefix)
        last = prefix[-1]
        seen = {p for p in all_primes if L % p == 0}
        if not packet_room(counts, last, seen, slots):
            return
        allowed = seen | {2} | {p for p in all_primes if p <= slots + 1}
        # Preserve only the former continuation's explicitly proved
        # through-rank40 high-four restrictions. Omitting the stronger
        # rank41--43 native pruning here safely retains extra branches.
        if counts[4] >= 9 and rank <= 40:
            if (counts[4] > 12 or not seen <= {2, 3, 5, 7}
                    or (counts[4] <= 11 and not counts[5]
                        and last > 7 and not counts[7])
                    or ((counts[4] == 12 or counts[4] >= 11 and counts[5])
                        and counts[7])):
                return
            allowed &= {2, 3, 5, 7}
        lo = max(last, isqrt(L * L // E) + 1)
        hi = isqrt((36 * slots + 1) * L * L // (36 * E))
        if lo > hi:
            return
        if prime_cut_prefix_failure(prefix, slots, tuple(sorted(allowed)), lo,
                Fraction(E, L*L), required_primes=seen,
                vacuum_upper=Fraction(1, max(L*L, 36*lo*lo))) is not None:
            return
        if slots == 2:
            for z, k, d, e in terminal(L, E, lo, hi, allowed, fixed_support=False):
                M = k * L
                row = [1, d, e] + [M // w for w in reversed(prefix)]
                actual = tuple(p for p in all_primes if M % p == 0)
                if any(sum(x % p == 0 for p in actual) < 2 for x in row[1:]):
                    continue
                if any(sum(x % p != 0 for x in row[1:])
                       < (3 if p in (2, 5) else 2) for p in actual):
                    continue
                yield row
            return
        for z in smooth_interval(lo, hi, allowed):
            g = gcd(L, z)
            scale = z // g
            prefix.append(z)
            yield from visit(L * scale, E * scale**2 - (L // g)**2, slots - 1)
            prefix.pop()

    yield from visit(L, E, rank - 1 - len(prefix))
