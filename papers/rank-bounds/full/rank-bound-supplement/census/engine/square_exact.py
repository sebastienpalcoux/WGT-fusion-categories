"""Exact continuation of a square-dimension denominator prefix, using Python integers."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, lcm

from denominators import smooth_interval, smooth_over, terminal
from primes import primes_up_to
from prime_cut_prefix import prime_cut_prefix_failure


def support(value, primes):
    return {p for p in primes if value % p == 0}


def packet_room(counts, last, primes, slots):
    for p in primes - {2}:
        packet = (p - 1) // 2
        runs = [n for b, n in counts.items() if b % p == 0 and b != last]
        extendible = (counts[last] if last % p == 0 else 0) + slots
        runs.append(extendible)
        if sum(runs) < 2 * packet:
            return False
        if max(runs) < 2 * packet and sum(n // packet for n in runs) < 3:
            return False
    return True


def witness_bound(prefix, primes, L, E, slots, high):
    """Bound the next denominator when exactly three primes are available."""
    fours = prefix.count(4)
    fives = 5 in prefix
    if slots == 2 and ((fours >= 4 and L % 8) or (fives and L % 125)):
        return 0
    if slots == 3 and fours >= 4 and L % 8 and (L % 4 or E % 2 == 0):
        return 0
    choices = []
    for p in primes:
        seen_power, cofactor = 1, L
        while cofactor % p == 0:
            seen_power *= p
            cofactor //= p
        floor = 8 if p == 2 and fours >= 4 else 125 if p == 5 and fives else p
        options = []
        for growing in (False, True):
            if not growing and seen_power < floor:
                continue
            power = max(floor, seen_power*p) if growing else seen_power
            known = set() if growing else {i for i, b in enumerate(prefix) if b % power == 0}
            missing4 = fours >= 4 and all(prefix[i] % 4 for i in known)
            missing5 = fives and all(prefix[i] % 25 for i in known)
            need = max((3 if p in (2, 5) else 2) - len(known), int(missing4 or missing5))
            if p == 2:
                need += (3 - len(known) - need) % 4
            if p == 3:
                need += (2 - len(known) - need) % 3
            if need > slots:
                continue
            mass = 0
            if need:
                largest = max(power, lcm(power, 4) if missing4 else power,
                              lcm(power, 25) if missing5 else power)
                def cost(divisor):
                    b = ((prefix[-1]-1)//divisor + 1)*divisor
                    return (L*L + b*b - 1)//(b*b)
                mass = (need-1)*cost(power) + cost(largest)
            options.append((need, mass, known))
        choices.append(options)
    best = 0
    for a, b, c in product(*choices):
        total = a[0] + b[0] + c[0]
        if total > slots or a[2] & b[2] or (a[2] | b[2]) & c[2]:
            continue
        reserved = 1 + a[1] + b[1] + c[1]
        bound = high
        if E > reserved:
            bound = min(high, isqrt((slots - total)*L*L // (E - reserved)))
        best = max(best, bound)
    return best


def tail(prefix, L, E, primes, low, high):
    # With two suffix positions only 3 can be a new odd prime.
    allowed = primes | {3}
    for b, k, d, e in terminal(L, E, low, high, allowed, fixed_support=False):
        M = k*L
        bs = prefix + (b, M//d)
        ps = primes | support(k, allowed)
        counts = Counter(bs)
        if (counts[4] > 13
                or counts[5] == 20 or counts[5] > 22):
            continue
        if not packet_room(counts, bs[-1], ps, 0):
            continue
        row = [1] + [M//z for z in reversed(bs)]
        if any(len(support(x, ps)) < 2 for x in row[1:]):
            continue
        if all(sum(x % p != 0 for x in row[1:]) >= 2 for p in ps):
            yield row


def complete(prefix, rank, *, _handoff=False):
    """Yield all rows under this prefix, in exactly the rank-search degree convention."""
    forced = tuple(prefix)
    if (not isinstance(rank, int) or not 9 <= rank <= 34 or not forced or
        len(forced) > rank - 3 or any(type(b) is not int or b < 1 for b in forced) or
        forced != tuple(sorted(forced)) or forced[0] < 4):
        raise ValueError("expected a sorted positive square-search prefix of length 1 to rank-3")
    first = forced[0]
    if first*first > rank-2:
        return

    primes = primes_up_to(rank - 1)
    max_fours = 13  # Uniform multiplicity bound for the degree M/4.
    # visit and its smooth-number iterator can both be active down the
    # search tree. Allow their joint depth even for ranks above native widths.
    import sys
    required_depth = 3*rank + 128
    if required_depth > sys.getrecursionlimit():
        sys.setrecursionlimit(required_depth)

    def visit(bs, L, E):
        h = rank - 1 - len(bs)
        counts = Counter(bs)
        forbidden_fives = {13: 4, 12: 6, 10: 9}
        if counts[4] > max_fours or (counts[4] in forbidden_fives and
                             counts[5] >= forbidden_fives[counts[4]]):
            return
        if counts[5] > 22 or (counts[5] == 20 and bs[-1] > 5):
            return
        # A sum of two dyadic roots has mean square zero or at least two.
        # At twelve fours and five fives this forces the two eight twists
        # to sum to zero. A closed singleton nine then outweighs the tail.
        if (counts[4]==12 and counts[5]==5 and counts[8]==2
                and counts[9]==1 and bs[-1]>9):
            return
        # With two slots k=M/L is odd; with three, even k requires E odd.
        if counts[4] == 13 and counts[5] == 3 and 6 in counts:
            return
        if 36*E > (36*h + 1)*L*L:
            return
        ps = support(L, primes)
        deficit = 0
        scale = 256 if L < 2**48 else 1
        scaled_ll = scale*L*L
        last = bs[-1]
        packet_mass = min(scale, 1+(scaled_ll-1)//(36*last*last))
        for p in ps - {2}:
            step = (p-1)//2
            capacity = sum(n - n%step for b, n in counts.items() if b%p == 0)
            if bs[-1]%p == 0:
                capacity += counts[bs[-1]]%step
            need = max(0, 2*step - capacity)
            deficit += need
            if h >= 6 and need:
                multiple = ((last-1)//p+1)*p
                packet_mass += need*(1+(scaled_ll-1)//(multiple*multiple))
        if deficit > h:
            return
        future = ps | {2} | {p for p in primes if p <= h + 1 - deficit}
        if len(future) < 3 or not packet_room(counts, bs[-1], ps, h):
            return
        # The packet union bound (Hall) reserves distinct suffix labels
        # for each outstanding prime. At least two odd primes occur.
        missing = max(0, 2-len(ps-{2}))
        for p in sorted(future-ps-{2})[:missing]:
            need = p-1
            deficit += need
            if h >= 6:
                floor = max(last, p)
                packet_mass += need*(1+(scaled_ll-1)//(floor*floor))
        if deficit > h:
            return
        upper = isqrt((36*h + 1)*L*L // (36*E))
        if h >= 6 and scale*E > packet_mass:
            upper = min(upper, isqrt((h-deficit)*scaled_ll // (scale*E-packet_mass)))
        high = witness_bound(bs, future, L, E, h, upper) if len(future) == 3 else upper
        if high < bs[-1]:
            return
        low = max(bs[-1], isqrt(L*L // E) + 1)
        # Future primes need not all occur. Every tested cut contains a
        # seen prime on each side, so it is nontrivial in every completion.
        if prime_cut_prefix_failure(bs, h, tuple(sorted(future)), low,
                Fraction(E, L*L), required_primes=ps,
                vacuum_upper=Fraction(1, max(L*L, 36*low*low))) is not None:
            return
        if h == 2:
            yield from tail(bs, L, E, ps, low, high)
            return
        if low > high:
            return
        if len(bs) < len(forced):
            b = forced[len(bs)]
            if not low <= b <= high or not smooth_over(b, future):
                return
            choices = (b,)
        else:
            choices = smooth_interval(low, high, future)
        for b in choices:
            g = gcd(b, L)
            k = b//g
            yield from visit(bs + (b,), k*L, k*k*E - (L//g)**2)

    # The native driver has already checked every ancestor up to its
    # first LCM width crossing. Replaying all those height/packet tests
    # for nearby siblings was expensive in the rank-32 laptop run.
    # Reconstruct the two exact integers at that crossing; if the
    # original request fixed a longer prefix, visit still validates it.
    # A narrow-width test build may hand off earlier, in which case
    # normal independent traversal is used from the root.
    if _handoff:
        L = first
        E = first*first-1
        for i, b in enumerate(forced[1:], start=2):
            g = gcd(L, b)
            k = b//g
            E = E*k*k - (L//g)**2
            L *= k
            if L > (1 << 58):
                if E <= 0:
                    return
                yield from visit(forced[:i], L, E)
                return
    yield from visit((first,), first, first*first - 1)


def resume(prefix, rank):
    """Complete a certified native width handoff (unbounded Python integers)."""
    yield from complete(prefix, rank, _handoff=True)


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", help="comma-separated denominators, sorted")
    parser.add_argument("--rank", type=int, required=True)
    args = parser.parse_args()
    for row in complete([int(b) for b in args.prefix.split(",")], args.rank):
        print(json.dumps(row, separators=(",", ":")))
