"""Exact continuation of a square-dimension denominator prefix, using Python integers."""
from bisect import bisect_left
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def support(value):
    return {p for p in PRIMES if value % p == 0}


def smooth_over(value, primes):
    for p in primes:
        while value % p == 0:
            value //= p
    return value == 1


def smooth_interval(low, high, primes):
    if high - low <= 35:
        for value in range(low, high + 1):
            if smooth_over(value, primes):
                yield value
        return
    primes = sorted(primes, reverse=True)

    def walk(i, value):
        if i == len(primes):
            if value >= low:
                yield value
            return
        p = primes[i]
        while value <= high:
            yield from walk(i + 1, value)
            value *= p

    yield from walk(0, 1)


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


@lru_cache(maxsize=32)
def nonunit_divisors(value):
    # Track the number of distinct prime factors while generating each divisor.
    divisors = [(1, 0)]
    for p in PRIMES:
        if value % p:
            continue
        old = tuple(divisors)
        power = 1
        while value % p == 0:
            value //= p
            power *= p
            divisors.extend((d * power, count + 1) for d, count in old)
    if value != 1:
        raise ValueError("prefix lcm has a prime factor outside the search")
    return sorted(d for d, count in divisors if count >= 2)


def tail(prefix, L, E, primes, low, high):
    """Solve E k²=1+d²+e² for d≤e dividing L and M=kL."""
    allowed = primes | {3}
    divisors = nonunit_divisors(L)
    if not divisors or low > high:
        return

    def candidate(k, d, e):
        if k < 1 or gcd(k, 10*d*e) != 1 or not smooth_over(k, allowed):
            return
        M = k*L
        if M//e < prefix[-1]:
            return
        bs = prefix + (M//e, M//d)
        ps = primes | support(k)
        if not packet_room(Counter(bs), bs[-1], ps, 0):
            return
        row = [1] + [M//b for b in reversed(bs)]
        if any(len(support(x)) < 2 for x in row[1:]):
            return
        if any(sum(x % p != 0 for x in row[1:]) < 2 for p in ps):
            return
        yield row

    if high - low <= max(32, 2*len(divisors)):
        for b in range(low, high + 1):
            if not smooth_over(b, allowed):
                continue
            g = gcd(b, L)
            k, e = b//g, L//g
            if gcd(k, 10) != 1 or len(support(e)) < 2:
                continue
            d2 = E*k*k - e*e - 1
            if d2 <= 0:
                continue
            d = isqrt(d2)
            if d*d == d2 and d <= e and L % d == 0 and len(support(d)) >= 2:
                yield from candidate(k, d, e)
        return

    buckets = defaultdict(list)
    for e in divisors:
        buckets[e*e % E].append(e)
    for d in divisors:
        bucket = buckets.get(-(1 + d*d) % E, ())
        for i in range(bisect_left(bucket, d), len(bucket)):
            e = bucket[i]
            k2 = (1 + d*d + e*e) // E
            k = isqrt(k2)
            if k*k == k2:
                yield from candidate(k, d, e)


def complete(prefix, rank=31):
    """Yield all rows under this prefix, in exactly the rank-search degree convention."""
    forced = tuple(prefix)
    if (not isinstance(rank, int) or not 17 <= rank <= 31 or not forced or
        len(forced) > rank - 3 or any(type(b) is not int or b < 1 for b in forced) or
        forced != tuple(sorted(forced)) or forced[0] not in (4, 5)):
        raise ValueError("expected a sorted positive square-search prefix of length 1 to rank-3")
    first = forced[0]
    if first*first >= rank:
        return

    def visit(bs, L, E):
        h = rank - 1 - len(bs)
        counts = Counter(bs)
        forbidden_fives = {13: 4, 12: 6, 10: 9}
        if counts[4] > 13 or (counts[4] in forbidden_fives and
                             counts[5] >= forbidden_fives[counts[4]]):
            return
        if counts[4] == 13 and counts[5] == 3 and 6 in counts:
            return
        if 36*E > (36*h + 1)*L*L:
            return
        ps = support(L)
        deficit = 0
        for p in ps - {2}:
            step = (p-1)//2
            capacity = sum(n - n%step for b, n in counts.items() if b%p == 0)
            if bs[-1]%p == 0:
                capacity += counts[bs[-1]]%step
            deficit += max(0, 2*step - capacity)
        if deficit > h:
            return
        future = ps | {2} | {p for p in PRIMES if p <= h + 1 - deficit}
        if len(future) < 3 or not packet_room(counts, bs[-1], ps, h):
            return
        upper = isqrt((36*h + 1)*L*L // (36*E))
        high = witness_bound(bs, future, L, E, h, upper) if len(future) == 3 else upper
        if high < bs[-1]:
            return
        low = max(bs[-1], isqrt(L*L // E) + 1)
        if h == 2:
            yield from tail(bs, L, E, ps, low, upper)
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

    yield from visit((first,), first, first*first - 1)


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", help="comma-separated denominators, sorted")
    parser.add_argument("--rank", type=int, default=31)
    args = parser.parse_args()
    for row in complete([int(b) for b in args.prefix.split(",")], args.rank):
        print(json.dumps(row, separators=(",", ":")))
