"""Fresh exact search for odd-dimensional simple perfect integral modular types.

Usage: python odd_rank.py 39 odd39.json
Add --audit-prime-codegrees to retain profiles with prime minimum codegree.
All denominators and all bounds use arbitrary-precision integers.
"""
import argparse
import json
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm


def primes_through(limit):
    return [p for p in range(3, limit + 1, 2)
            if all(p % j for j in range(3, isqrt(p) + 1, 2))]


def odd_smooth_between(lower, upper, primes):
    """Enumerate each odd integer in the interval supported on these primes."""
    if upper - lower < 96:
        for b in range(lower | 1, upper + 1, 2):
            rest = b
            for p in primes:
                while rest % p == 0:
                    rest //= p
            if rest == 1:
                yield b
        return

    def products(i, value):
        if i == len(primes):
            if value >= lower:
                yield value
            return
        while value <= upper:
            yield from products(i + 1, value)
            value *= primes[i]

    yield from products(0, 1)


def enumerate_odd(rank_limit, audit_prime_codegrees=False):
    primes = primes_through(2 * rank_limit + 1)
    half_orbit = [(p - 1) // 2 for p in primes]
    indices = range(len(primes))

    def mask(indices_to_include):
        return sum(1 << i for i in indices_to_include)

    @lru_cache(maxsize=65536)
    def prime_factors(n):
        return mask(i for i in indices if n % primes[i] == 0)

    # Bit i refers to primes[i]. A run earns a new same-degree packet
    # each time its length reaches a multiple of h=(p-1)/2.
    earns = [mask(i for i in indices if n % half_orbit[i] == 0)
             for n in range(rank_limit + 1)]
    room_one = [mask(i for i in indices if half_orbit[i] <= n)
                for n in range(rank_limit + 1)]
    room_two = [mask(i for i in indices if 2 * half_orbit[i] <= n)
                for n in range(rank_limit + 1)]
    visited = 0
    orbit_pruned = 0
    survivors = []

    def record(q, M, denominators):
        degrees = [M // b for b in denominators]  # once per dual pair
        if any(prime_factors(d).bit_count() < 2 for d in degrees):
            return
        if any(gcd(*(degrees[:i] + degrees[i+1:])) != 1 for i in range(len(degrees))):
            return
        counts = Counter(denominators)
        # At termination use the one-packet and total-label conditions;
        # these are weaker than the two-packet pruning used inside the tree.
        for p in primes:
            if q*M % p:
                continue
            active = [n for b, n in counts.items() if q*b % p == 0]
            if not active or max(active) < (p-1)//2:
                return
            if q % p and sum(active) < p-1:
                return
        dimension = q*M*M
        for x in set(degrees):
            available = [y for y in set(degrees)
                         if gcd(x,y) > 1 and lcm(x,y)**2 + x*x <= dimension]
            if gcd(*available, 0) != 1:
                return
        doubled = [d for d in degrees for _ in range(2)]
        assert 1 + sum(d*d for d in doubled) == dimension
        survivors.append([q, M, doubled])

    def visit(q, qmask, cap, L, E, last, prefix, seen, one, two, run, active):
        nonlocal visited, orbit_pruned
        visited += 1
        slots = (cap-1)//2 - len(prefix)
        if slots <= 0:
            return
        # q-primes need one packet. Other primes in M need two. If the
        # missing packets do not fit in fresh slots, the open run must continue.
        need_two = seen & ~qmask & ~one
        need_one = (qmask & ~one) | (seen & ~qmask & one & ~two)
        forced = (need_two & ~room_two[slots]) | (need_one & ~room_one[slots])
        pending = forced
        while pending:
            bit = pending & -pending
            h = half_orbit[bit.bit_length()-1]
            number_needed = 2 if need_two & bit else 1
            if not active & bit or slots + run % h < number_needed*h:
                orbit_pruned += 1
                return
            pending ^= bit

        lower = max(last, isqrt(2*L*L//E) + 1)
        if forced and lower > last:
            return
        # The final pair cannot introduce a new maximal prime height.
        final_rank = 2*len(prefix) + 3
        if final_rank % 8 == q % 8:
            d = isqrt((E-1)//2)
            if d >= 15 and 1+2*d*d == E and L % d == 0 and L//d >= last:
                record(q, L, prefix+(L//d,))
        if slots == 1:
            return
        upper = isqrt((450*slots+1)*L*L//(225*E))
        if lower > upper:
            return
        if forced:
            choices = (last,) if lower <= last <= upper else ()
        else:
            allowed = [p for p in primes if p <= slots+1 or q*L % p == 0]
            choices = odd_smooth_between(lower, upper, allowed)
        for b in choices:
            g = gcd(L,b)
            k = b//g
            new_lcm = L*k
            new_residual = E*k*k - 2*(L//g)**2
            if slots == 2:
                if (final_rank+2) % 8 == q % 8:
                    d = isqrt((new_residual-1)//2)
                    if (d >= 15 and 1+2*d*d == new_residual and
                            new_lcm % d == 0 and new_lcm//d >= b):
                        record(q, new_lcm, prefix+(b,new_lcm//d))
                continue
            new_run = run+1 if prefix and b == last else 1
            new_active = qmask | prime_factors(b)
            earned = new_active & earns[new_run]
            visit(q, qmask, cap, new_lcm, new_residual, b, prefix+(b,),
                  seen | new_active, one | earned, two | (one & earned),
                  new_run, new_active)

    for q in range(1, rank_limit, 2):
        if any(q % (p*p) == 0 for p in primes):
            continue
        lower = isqrt(5//q) + 1  # all ordinary codegrees are at least six
        if not audit_prime_codegrees and q in primes:
            lower = max(lower, 3)  # a denominator one would give prime codegree q
        cap = rank_limit - (rank_limit-q) % 8
        qmask = prime_factors(q)
        visit(q, qmask, cap, 1, q, lower, (), qmask, 0, 0, 0, 0)
    return {"max_rank":rank_limit, "retain_prime_codegrees":audit_prime_codegrees,
            "visited":visited, "orbit_pruned":orbit_pruned,
            "row_feasible":sorted(survivors)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rank", type=int)
    parser.add_argument("output")
    parser.add_argument("--audit-prime-codegrees", action="store_true")
    args = parser.parse_args()
    data = enumerate_odd(args.rank, args.audit_prime_codegrees)
    with open(args.output, "w", encoding="utf-8") as output:
        json.dump(data, output, separators=(",", ":"))
        output.write("\n")
