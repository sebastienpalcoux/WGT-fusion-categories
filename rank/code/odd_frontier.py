"""Fresh exact search for odd-dimensional simple perfect integral modular types.

Usage: python odd_frontier.py 57 odd57.json
Add --audit-prime-codegrees to retain profiles with prime minimum codegree.
All denominators and all bounds use arbitrary-precision integers.
"""
import argparse
import json
from bisect import bisect_left
from itertools import combinations
from pathlib import Path
import sys
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm, prod


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


def enumerate_odd(rank_limit, audit_prime_codegrees=False, omit_square_class_one=False, selected_support=None):
    primes = primes_through(2 * rank_limit + 1)
    half_orbit = [(p - 1) // 2 for p in primes]
    indices = range(len(primes))
    simple_block_primes = [p for p in primes if not any((p-1)%(t*t)==0 for t in primes)]

    @lru_cache(maxsize=32768)
    def block_multiple(q,b):
        result=1
        for p in simple_block_primes:
            if q%p==0 and b%p:
                result *= p
        # Complex conjugation fixes no nonunit. In Gal(Q(zeta_f)/Q)
        # it is a 2^(s-1)-st power, s=min v2(p-1) for p dividing f.
        two_parts = [(p-1)&(1-p) for p in primes if q*b%p==0]
        if two_parts:
            result *= min(two_parts)//2
        return result

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
    @lru_cache(maxsize=4096)
    def divisors(value):
        ds = [1]
        n = value
        for p in primes:
            if n % p:
                continue
            old = ds[:]
            power = 1
            while n % p == 0:
                n //= p
                power *= p
                ds.extend(d*power for d in old)
            if n == 1:
                break
        assert n == 1
        return [d for d in ds if d >= 15 and prime_factors(d).bit_count() >= 2]

    @lru_cache(maxsize=4096)
    def fiber_size(q, omit=0):
        other = [p for p in primes if q % p == 0 and p != omit]
        if not other:
            return 1
        result = (other[0]+3)//2
        for p in other[1:]:
            result *= (p+1)//2
        return result

    @lru_cache(maxsize=4096)
    def packet_need(q,p):
        return (1 if q%p==0 else 2)*fiber_size(q,p)

    target_mask = 0
    target_primes = ()
    heights = ()
    smooth_small = ()
    scale = 1 << 64
    @lru_cache(maxsize=65536)
    def upper_square(n):
        return (scale+n*n-1)//(n*n)

    @lru_cache(maxsize=20000)
    def actual_height(L,p,hp):
        while L%(hp*p)==0:
            hp *= p
        return hp

    @lru_cache(maxsize=8192)
    def new_prime_mask(q,slots):
        return mask(i for i in indices if q%primes[i] == 0 or
                    packet_need(q,primes[i])*half_orbit[i] <= slots)

    def nonsquare_prefix(q,qmask,L,E,lower,prefix,seen,slots,cap):
        possible = seen | new_prime_mask(q,slots)
        # A prospective dimension prime cannot consume more normalized
        # mass than its required eligible denominators permit.
        candidates = possible & ~seen
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            i = bit.bit_length()-1
            p = primes[i]
            n = packet_need(q,p)*half_orbit[i]
            b0 = lower | 1
            bp = p*(((lower+p-1)//p)|1)
            bound = 2*((slots-n)*upper_square(b0)+n*upper_square(bp)) + upper_square(15*bp)
            if E*scale > L*L*bound:
                possible ^= bit
        if possible.bit_count()<3:
            return False
        forced=possible if possible.bit_count()==3 else seen
        floors=[lower]*slots
        witnesses=[]
        needs=[]
        position=0
        while forced:
            bit=forced & -forced
            forced ^= bit
            i=bit.bit_length()-1
            p=primes[i]
            if not qmask & bit:
                needed=packet_need(q,p)*half_orbit[i]-sum(b%p==0 for b in prefix)
                if needed>slots:return False
                if needed>0:
                    needs.append((p,needed))
                    quot=(lower+p-1)//p
                    bound=(quot | 1)*p
                    for j in range(slots-needed,slots):
                        floors[j]=max(floors[j],bound)
            if possible.bit_count()==3:
                need_exponent=(6,5,4)[position]
                h=(need_exponent-int(q%p==0)+1)//2
                height=actual_height(L,p,p**h)
                missing=(4 if p==3 else 2)-sum(b%height==0 for b in prefix)
                witnesses.extend([height]*max(0,missing))
            position+=1
        if cap>=75:
            # Mark n and m necessary future denominators. The intersection
            # size j ranges over an integer interval and its mass bound is affine.
            def ceiling_weight(t):
                return upper_square(t*(((lower+t-1)//t)|1))
            for (p,n),(r,m) in combinations(needs,2):
                lowj=max(0,n+m-slots)
                highj=min(n,m)
                def paired(j):
                    return ((slots-n-m+j)*ceiling_weight(1)+(n-j)*ceiling_weight(p)+
                            (m-j)*ceiling_weight(r)+j*ceiling_weight(p*r))
                mass=2*max(paired(lowj),paired(highj))+upper_square(15*(lower|1))
                if E*scale > L*L*mass:
                    return False
        if len(witnesses)>slots:return False
        for j,z in enumerate(sorted(witnesses,reverse=True)):
            at=slots-1-j
            floors[at]=max(floors[at],z)
        mass=2*sum(upper_square(t) for t in floors)+upper_square(15*floors[-1])
        return E*scale<=L*L*mass

    def support_mass_possible(q,L,E,lower,slots,prefix,one,two,run,active):
        # Sorted lower bounds for the h still-unseen pair denominators.
        floors = [lower]*slots
        free_floors = []
        for (p,hp,i) in heights:
            hp = actual_height(L,p,hp)
            h = half_orbit[i]
            bit = 1 << i
            packets = 2 - bool(one & bit) - bool(two & bit)
            count = max(0, packets*h - (run%h if active & bit else 0))
            free = (4 if p==3 else 2) - sum(b%hp==0 for b in prefix)
            if free>0:
                free_floors.extend([hp]*free)
            for mult,n in ((p,count),(hp,free)):
                if n<=0: continue
                if n>slots: return False
                quotient = max(1,(lower+mult-1)//mult)
                if quotient<=1000:
                    quotient=smooth_small[bisect_left(smooth_small,quotient)]
                bound=mult*quotient
                for j in range(slots-n,slots):
                    if bound>floors[j]:floors[j]=bound
        # A nonunit degree retains at least two primes, so can be free
        # of at most k-2 distinct support primes. Apply at each threshold.
        k = len(heights)-2
        free_floors.sort(reverse=True)
        if len(free_floors)>k*slots:
            return False
        for j in range(0,len(free_floors),k):
            ix=slots-1-j//k
            floors[ix]=max(floors[ix],free_floors[j])
        # Upper round each summand before comparing by cross-multiplication.
        mass = 2*sum(upper_square(t) for t in floors) + upper_square(15*floors[-1])
        return E*scale <= L*L*mass

    visited = 0
    orbit_pruned = 0
    survivors = []

    def record(q, M, denominators):
        if target_mask and prime_factors(q*M) != target_mask:
            return
        if any(count % block_multiple(q, b) for b, count in Counter(denominators).items()):
            return
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
            h = (p-1)//2
            if sum(n//h for n in active) < packet_need(q,p):
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
        if prefix and q*last*last in (25,27) and 2*run >= q*last*last-1:
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
        if prefix and q*last*last in (25,27) and 2*(run+1)>=q*last*last-1 and cap>=75:
            lower=max(lower,last+2)
        if q>1 and cap>=67 and not nonsquare_prefix(q,qmask,L,E,lower,prefix,seen,slots,cap):
            return
        if target_mask and not support_mass_possible(q,L,E,lower,slots,prefix,one,two,run,active):
            return
        open_multiple = block_multiple(q,last) if prefix else 1
        debt = (-run) % open_multiple
        if debt:
            # Every missing copy of the open run must precede any larger b.
            if (q*last*last in (25,27) and 2*(run+debt)>=q*last*last-1) or debt > slots or E*last*last <= 2*debt*L*L:
                return
            forced |= 1
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
        if slots == 2:
            if (final_rank+2) % 8 != q % 8:
                return
            # Both final degrees divide L. Near a narrow denominator interval
            # the direct Egyptian search below is cheaper than factorization.
            upper2 = isqrt(901*L*L//(225*E))
            if upper2 < lower:
                return
            if upper2-lower > 40:
                eligible = divisors(L)
                buckets = {}
                for e in eligible:
                    buckets.setdefault(2*e*e % E,[]).append(e)
                for d in eligible:
                    for e in buckets.get((-1-2*d*d)%E, ()):
                        if e>d:
                            continue
                        k2 = (1+2*d*d+2*e*e)//E
                        k = isqrt(k2)
                        if k*k != k2 or gcd(k,d*e) != 1:
                            continue
                        M = k*L
                        b,c = M//d,M//e
                        if b < lower or (forced and b != last):
                            continue
                        record(q,M,prefix+(b,c))
                return
        upper = isqrt((450*slots+1)*L*L//(225*E))
        if lower > upper:
            return
        if forced:
            choices = (last,) if lower <= last <= upper else ()
        else:
            allowed = [p for p in (target_primes or primes) if p <= slots+1 or q*L % p == 0]
            choices = odd_smooth_between(lower, upper, allowed)
        for b in choices:
            if not prefix and q*b*b == 9:
                continue
            if not prefix or b != last:
                required = block_multiple(q,b)
                if required > slots or E*b*b <= 2*required*L*L:
                    continue
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
        if (q == 1 and omit_square_class_one) or (q != 1 and selected_support):
            continue
        if any(q % (p*p) == 0 for p in primes):
            continue
        lower = isqrt(5//q) + 1  # all ordinary codegrees are at least six
        if not audit_prime_codegrees and q in primes:
            lower = max(lower, 3)  # a denominator one would give prime codegree q
        cap = rank_limit - (rank_limit-q) % 8
        qmask = prime_factors(q)
        if any(fiber_size(q,p)*(1+(p-1)) > cap for p in primes if q%p==0):
            continue
        if q == 1 and cap >= 57:
            # The rank envelope was separately established in this manuscript.
            # Locate its module both beside this script and in the final package.
            for base in (Path(__file__).resolve().parent,
                         Path(__file__).resolve().parents[1]/'code',
                         Path(__file__).resolve().parents[1]/'package'/'code'):
                if (base/'rank_family.py').exists():
                    sys.path.insert(0,str(base)); break
            from rank_family import root_possible
            pool = [p for p in primes if 2*(p-1) <= cap-1]
            for count in range(3,len(pool)+1):
                for support in combinations(pool,count):
                    if selected_support and support != selected_support:
                        continue
                    exps=(3,3,2) if count==3 else (1,)*count
                    for first in odd_smooth_between(lower,isqrt(cap-1),support):
                        if first*first<8 or first==3:continue
                        if not root_possible(support,exps,q,cap,first):
                            continue
                        target_primes=support
                        target_mask=prime_factors(prod(support))
                        heights=tuple((p,p**h,primes.index(p)) for p,h in zip(support,exps))
                        sm=sorted(odd_smooth_between(1,3000,support))
                        smooth_small=sm+[3001]
                        bit=prime_factors(first)
                        earned=bit & earns[1]
                        visit(q,0,cap,first,first*first-2,first,(first,),
                              target_mask,earned,0,1,bit)
            target_primes=(); target_mask=0
        else:
            visit(q, qmask, cap, 1, q, lower, (), qmask, 0, 0, 0, 0)
    result = {"max_rank":rank_limit, "retain_prime_codegrees":audit_prime_codegrees,
              "visited":visited, "orbit_pruned":orbit_pruned,
              "row_feasible":sorted(survivors)}
    if omit_square_class_one:
        result["omitted_square_class"] = 1
    if selected_support:
        result["square_class"] = 1
        result["exact_prime_support"] = selected_support
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rank", type=int)
    parser.add_argument("output")
    parser.add_argument("--audit-prime-codegrees", action="store_true")
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--omit-square-class-one", action="store_true")
    group.add_argument("--square-support", metavar="P1,P2,P3,...",
                       help="Check only q=1 on this exact prime support; minimum cutoff 57.")
    args = parser.parse_args()
    support = None if not args.square_support else tuple(map(int,args.square_support.split(',')))
    if support and (args.rank<57 or len(support)<3 or tuple(sorted(set(support)))!=support
                    or not set(support)<=set(primes_through(2*args.rank+1))):
        parser.error("supply an increasing list of at least three odd primes at cutoff at least 57")
    data = enumerate_odd(args.rank, args.audit_prime_codegrees, args.omit_square_class_one, support)
    with open(args.output, "w", encoding="utf-8") as output:
        json.dump(data, output, separators=(",", ":"))
        output.write("\n")
