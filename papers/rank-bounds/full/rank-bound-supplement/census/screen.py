"""Terminal necessary conditions for the based-simple rank search."""
from collections import Counter, defaultdict
from functools import lru_cache
from heapq import heappop, heappush
from itertools import product
from math import gcd, isqrt, prod

from engine.defect import defect_tiles
from multiplication import multiplication_obstruction
from trace import trace_obstruction


OBSTRUCTIONS = frozenset({"prime_power", "defect_one", "degree",
                          "crossed", "activity", "balance", "orbit", "dyadic",
                          "multiplication"})


def validate_profile(degrees, q, B, rank):
    """Validate an emitted arithmetic profile. Bad output is an error, never a rejection."""
    if (not isinstance(degrees, (list, tuple)) or len(degrees) != rank
            or any(type(d) is not int for d in degrees) or degrees[0] != 1
            or degrees[1] < 6 or any(a > b for a, b in zip(degrees, degrees[1:]))
            or gcd(*degrees[1:]) != 1):
        raise ValueError("invalid sorted simple-degree profile")
    M = B * degrees[-1]
    D = q * M * M
    if (D != sum(d*d for d in degrees) or any(M % d for d in degrees)):
        raise ValueError("profile violates dimension or denominator identity")
    return D


@lru_cache(maxsize=8192)
def factor(n):
    factors = []
    p = 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            factors.append((p, e))
        p = 3 if p == 2 else p+2
    if n > 1:
        factors.append((n, 1))
    return factors


@lru_cache(maxsize=32768)
def degree_masks(d, D):
    """Prime support of d, and primes of defect at least two at this degree."""
    support = defect = 0
    f = D//(d*d)
    for i, (p, _) in enumerate(factor(D)):
        if d % p == 0:
            support |= 1 << i
        if f % (p*p) == 0:
            defect |= 1 << i
    return support, defect


def principal_degrees_fail(values, D, factors):
    all_primes = (1 << len(factors))-1
    classes = {}
    for d in values:
        support, defect = degree_masks(d, D)
        classes[support] = classes.get(support, 0) | defect
    common = [all_primes if e >= 2 else 0 for p, e in factors]
    for support, admissible in classes.items():
        # An opposite, coprime degree must be p-divisible for a principal row.
        for other in classes:
            if not support & other:
                admissible &= other
        while admissible:
            bit = admissible & -admissible
            common[bit.bit_length()-1] &= support
            admissible ^= bit
    return any(common)


def large_coin_oracle(values):
    """Exact membership using bounded memory, even for very large moduli."""
    coins = [values[0]]
    divisors = [values[0]]
    cutoffs = [0]
    for a in values[1:]:
        old_g, old_c = divisors[-1], cutoffs[-1]
        g = gcd(old_g, a)
        if g == old_g and a >= old_c:
            continue  # Already represented by the preceding coins.
        # Every multiple of g from c onward is represented: subtract k*a
        # with 0 <= k < old_g/g to land on a multiple of old_g >= old_c.
        c = old_c + (old_g//g-1)*a-old_g+g
        coins.append(a)
        divisors.append(g)
        cutoffs.append(c)

    periods = [1]
    inverses = [0]
    for i in range(1, len(coins)):
        period = divisors[i-1]//divisors[i]
        periods.append(period)
        inverses.append(pow(coins[i]//divisors[i], -1, period))

    @lru_cache(maxsize=4096)
    def represented(i, n):
        if n < 0 or n % divisors[i]:
            return False
        if n >= cutoffs[i]:
            return True
        p = periods[i]
        # These and only these multiplicities leave a multiple of the
        # preceding gcd. The recursive depth is bounded by the rank.
        k = (n//divisors[i]*inverses[i]) % p
        while k*coins[i] <= n:
            if represented(i-1, n-k*coins[i]):
                return True
            k += p
        return False

    last = len(coins)-1
    return lambda n: represented(last, n), cutoffs[last], represented.cache_clear


def coin_oracle(values):
    """Small moduli use the least representable element in each residue."""
    m = values[0]
    if m > 2048:
        return large_coin_oracle(values)
    # In one residue only the lightest coin matters: all other differences
    # are nonnegative multiples of m. The input degrees are increasing.
    weights = {}
    for a in values[1:]:
        residue = a % m
        if residue and residue not in weights:
            weights[residue] = a
    dist = [None]*m
    dist[0] = 0
    queue = [(0, 0)]
    while queue:
        n, i = heappop(queue)
        if dist[i] != n:
            continue
        for offset, a in weights.items():
            j = i+offset
            if j >= m:
                j -= m
            if dist[j] is None or n+a < dist[j]:
                dist[j] = n+a
                heappush(queue, (n+a, j))
    return lambda n: n >= dist[n % m], max(dist)-m+1, None


def degrees_fail(values, D, factors):
    m = values[0]
    if any(m*a < b for a, b in zip(values, values[1:])) or gcd(*values) != 1:
        return True
    if principal_degrees_fail(values, D, factors):
        return True
    possible, conductor, cleanup = coin_oracle(values)
    try:
        for a in values:
            if a*a-1 >= conductor:
                break
            if not possible(a*a-1):
                return True
        # Once every consecutive pair is joined, the degree graph is
        # connected. The least residual for the pair (u,v), u<v, is a*u-v.
        threshold = max(((v+conductor+u-1)//u for u, v in zip(values, values[1:])), default=0)
        for a in values:
            if a >= threshold:
                break
            reached = {a}
            queue = [a]
            while queue:
                u = queue.pop()
                for v in values:
                    if v in reached:
                        continue
                    if not (possible(a*u-v) or a == u and possible(a*u-v-1)):
                        continue
                    if not (possible(a*v-u) or a == v and possible(a*v-u-1)):
                        continue
                    reached.add(v)
                    queue.append(v)
            if len(reached) != len(values):
                return True
        return False
    finally:
        if cleanup is not None:
            cleanup()


@lru_cache(maxsize=2048)
def balanced_capacity(groups, p):
    """Total mass in units of one packet, allowing higher primary levels."""
    weights = []
    for d, n, higher in groups:
        opts = []
        for k in range(-n, n+1):
            left = n-abs(k)
            unused = left % 2 and not (higher and left >= p)
            opts.append((k*d*d, (n-int(unused))*d*d))
        weights.append(opts)
    weights.sort(key=lambda opts: opts[-1][0], reverse=True)
    div = [0]*(len(weights)+1)
    bound = [0]*(len(weights)+1)
    for i in range(len(weights)-1, -1, -1):
        div[i] = gcd(div[i+1], *(w for w, _ in weights[i]))
        bound[i] = bound[i+1]+weights[i][-1][0]

    @lru_cache(maxsize=8192)
    def best(i, target):
        if i == len(weights):
            return 0 if target == 0 else -1
        if abs(target) > bound[i] or target % div[i]:
            return -1
        maximum = -1
        for signed, active in weights[i]:
            tail = best(i+1, target-signed)
            if tail >= 0:
                maximum = max(maximum, tail+active)
        return maximum
    answer = best(0, 0)
    best.cache_clear()
    return answer


@lru_cache(maxsize=32)
def ramified_primes(q):
    return tuple(p for p,e in factor(q))


@lru_cache(maxsize=32)
def fiber_floor(q):
    """Lower bound for the size of a complementary odd-phase fiber."""
    factors = [(p+1)//2 for p in ramified_primes(q) if p != 2]
    base = prod(factors)
    return max([(2 if q % 2 == 0 else 1)*base] + [(a+1)*(base//a) for a in factors])


def tiles_in_fiber(degree_counts, q, M):
    """Defect one at each ramified prime is equivalent to p not dividing M/d."""
    counts = {M//d: n for d, n in degree_counts if n}
    return all(defect_tiles(counts, p) for p in ramified_primes(q) if p in (2,3,5))


def square_kernel(factors):
    return prod(p for p, exponent in factors if exponent % 2)


def prime_power_rows(f, factors):
    """Return {row degree: prime} using factors of the already factored dimension."""
    result = {}
    for x, codegree in f.items():
        for p, _ in factors:
            if codegree % p:
                continue
            other = codegree
            while other % p == 0:
                other //= p
            if other == 1:
                result[x] = p
            break
    return result


def prime_power_failure(mult, f, pure_rows, rank):
    for x, p in pure_rows.items():
        codegree, m = f[x], mult[x]
        if (m > codegree - (1 + isqrt(p - 1))
                or (codegree != p and m == codegree - p)):
            return True
        if p == 2 and x % 2:
            if codegree >= 16 and m > codegree - 8:
                return True
            if codegree == 8 and m > 2:
                return True
        # A prime-to-p degree meets every one of these principal rows.
        # Its entry is divisible by y/gcd(x,y), and has mean squared norm
        # at least that square. A nonunit coprime entry must vanish.
        used = 0
        for y, copies in mult.items():
            if y % p:
                g = gcd(x, y)
                if y != 1 and g == 1:
                    return True
                used += copies*(y//g)**2
        if used > codegree:
            return True
    return False


def phase_vectors_exist(pool, h, p, q, M, D, least, minimum, rank):
    """Two primitive phase fibers of equal mass, after ruling out higher orders.

    pool gives each eligible degree and its number of available h-packets.
    Each vector records one representative fiber; its h conjugates have the
    same degrees. Both vectors must fit simultaneously. All optional limits
    in the caller skip tests, and never reject a candidate.
    """
    weights = [d*d for d, _ in pool]
    limits = [n for _, n in pool]
    by_mass = defaultdict(list)
    max_size = (rank - minimum) // h - minimum
    for vector in product(*(range(n+1) for n in limits)):
        size = sum(vector)
        if not minimum <= size <= max_size:
            continue
        mass = sum(k*w for k, w in zip(vector, weights))
        # With no higher orders the total active mass is 2h times this mass.
        if 2*p*least*mass < D:
            continue
        if q % 2 == 0 and size == 2 and 2 not in vector:
            # Lemma on a dyadically invariant two-label phase fiber.
            continue
        if q > 1:
            occupied = [(d,k) for (d,_),k in zip(pool, vector) if k]
            if not tiles_in_fiber(occupied, q, M):
                continue
            g = gcd(*(d for d,_ in occupied))
            # For every odd ramified prime t, the common mass is divisible
            # by t*g². Since q is square-free these conditions combine.
            odd_ramified = q//2 if q % 2 == 0 else q
            if mass % (odd_ramified*g*g):
                continue
        matching = by_mass[mass]
        matching.append(vector)
        if any(all(a+b <= available for a, b, available in zip(vector, other, limits))
               for other in matching):
            return True
    return False


def dyadic_fail(degrees, D, mult):
    """The class traces and mean squared norms of a codegree-sixteen row."""
    x = degrees[-1]
    if D != 16*x*x:
        return False
    data = []
    for y, count in sorted(mult.items()):
        g = gcd(x, y)
        if y == 1 or g == 1 or D//(y*y) % 4:
            continue
        t = y//g
        # If t is even, x/g is odd. The two principal-row congruences
        # thus prescribe y mod 2 unless x is odd and y is even.
        parity = y % 2 if x % 2 == 0 or y % 2 else None
        options = []
        for B in range(15//(t*t)+1):
            if parity is not None and (B-count*parity) % 2:
                continue
            maximum = isqrt(B*count)
            for T in range(-maximum, maximum+1):
                if parity is not None and (T-count*parity) % 2:
                    continue
                j, r = divmod(abs(T), count)
                def energy(k):
                    return k*k+int(parity is not None and (k-parity) % 2)
                if B < (count-r)*energy(j)+r*energy(j+1):
                    continue
                scale = 256*t*y
                eta = 2 if parity == 0 else 1
                cap = 1+min(isqrt(scale**2*count*B), isqrt((scale*B)**2//eta))
                options.append((t*y*T, t*t*B, cap))
        data.append((t*y, options))
    div = [0]*(len(data)+1)
    for i in range(len(data)-1, -1, -1):
        div[i] = gcd(div[i+1], data[i][0])

    @lru_cache(maxsize=8192)
    def best(i, target, left):
        if i == len(data):
            return 0 if target == left == 0 else -1
        if target % div[i]:
            return -1
        maximum = -1
        for w, cost, cap in data[i][1]:
            if cost <= left:
                tail = best(i+1, target-w, left-cost)
                if tail >= 0:
                    maximum = max(maximum, tail+cap)
        return maximum
    failed = best(0, -1, 15) < 256*(4*x-1)
    best.cache_clear()
    return failed


def exclusion_reason(degrees, dimension=None):
    """Return the first obstruction, or None if every implemented test passes."""
    mult = Counter(degrees)
    values = sorted(d for d in mult if d != 1)
    D = dimension if dimension is not None else sum(d*d*n for d,n in mult.items())
    factors = factor(D)
    q = square_kernel(factors)
    M = isqrt(D//q)
    rank = len(degrees)
    if q > 1 and not tiles_in_fiber(mult.items(), q, M):
        return 'defect_one'
    f = {d: D//(d*d) for d in values}
    pure_rows = prime_power_rows(f, factors)
    if prime_power_failure(mult, f, pure_rows, rank):
        return 'prime_power'
    if degrees_fail(values, D, factors):
        return 'degree'

    pure_primes = set(pure_rows.values())
    fiber = fiber_floor(q)
    for p, e in factors:
        h = max(1, (p-1)//2)
        # A degree prime to t has purely t-primary twist whenever a t-power
        # codegree occurs. This also handles several such primes at once.
        groups = [(d, mult[d]//h, f[d] % (p*p) == 0)
                  for d in values if f[d] % p == 0 and mult[d] >= h
                  and all(d % t == 0 for t in pure_primes if t != p)]
        cap = sum(n*d*d for d, n, _ in groups)
        threshold = (p-1)*D
        scale = 2*p*values[0]*h
        if scale*cap < threshold:
            return 'activity'
        if p > 2 and e % 2 == 0:
            no_higher = (rank < fiber*(1+p*h)
                         or (q == 3 and rank < 4+3*p*h)
                         or not any(higher and n >= p for _, n, higher in groups))
            if no_higher:
                groups = [(d,n,False) for d,n,_ in groups]
            # Choosing zero imbalance independently at every degree gives
            # this feasible mass. Only smaller values need the optimization.
            neutral = cap - sum(d*d for d, n, higher in groups
                                if n % 2 and not (higher and n >= p))
            if scale*neutral < threshold and scale*balanced_capacity(tuple(groups), p) < threshold:
                return 'balance'
            if no_higher:
                pool = tuple((d,n) for d,n,_ in groups)
                # Keep the exponential optional phase test bounded. For the
                # mandatory q=3 test at p=5,7 the bound always fits at rank 32.
                if (prod(n+1 for d,n in pool) <= 32768
                        and not phase_vectors_exist(pool,h,p,q,M,D,values[0],fiber,rank)):
                    return 'crossed'
    if trace_obstruction(degrees, dimension=D, multiplicities=mult):
        return 'orbit'
    if dyadic_fail(degrees, D, mult):
        return 'dyadic'
    # The same Frobenius-reciprocity equations used in the archived
    # smallest-degree test apply for every simple multiplier. Go farther
    # only if the cheaper screens leave this list alive.
    for x in values:
        if x*values[0] > 3001:
            break  # All aggregate targets exceed the solver's bitset cutoff.
        if multiplication_obstruction(degrees, x):
            return 'multiplication'
    return None
