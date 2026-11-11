"""Terminal necessary conditions for a simple perfect integral modular category."""
from collections import Counter, defaultdict
from functools import lru_cache
from heapq import heappop, heappush
from itertools import product
from math import gcd, isqrt

from multiplication import multiplication_obstruction
from trace import trace_obstruction


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


@lru_cache(maxsize=2048)
def possible_phase_fibers(pool):
    """For D=3M² compare one fiber in each order-p square class."""
    by_mass = defaultdict(list)
    for packet in product(*(range(n+1) for d, n in pool)):
        if not any(packet):
            continue
        common = gcd(*(d for (d, _), k in zip(pool, packet) if k))
        mass = sum(k*d*d for (d, _), k in zip(pool, packet))
        if mass % (3*common*common):
            continue
        by_mass[mass].append(packet)
        for other in by_mass[mass]:
            if all(a+b <= n for a, b, (_, n) in zip(packet, other, pool)):
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
        parity = 0 if x % 2 == 0 and (x//g) % 2 else None
        if 2*mult[x] > 16 and t % 2:
            parity = y % 2
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
    if degrees_fail(values, D, factors):
        return 'degree'
    f = {d: D//(d*d) for d in values}
    if 3*isqrt(D//3)**2 == D:
        for p, _ in factors:
            # At order p² the Galois orbit has p(p-1)/2 phases.
            # Each phase, including the trivial phase, has at least 3 labels.
            if (p in (5, 7) and len(degrees) < 3*(1+p*(p-1)//2)
                    and not possible_phase_fibers(tuple((d, mult[d]//((p-1)//2)) for d in values
                        if mult[d] >= (p-1)//2 and f[d] % p == 0))):
                return 'crossed'
    for p, e in factors:
        h = max(1, (p-1)//2)
        groups = [(d, mult[d]//h, f[d] % (p*p) == 0)
                  for d in values if f[d] % p == 0 and mult[d] >= h]
        cap = sum(n*d*d for d, n, _ in groups)
        threshold = (p-1)*D
        scale = 2*p*values[0]*h
        if scale*cap < threshold:
            return 'activity'
        if p > 2 and e % 2 == 0:
            # Choosing zero imbalance independently at every degree gives
            # this feasible mass. Only smaller values need the optimization.
            neutral = cap - sum(d*d for d, n, higher in groups
                                if n % 2 and not (higher and n >= p))
            if scale*neutral < threshold and scale*balanced_capacity(tuple(groups), p) < threshold:
                return 'balance'
    if trace_obstruction(degrees, dimension=D, multiplicities=mult):
        return 'orbit'
    if dyadic_fail(degrees, D, mult):
        return 'dyadic'
    if multiplication_obstruction(degrees):
        return 'multiplication'
    return None
