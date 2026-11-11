"""Terminal necessary conditions for a simple perfect integral modular category."""
from collections import Counter, defaultdict
from functools import lru_cache
from heapq import heappop, heappush
from itertools import product
from math import gcd, isqrt

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


def degrees_fail(values, codegrees, factors):
    for p, e in factors:
        if e >= 2:
            free = [d for d in values if d % p]
            pool = [y for y in values if codegrees[y] % (p*p) == 0
                    and all(gcd(x, y) > 1 for x in free)]
            if gcd(*pool) != 1:
                return True

    # The least representable integer in each class modulo the smallest coin.
    m = values[0]
    if gcd(*values) != 1:
        return True
    dist = [None]*m
    dist[0] = 0
    queue = [(0, 0)]
    while queue:
        n, i = heappop(queue)
        if dist[i] != n:
            continue
        for a in values[1:]:
            j = (i+a) % m
            if dist[j] is None or n+a < dist[j]:
                dist[j] = n+a
                heappush(queue, (n+a, j))

    def possible(n):
        return n >= 0 and dist[n % m] <= n

    if any(not possible(a*a-1) for a in values):
        return True
    for a in values:
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

    @lru_cache(None)
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
    return best(0, 0)


def possible_phase_fibers(degrees, p):
    """p=5 or 7 divides M, and D=3M². Compare the two square classes."""
    dimension = sum(d*d for d in degrees)
    step = (p-1)//2
    pool = [(d, m//step) for d, m in Counter(degrees).items()
            if d > 1 and m >= step and dimension % (p*d*d) == 0]
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

    @lru_cache(None)
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
    return best(0, -1, 15) < 256*(4*x-1)


def exclusion_reason(degrees):
    """Return the first obstruction, or None if every implemented test passes."""
    mult = Counter(degrees)
    values = sorted(set(degrees)-{1})
    D = sum(d*d for d in degrees)
    factors = factor(D)
    f = {d: D//(d*d) for d in values}
    if degrees_fail(values, f, factors):
        return 'degree'
    if isqrt(D)**2 != D:
        if len(degrees) <= 30 and 3*isqrt(D//3)**2 == D:
            if any(not possible_phase_fibers(degrees, p) for p, _ in factors if p in (5, 7)):
                return 'crossed'
        return None
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
            if scale*balanced_capacity(groups, p) < threshold:
                return 'balance'
    if trace_obstruction(degrees):
        return 'orbit'
    if dyadic_fail(degrees, D, mult):
        return 'dyadic'
    return None
