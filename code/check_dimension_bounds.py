from heapq import heappop, heappush
from math import prod

def covered(primes, exponents):
    e = sorted(exponents)
    if len(e) == 3:
        a, b, c = exponents
        return (a <= 5 or b <= 4 or c <= 3 or
                all(x <= y for x, y in zip(e, (4, 5, 7))))
    if len(e) == 4:
        d, c, b, a = e
        dyadic = primes[0] == 2 and exponents[0] >= 4
        return b <= 3 or (d == 2 and c <= 3 and
               (a//2 + b//2 <= 5 or (dyadic and a <= 9 and b <= 5)))
    if len(e) == 5:
        return e[-1] <= 3 or (e[-1] == 4 and e[-2] <= 3 and e[0] == 2)
    return False

def first(primes, center):
    initial = (1 if center else 2,) * len(primes)
    queue = [(prod(primes)**initial[0], initial)]
    seen = {initial}
    while queue:
        n, e = heappop(queue)
        if not covered(primes, tuple(2*x for x in e) if center else e):
            return n
        for j, p in enumerate(primes):
            f = e[:j] + (e[j] + 1,) + e[j+1:]
            if f not in seen:
                seen.add(f)
                heappush(queue, (n*p, f))

primes = (2, 3, 5, 7, 11, 13, 17)
for center in (False, True):
    minima = [min(first(primes[i:i+s], center) for s in range(3, 7))
              for i in (0, 1)]
    print(min(minima), minima[1])
