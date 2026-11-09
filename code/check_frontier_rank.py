"""Optional exact verification of the two rank-certificate tables, Python 3.9+."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import gcd, isqrt

D = 5400**2
G = 900**2
M0 = 6703248
degrees = [d for d in range(675, 9, -1) if 5400 % d == 0
           and sum(d % p == 0 for p in (2, 3, 5)) >= 2]
quartets = sorted(d for d in degrees if 360 % d == 0)
active = (24, 72, 120, 216, 360)
joint = (24, 72, 120, 360)
ten_sums = {0}
for _ in range(10):
    ten_sums = {s + (d // 24)**2 for s in ten_sums for d in active}
ten_sums = sorted(ten_sums, reverse=True)
plans = defaultdict(list)
for u, v in combinations_with_replacement(joint, 2):
    for h, i in combinations_with_replacement(active, 2):
        plans[(4*u*u + 4*v*v + h*h + i*i)//576].append((u, v, h, i))

# Bit t records whether the suffix of allowed decrements can sum to t.
limit = M0 + (100-66)*G - D
mask = (1 << (limit+1)) - 1
suffix = {}
bits = 1
for j in reversed(range(len(degrees))):
    step = G - degrees[j]**2
    while step <= limit:
        bits = (bits | (bits << step)) & mask
        step *= 2
    suffix[j] = bits.to_bytes((limit+8)//8, 'little')

def tails(rem, start=0):
    if rem == 0:
        yield ()
        return
    if not (suffix[start][rem//8] >> (rem%8)) & 1:
        return
    for j in range(start, len(degrees)):
        cost = G - degrees[j]**2
        if cost > rem:
            break
        for tail in tails(rem-cost, j):
            yield (degrees[j],) + tail

def profiles(rank):
    budget = M0 + (rank-66)*G - D
    for b,c in combinations_with_replacement(quartets, 2):
        for o,w,y in product((15,45,75,135,225,675),(10,20,50,100),(12,18,36)):
            base = budget - 4*(2*360**2-b*b-c*c) - (675**2-o*o)
            base -= (100**2-w*w) + (36**2-y*y)
            for s in ten_sums:
                rem = base - 2304*(2250-s)
                if rem < 0:
                    break
                for tail in tails(rem):
                    if len(tail) <= rank-66 and tail.count(54)%2 == 0:
                        if (tail.count(18)+(y==18))%2 == 0:
                            yield b,c,o,w,y,s,tail

@lru_cache(None)
def sum_of_squares(weights, rem):
    if not weights:
        return rem == 0
    q = weights[0]
    return any(sum_of_squares(weights[1:], rem-q*x*x)
               for x in range(isqrt(rem//q)+1))

@lru_cache(None)
def disjoint_fibers(groups, active_mass):
    good = []
    for bits in range(1, 1 << len(groups)):
        members = [dm for j,dm in enumerate(groups) if bits & (1 << j)]
        if sum(m for d,m in members) < 8:
            continue
        weights = tuple(sorted((m*d*d for d,m in members), reverse=True))
        total = active_mass + sum(weights)
        if total % (25*gcd(*(d for d,m in members))**2):
            continue
        if sum_of_squares(weights, total//25):
            if any(bits & earlier == 0 for earlier in good):
                return True
            good.append(bits)
    return False

for rank in range(95, 101):
    count, maximum = 0, None
    results = []
    for b,c,o,w,y,s,tail in profiles(rank):
        count += 1
        eligible = sum(d in active for d in tail)
        maximum = eligible if maximum is None else max(maximum, eligible)
        if not plans[s]:
            continue
        spare = [*tail, o]
        if y == 18:
            spare.remove(18)
        groups = [(d,4) for d in (12,15,b,c)]
        for d,m in Counter(d for d in spare if 360 % d == 0).items():
            if d in (18,54):
                assert m%2 == 0
                groups.extend([(d,2)]*(m//2))
            else:
                groups.extend([(d,1)]*m)
        for u,v,h,i in plans[s]:
            if disjoint_fibers(tuple(sorted(groups)), 8*(u*u+v*v)):
                results.append((b,c,o,w,y,s,tail,u,v,h,i))
    print(rank, count, maximum, len(results), flush=True)
    for result in results:
        print(result)
