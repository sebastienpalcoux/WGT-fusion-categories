"""Check the seven displayed trace obstructions and the two fibers of each survivor."""
from collections import Counter
from itertools import product
import json
from pathlib import Path

certificates = json.loads((Path(__file__).resolve().parents[1] / "certificates" / "certificates.json").read_text())


def centered(n, prime):
    r = n % prime
    return r if 2*r <= prime else r-prime


def orbit_floor(t, g):
    if t <= 2*g:
        return 1
    return min((t-1)//e for e in range(1,g+1) if (t-1) % e == 0)


def check_trace(row):
    h,k,signed_q,u,v,j_u,j_v,t,cap,printed,needed = row["data"]
    q = abs(signed_q)
    eps = 1 if signed_q > 0 else -1
    w,j = (v,j_u) if t == u else (u,j_v)
    choices = []
    for i in range(1,h+1):
        d = (2**i)*w
        r = centered(d,t)
        a = centered(-eps*j*d,t)
        if j*abs(r)*(t-abs(r)) >= t*t:
            continue
        free = a == -eps*j*r and j*(t-1)*r*r < t*t
        lower = max(orbit_floor(t,w),orbit_floor(q,w),1 if free else 2)
        choices.append([(0,0)]+[(m*a*a,m*d*a) for m in range(lower,cap//(a*a)+1)])
    achievable = {sum(x[1] for x in v) % q for v in product(*choices)
                  if 0 < sum(x[0] for x in v) <= cap}
    assert achievable == set(printed) and needed == eps*j % q and needed not in achievable


def check_fibers(profile):
    b,c,o,y = (profile[k] for k in ("b","c","o","y"))
    u,v,h,i = profile["plan"]
    assert profile["rank"] == 100 and u == v == 360 and h != i
    assert 4*u*u+4*v*v+h*h+i*i == 576*profile["S"]
    spare = profile["tail"]+[o]
    if y == 18:
        spare.remove(18)
    available = Counter((degree,4) for degree in (12,15,b,c))
    for degree, number in Counter(d for d in spare if 360 % d == 0).items():
        size = 2 if degree in (18,54) else 1
        assert number % size == 0
        available[degree,size] += number//size
    witness = profile["witness"]
    used = Counter()
    for side in ("left","right"):
        groups = witness[side]
        used.update(map(tuple,groups))
        assert sum(m for d,m in groups) >= 8
        weights, coordinates = witness[side+"_quadratic"]
        assert weights == sorted(m*d*d for d,m in groups)
        assert len(weights) == len(coordinates) and all(isinstance(x,int) and x>=0 for x in coordinates)
        assert 8*(u*u+v*v)+sum(weights) == 25*sum(w*x*x for w,x in zip(weights,coordinates))
    assert all(m <= available[group] for group,m in used.items())


traces = [row for block in certificates["dyadic"] for row in block["trace_rows"]]
profiles = certificates["rank100"]["survivors"]
assert len(traces) == 7 and len(profiles) == 5
for row in traces:
    check_trace(row)
for profile in profiles:
    check_fibers(profile)
print("7 trace obstructions and 5 pairs of disjoint fiber witnesses verified")
