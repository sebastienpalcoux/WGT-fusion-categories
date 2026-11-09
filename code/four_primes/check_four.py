#!/usr/bin/env python3
"""Trace and Gram check on four prime positions in n=P^h Q^k u v.

P is the least prime. The other three prime positions have no input
assignments. Bounds on their ranges are proved in the four-prime reductions.
The supported odd case has h>=3 and k=2. The even case allows h,k>=2.
"""

import argparse
from functools import lru_cache
from itertools import product
from math import isqrt
from time import monotonic


@lru_cache(maxsize=100000)
def prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


@lru_cache(maxsize=512)
def primes(lo, hi):
    return tuple(x for x in range(lo, hi + 1) if prime(x))


def centered(n, t):
    return (n + t // 2) % t - t // 2


def multiplicity_floor(t, g):
    if t <= 2*g:
        return 1
    largest = max(d for d in range(1, min(t - 1, g) + 1) if (t - 1) % d == 0)
    return (t - 1) // largest


def coordinates(t, j, sign, d):
    r = centered(d, t)
    z = centered(-sign*j*d, t)
    a = abs(r)
    admissible = j*a*(t - a) < t*t
    unramified_possible = z == -sign*j*r and j*(t - 1)*r*r < t*t
    return z, admissible, unramified_possible


@lru_cache(maxsize=1024)
def power_classes(prime_modulus, exponent):
    buckets = {}
    for residue in range(1, prime_modulus):
        buckets.setdefault(pow(residue, exponent, prime_modulus), []).append(residue)
    return buckets


def congruent_primes(lower, upper, modulus, residues):
    for residue in residues:
        first = lower + (residue-lower) % modulus
        for q in range(first, upper+1, modulus):
            if prime(q):
                yield q


def orbit_deficit(t, j):
    return 3*j <= t+1 and (t-1) % j == 0 and ((t-1)//j) % 2 == 1


def vacuum_cases(P, h, k, orbit_window=False):
    odd = P != 2
    step = 2 if odd else 1
    if not odd:
        upper = 2**(h + 1)
    elif P == 3:
        upper = (9*P**(h-1) + 2) // 4
    else:
        upper = P**(h-1)*(2*(P-1)**2+P**2) // (2*(P-1)**2)

    for u in primes(max(5, P + 1), upper):
        if odd and u in (5,7,13):
            continue
        for ju in range(step, u, step):
            if (u-1) % ju and (u+1) % ju:
                continue
            if orbit_window and not orbit_deficit(u, ju):
                continue
            quota = max(2, (u+1)//ju - 1)
            v_top = P**h*quota + 1
            if orbit_window:
                v_top = min(v_top, P*u//(P-1))
            for v in primes(u+1, v_top):
                if odd and v in (5,7,13):
                    continue
                top = v*quota + 1
                if v*(P-1) > P*u:
                    top = min(top, P*u//(P-1))
                inverse = pow(P**h*v*ju,-1,u)
                for sign in (-1,1):
                    classes = power_classes(u,k).get((-sign*inverse) % u, ())
                    for q in congruent_primes(P+1,top,u,classes):
                        if q == v:
                            continue
                        if odd and u*(P-1) > P*min(P**h, q**k):
                            continue
                        core = P**h*q**k
                        jv = (-sign*pow(core*u,-1,v)) % v
                        if jv % step or ((v-1)%jv and (v+1)%jv):
                            continue
                        if orbit_window and not orbit_deficit(v, jv):
                            continue
                        n = core*u*v
                        if 0 < n - sign - core*v*ju - core*u*jv < n:
                            yield q,u,v,ju,jv,sign


def axis_options(base, height, other, t, w, j, sign, step):
    for i in range(1, height+1):
        d = base**i*w
        z, admissible, unramified = coordinates(t,j,sign,d)
        if admissible:
            lower = max(step, multiplicity_floor(t,w), multiplicity_floor(other,w),
                        1 if unramified else 2)
            yield d, z, lower + (-lower) % step


def axis_trace_costs(options, modulus, cap, step):
    """Least squared cost for each attainable trace of a nonempty axis."""
    states = {}
    for d, z, lower in options:
        square = z*z
        updated = dict(states)
        for previous, spent in [(0,0), *states.items()]:
            for m in range(lower, (cap-spent)//square+1, step):
                residue = (previous+m*d*z) % modulus
                updated[residue] = min(updated.get(residue,cap+1), spent+m*square)
        states = updated
    return states


def axis_floor(base, height, other, t, w, j, sign, step, trace):
    cap = j*(t+1-j)
    options = axis_options(base,height,other,t,w,j,sign,step)
    if trace:
        return axis_trace_costs(options,other,cap,step).get(sign*j % other, cap+1)
    return min([cap+1, *(lower*z*z for _,z,lower in options)])


def common_terms(P,h,k,case,caps):
    """Common powers meeting both residual axis budgets."""
    q,u,v,ju,jv,sign = case
    step = 1 if P == 2 else 2
    for i,l in product(range(1,h+1),range(1,k+1)):
        x,y = P**i,q**l
        d,g = x*y,min(x,y)
        zu, admissible_u, unram_u = coordinates(u,ju,sign,d)
        zv, admissible_v, unram_v = coordinates(v,jv,sign,d)
        if not admissible_u or not admissible_v:
            continue
        usable = int(unram_u and (P-1)*u <= P*g) + int(unram_v and (P-1)*v <= P*g)
        if usable == 0:
            continue
        lower = max(step, 3-usable, multiplicity_floor(u,g), multiplicity_floor(v,g))
        lower += (-lower) % step
        if lower*zu*zu <= caps[0] and lower*zv*zv <= caps[1]:
            yield i,l,zu,zv,lower


def gram_possible(P,h,k,case,trace):
    q,u,v,ju,jv,sign = case
    step = 1 if P == 2 else 2
    caps = []
    for t,w,j in ((u,v,ju),(v,u,jv)):
        cap = j*(t+1-j)
        cap -= axis_floor(P,h,q,t,w,j,sign,step,trace)
        cap -= axis_floor(q,k,P,t,w,j,sign,step,trace)
        if cap <= 0:
            return False
        caps.append(cap)

    choices = [(zu,zv,lower) for _,_,zu,zv,lower in common_terms(P,h,k,case,caps)]
    if not any(zu*zv < 0 for zu,zv,_ in choices):
        return False

    @lru_cache(None)
    def fill(i, remain_u, remain_v, target):
        if target == 0:
            return True
        if i == len(choices):
            return False
        zu,zv,lower = choices[i]
        if fill(i+1,remain_u,remain_v,target):
            return True
        upper = min(remain_u//(zu*zu), remain_v//(zv*zv))
        for m in range(lower, upper+1, step):
            if fill(i+1,remain_u-m*zu*zu,remain_v-m*zv*zv,target-m*zu*zv):
                return True
        return False

    return fill(0,*caps,-ju*jv)


def check(P,h,k,orbit_window=False):
    if not prime(P) or h<2 or k<2 or (P!=2 and (h<3 or k!=2)):
        raise ValueError("choose a prime P and heights h,k>=2; odd P requires h>=3, k=2")
    if orbit_window and P == 2:
        raise ValueError("the odd-orbit theorem requires odd P")
    vacuum, gram, open_cases = 0,0,[]
    for case in vacuum_cases(P,h,k,orbit_window):
        vacuum += 1
        if gram_possible(P,h,k,case,False):
            gram += 1
            if gram_possible(P,h,k,case,True):
                open_cases.append(case)
    return vacuum,gram,open_cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("P", type=int)
    parser.add_argument("h", type=int)
    parser.add_argument("k", type=int)
    parser.add_argument("--orbit-window", action="store_true",
                        help="use the proved odd bounds v <= P*u/(P-1) and on both orbit deficits")
    args = parser.parse_args()
    start = monotonic()
    try:
        vacuum,gram,residual = check(args.P,args.h,args.k,args.orbit_window)
    except ValueError as error:
        parser.error(str(error))
    print(f"P={args.P}, h={args.h}, k={args.k}: {vacuum} signed vacuum classes, "
          f"{gram} after the two norm budgets and Gram test, {len(residual)} after the two trace minima.")
    if residual:
        print("First residual (Q,u,v,j_u,j_v,sign):", residual[0])
    print(f"Elapsed: {monotonic()-start:.3f} seconds")
    raise SystemExit(2 if residual else 0)
