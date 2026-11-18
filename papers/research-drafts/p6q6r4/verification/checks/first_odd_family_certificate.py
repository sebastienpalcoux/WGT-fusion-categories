#!/usr/bin/env python3
"""Exact arithmetic audit of the first residual odd parameter in p^6 q^6 r^4.

This verifies the finite parameter enumeration and a witness for the combined
symbolic necessary conditions.  It neither constructs a category nor certifies
that every other structural obstruction in the manuscript has a realization.
Only the standard library is used; no search is bounded by a time/work limit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path


def iroot(n: int, k: int) -> int:
    lo, hi = 0, 1
    while hi**k <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**k <= n:
            lo = mid
        else:
            hi = mid
    return lo


def prime(n: int) -> bool:
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def u(h: int) -> int:
    j = 1
    while j * (j + 1) < h:
        j += 1
    return j


def audit() -> dict:
    p, q, r = 5, 7, 23
    dimension = p**6 * q**6 * r**4
    assert dimension == 514422090765625 == 22680875**2

    # For an ordered odd triple, p >= 3 and q >= 5; this bounds every r.
    r_bound = iroot(dimension // (3**6 * 5**6), 4)
    primes = [n for n in range(3, r_bound + 1, 2) if prime(n)]
    cases = []
    for a in primes:
        if a**16 > dimension:
            break
        for b in primes:
            if b <= a:
                continue
            if a**6 * b**10 > dimension:
                break
            for c in primes:
                if c <= b:
                    continue
                d = a**6 * b**6 * c**4
                if d > dimension:
                    break
                m = min(c, a * a) // b
                E, H = m * (m + 1) // 2, (c - 1) // 2
                scalar_excluded = (m == 0) or all(
                    (H - 1) * u(h) > h * (h - 1)
                    for h in range(1, E + 1)
                )
                cases.append({
                    "p": a, "q": b, "r": c, "dimension": d,
                    "m": m, "E": E, "H": H,
                    "scalar_excluded": scalar_excluded,
                })
    cases.sort(key=lambda x: x["dimension"])
    residual = [x for x in cases if not x["scalar_excluded"]]
    assert len(cases) == 41
    assert len(residual) == 1
    assert (residual[0]["p"], residual[0]["q"], residual[0]["r"]) == (p, q, r)
    assert all(x["m"] <= 2 for x in cases[:-1])

    # An explicit witness to ALL five conditions of enhanced-criterion.
    # Quotient group G = Z/11; A is the cyclic (11,5,2) difference set.
    H, E, h, k, delta = 11, 6, 5, 2, 0
    A = {1, 3, 4, 5, 9}
    correlations = [len(A & {(a + g) % H for a in A}) for g in range(H)]
    assert correlations == [5] + [2] * 10
    assert u(h) == k and h * (h - 1) - (H - 1) * k == delta
    v = next(j for j in range(1, 20) if j * (j + 1) >= 2 * h)
    conductor_bound = v * (2 * h - v) - (h - k)
    assert delta < conductor_bound
    quotient_checks = []
    for b in range(1, H + 1):
        if H % b:
            continue
        x = [sum(a % b == i for a in A) for i in range(b)]
        e = [0] * b
        lhs = [sum(x[i] * x[(i + j) % b] for i in range(b)) for j in range(b)]
        rhs = [(H // b) * k + (h - k if j == 0 else 0) + e[j]
               for j in range(b)]
        assert lhs == rhs and all(0 <= c <= H // b for c in x)
        quotient_checks.append({"b": b, "x": x, "e": e, "correlation": lhs})

    # Only the prime 3 has odd valuation in h-k=3.  Modulo 11, its
    # multiplicative subgroup has odd order five, so it does not contain -1.
    powers = []
    value = 1
    while True:
        value = value * 3 % 11
        powers.append(value)
        if value == 1:
            break
    assert powers == [3, 9, 5, 4, 1] and 10 not in powers
    assert h < k * (k + 1)
    assert not (h == k * (k + 1))
    assert not (h == k * (k + 1) - 1 and k >= 4)
    assert r < 29  # k_word=4 does not meet the short-support threshold.

    return {
        "status": "PASS_FIRST_ODD_UNIFORM_CRITERION_FRONTIER",
        "scope": "Ordered p^6 q^6 r^4 family only; arithmetic frontier, not category existence.",
        "candidate": {"p": p, "q": q, "r": r, "dimension": dimension,
                      "square_root": isqrt(dimension)},
        "complete_bound": {"r_max": r_bound, "ordered_odd_triples_at_most_candidate": 41,
                           "smaller_triples_excluded_by_scalar_test": 40},
        "mathematical_inputs": ["family:exact-corner", "ar5:thm:uniform-family",
                                "r20g:short-support", "ar5:cor:enhanced-criterion"],
        "combined_criterion_witness": {
            "m": 3, "E": E, "H": H, "h": h, "k": k, "delta": delta,
            "mask_in_Z11": sorted(A), "autocorrelations": correlations,
            "prime_level_bound": conductor_bound, "quotient_checks": quotient_checks,
            "norm": h - k, "powers_of_3_mod_11": powers,
            "critical_endpoint": False, "near_critical_rule_applicable": False,
        },
        "cases": cases,
        "limitations": [
            "No claim about the first unresolved odd dimension among all braided categories.",
            "No claim of existence of a modular table, fusion ring, or category at the candidate.",
            "Conditional crossing and local-block constraints remain necessary at the candidate.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit()
    report["validator_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
