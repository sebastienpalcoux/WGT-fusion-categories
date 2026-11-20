"""Validated terminal necessary conditions for ranks 1 through 44.

This module imports no old driver, receipt or published census.  Invalid
producer output raises ValueError and cannot be recorded as an exclusion.
The square frontier restrictions below are proved throughout this range.
Resource-limited optional tests retain the candidate.
"""
from collections import Counter
from math import gcd, isqrt, prod
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
for directory in (ROOT, ROOT / "engine_src", ROOT / "phase-mass"):
    value = str(directory)
    if value not in sys.path:
        sys.path.insert(0, value)

from screen import OBSTRUCTIONS as BASE_OBSTRUCTIONS, exclusion_reason, factor
from engine.primes import primes_up_to
from block import block_failure
from arithmetic.block_trace import block_trace_failure
from phase_norms import phase_columns_possible
from prime_square import prime_square_failure
from prime_cut import prime_cut_failure
from principal_overlap import phase_failure
from degree_slices import slice_failure
from fusion_support import fusion_failure
from screen import tiles_in_fiber


RANK = 44
OBSTRUCTIONS = BASE_OBSTRUCTIONS | {
    "block_entry", "simultaneous_block_trace", "frontier", "phase_norm",
    "one_packet", "four_prime", "effective_prime_height",
    "four_prime_isolation", "prime_cut", "prime_square_trace", "degree_slices", "principal_overlap", "fusion_phases",
}


def validate_row(row, task, rank=RANK):
    """Return D after checking the exact row and its task coordinates.

    Tasks have squarefree q, positive B and optional sorted denominator
    prefix.  A nonsquare task additionally fixes the top multiplicity m.
    Only mathematical candidates emitted inside that task may be screened.
    """
    if type(rank) is not int or not 1 <= rank <= RANK:
        raise ValueError("terminal checker supports ranks 1 through 44")
    if not isinstance(task, dict):
        raise ValueError("expected a task object")
    b, q = task.get("B"), task.get("q")
    if (type(b) is not int or type(q) is not int or b < 2 or q < 1
            or any(e != 1 for _, e in factor(q))):
        raise ValueError("invalid denominator or squarefree kernel")
    f = q * b * b
    if not 8 <= f <= rank - 2 or f in (8, 9):
        raise ValueError("task outside the proved starting cases through rank 44")
    if (not isinstance(row, (list, tuple)) or len(row) != rank
            or any(type(d) is not int for d in row) or row[0] != 1
            or row[1] < 6 or any(a > c for a, c in zip(row, row[1:]))
            or gcd(*row[1:]) != 1):
        raise ValueError("invalid sorted simple-degree profile")
    M = b * row[-1]
    D = q * M * M
    if sum(d * d for d in row) != D or any(M % d for d in row):
        raise ValueError("profile violates dimension or denominator identity")
    prefix = task.get("prefix", [b])
    if (not isinstance(prefix, (list, tuple)) or not prefix
            or len(prefix) > rank - 3 or prefix[0] != b
            or any(type(z) is not int or z < b for z in prefix)
            or any(a > c for a, c in zip(prefix, prefix[1:]))):
        raise ValueError("invalid task denominator prefix")
    if any(d * z != M for d, z in zip(reversed(row), prefix)):
        raise ValueError("degree list outside its denominator prefix")
    if q != 1:
        m = task.get("m")
        if type(m) is not int or not 1 <= m < rank or row.count(row[-1]) != m:
            raise ValueError("nonsquare row outside its top multiplicity")
    return D


def _prime_power(n, p):
    if n < 2 or n % p:
        return False
    while n % p == 0:
        n //= p
    return n == 1


def square_frontier_failure(row, M):
    """Support and present-prime restrictions proved through rank 44."""
    if not 1 <= len(row) <= RANK:
        raise ValueError("square frontier restrictions stop at rank 44")
    support = {p for p, _ in factor(M)}
    if len(support) != 4:
        return True
    counts = Counter(M // d for d in row[1:])
    if counts[5] and (M % 25 == 0 or counts[5] != 5):
        return True
    return any(counts[p] and (M % (p * p) == 0 or counts[p] < p - 1)
               for p in primes_up_to(len(row)) if p >= 7)


def one_packet_failure(row, M):
    """Uniform conditional norm test at a valuation-one class of size p."""
    counts = Counter(row[1:])
    for p, exponent in factor(M):
        n = M // p
        if p < 5 or exponent != 1 or n <= p + 1 or counts[n] != p:
            continue
        norm = 0
        for d, multiplicity in counts.items():
            if d == n or d % p == 0:
                continue
            if n % d:
                raise ValueError("p-free degree does not divide the p-free part of M")
            residue = pow(n // d, -1, p)
            norm += multiplicity * min(residue, p - residue) ** 2
        if norm > p - 1:
            return True
    return False


def phase_norm_failure(row, M):
    """Rank-independent integral phase norms with conservative work skips."""
    rank = len(row)
    counts = Counter(row[1:])
    factors = factor(M)
    pure = {p for p, _ in factors
            if any(_prime_power(M // d, p) for d in counts)}
    for p, _ in factors:
        if p == 2:
            continue
        h = (p - 1) // 2
        p_free = [d for d in counts if d % p]
        groups = tuple(
            (d, n // h, rank >= 1 + p * h,
             int(all(gcd(d, t) > 1 for t in p_free)))
            for d, n in sorted(counts.items())
            if n >= h and (M // d) % p == 0
            and all(d % t == 0 for t in pure if t != p)
        )
        choices = prod((n + 1) * (1 + n // p if high else 1)
                       for _, n, high, _ in groups)
        if choices > 100000:
            continue
        if not phase_columns_possible(groups, p, M, 1, row[1]):
            return True
    return False


def reason_for(row, task, rank=RANK, *, adapter_failure=None):
    """First exact obstruction, or None (not a realizability certificate).

    adapter_failure, if supplied, is called as (row, task, rank) after
    principal-block admission.  It must return None or a known reason;
    errors and unfinished computations must propagate to the runner.
    """
    D = validate_row(row, task, rank)
    b, q = task["B"], task["q"]
    M = b * row[-1]
    factors = tuple(factor(D))
    if slice_failure(row, D, [p for p, _ in factors]):
        return "degree_slices"
    if phase_failure(row, D, factors):
        return "principal_overlap"
    if fusion_failure(row):
        return "fusion_phases"
    # Two distinct valuation-two primes on four-prime support force a square.
    if q != 1 and len(factors) == 4 and sum(e == 2 for _, e in factors) >= 2:
        return "four_prime_isolation"
    if prime_square_failure(row, q, M) is not None:
        return "prime_square_trace"
    if block_failure(row, D):
        return "block_entry"
    if block_trace_failure(row, D) is not None:
        return "simultaneous_block_trace"
    if q > 1 and not tiles_in_fiber(Counter(row).items(), q, M):
        return "defect_one"
    if prime_cut_failure(row, M, tuple(p for p, _ in factors)) is not None:
        return "prime_cut"
    if adapter_failure is not None:
        reason = adapter_failure(row, task, rank)
        if reason is not None:
            if reason not in OBSTRUCTIONS:
                raise ValueError("unknown adapter exclusion reason: " + str(reason))
            return reason
    if M % 2:
        return "frontier"
    if q == 1:
        if square_frontier_failure(row, M):
            return "frontier"
        if one_packet_failure(row, M):
            return "one_packet"
    reason = exclusion_reason(row, dimension=D)
    if reason is not None:
        return reason
    if q == 1 and phase_norm_failure(row, isqrt(D)):
        return "phase_norm"
    return None
