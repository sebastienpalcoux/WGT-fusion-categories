"""Exact necessary conditions from simultaneous principal-block traces.

Install in rank-bound-supplement/arithmetic.  This module does not change
the original enumeration engine or write certificates.  Its entry point
returns a certificate on exclusion and None on retention, including every
resource limit.  Inputs belong to the manuscript's based-simple setting.
"""
from collections import Counter
from math import gcd, isqrt
from pathlib import Path
import sys

SUPPLEMENT = Path(__file__).resolve().parent.parent
for directory in (SUPPLEMENT, SUPPLEMENT / "engine_src"):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
from block import _class_cost
from degree_rows import residue_cost


class _ResourceLimit(Exception):
    pass


class _Work:
    def __init__(self, maximum):
        self.maximum = maximum
        self.used = 0

    def tick(self):
        if self.used >= self.maximum:
            raise _ResourceLimit("operation limit")
        self.used += 1


def _factor(n, work):
    """Trial division subject to the same operation guard as the search."""
    factors = []
    p = 2
    while p * p <= n:
        work.tick()
        exponent = 0
        while n % p == 0:
            work.tick()
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append((n, 1))
    return factors


def _valuation(n, p, work):
    exponent = 0
    while n % p == 0:
        work.tick()
        n //= p
        exponent += 1
    return exponent


def _is_power(n, p, work):
    if n <= 1:
        return False
    while n % p == 0:
        work.tick()
        n //= p
    return n == 1


def block_trace_failure(row, dimension=None, *, max_codegree=64,
                        max_modulus=16384, max_states=16384,
                        max_operations=200000, diagnostics=None):
    """Return an exact certificate, or None when this test retains the type.

    diagnostics, when supplied, receives at most 32 bounded status records.
    Either work guard returns None for the entire call.  Codegree/modulus
    cutoffs and empty class-option sets skip that branch conservatively.
    None never certifies realizability.

    Force a whole degree class into P_ell only if its codegree is an
    ell-power, or an ell-power row exists and its degree is prime to ell.
    All other eligible classes allow every membership count from zero to
    their full multiplicity.  Thus partial principal membership is retained.
    """
    if (not isinstance(row, (tuple, list)) or len(row) < 2 or row[0] != 1
            or any(type(x) is not int or x <= 1 for x in row[1:])
            or any(a > b for a, b in zip(row, row[1:]))):
        raise ValueError("expected a sorted perfect degree multiset")
    actual_dimension = sum(x * x for x in row)
    if dimension is None:
        dimension = actual_dimension
    if (type(dimension) is not int or dimension != actual_dimension
            or any(dimension % (x * x) for x in row)):
        raise ValueError("dimension or formal-codegree identity is invalid")
    for value in (max_codegree, max_modulus, max_states, max_operations):
        if type(value) is not int or value < 0:
            raise ValueError("resource bounds must be nonnegative integers")

    def note(status, **details):
        if diagnostics is not None and len(diagnostics) < 32:
            diagnostics.append({"status": status, **details})

    work = _Work(max_operations)
    if max_states == 0:
        note("retained", reason="state limit", operations=0)
        return None
    classes = sorted(Counter(row[1:]).items())
    codegrees = {x: dimension // (x * x) for x, _ in classes}
    try:
        factors = _factor(dimension, work)
        source_rows = []
        for y, multiplicity in classes:
            work.tick()
            f = codegrees[y]
            if f > max_codegree:
                continue
            for p, _ in factors:
                if _is_power(f, p, work):
                    source_rows.append((y, multiplicity, f, p))
                    break

        for ell, _ in factors:
            pure = {x for x, _ in classes
                    if _is_power(codegrees[x], ell, work)}
            fixed = set(pure)
            if pure:
                fixed.update(x for x, _ in classes if x % ell)
            possible = set()
            for x, _ in classes:
                work.tick()
                if codegrees[x] % (ell * ell):
                    continue
                remaining = codegrees[x] - 1
                for z, count in classes:
                    work.tick()
                    # Necessary norm cost for ONE degree-x row to be in
                    # P_ell; the z-columns are the WHOLE degree class.
                    # No claim that the whole x- or z-class is principal.
                    cost = _class_cost(z, count, x, ell)
                    if cost < 0 or cost > remaining:
                        break
                    remaining -= cost
                else:
                    possible.add(x)
            if not fixed <= possible:
                # Already excluded by the old block stage; do not reclassify.
                note("skipped", reason="old block condition", block_prime=ell)
                continue

            for y, multiplicity_y, f, p in source_rows:
                work.tick()
                if p == ell:
                    continue
                modulus = ell ** _valuation(dimension // y, ell, work)
                if modulus > max_modulus:
                    note("skipped", reason="modulus limit", block_prime=ell)
                    continue
                terms = []
                branch_skipped = False
                for x, multiplicity in classes:
                    work.tick()
                    if x not in possible:
                        continue
                    g = gcd(x, y)
                    u = x // g
                    # Coprime vanishing, full row norm, P_p support/height.
                    if g == 1 or u * u >= f or codegrees[x] % (p * p):
                        continue
                    options = {}
                    counts = (multiplicity,) if x in fixed else range(multiplicity + 1)
                    for count in counts:
                        work.tick()
                        if count == 0:
                            options[0] = 0
                            continue
                        bound = min(isqrt(count * (f - 1)) // u, count * g - 1)
                        constrained = u % p != 0
                        for trace in range(-bound, bound + 1):
                            work.tick()
                            if constrained and (trace - count * g) % p:
                                continue
                            t, remainder = divmod(trace, count)
                            energy = ((count - remainder) * t * t
                                      + remainder * (t + 1) ** 2)
                            if constrained:
                                for shift in (t - 1, t, t + 1):
                                    work.tick()
                                    cost = residue_cost(
                                        p, (g - shift) % p,
                                        min(count, multiplicity_y))
                                    lower = (count * cost + p - 2) // (p - 1)
                                    energy = max(energy, lower + 2 * shift * trace
                                                 - count * shift * shift)
                                energy += (count * g * g - energy) % p
                            energy *= u * u
                            if energy <= f - 1:
                                residue = (x * u * trace) % modulus
                                options[residue] = min(options.get(residue, f), energy)
                    if not options:
                        note("skipped", reason="empty class options",
                             row_degree=y, block_prime=ell, class_degree=x)
                        branch_skipped = True
                        break
                    terms.append(sorted(options.items()))
                if branch_skipped:
                    continue

                states = {0: 0}
                for options in sorted(terms, key=len):
                    next_states = {}
                    for residue, cost in states.items():
                        for delta, extra in options:
                            work.tick()
                            used = cost + extra
                            if used <= f - 1:
                                combined = (residue + delta) % modulus
                                next_states[combined] = min(
                                    next_states.get(combined, f), used)
                                if len(next_states) > max_states:
                                    raise _ResourceLimit("state limit")
                    states = next_states
                if modulus - 1 not in states:
                    certificate = {
                        "row_degree": y, "codegree": f, "row_prime": p,
                        "block_prime": ell, "modulus": modulus,
                        "reachable_residues": len(states),
                        "operations": work.used,
                    }
                    note("excluded", **certificate)
                    return certificate
    except _ResourceLimit as error:
        note("retained", reason=str(error), operations=work.used)
        return None
    note("retained", reason="no obstruction within scope", operations=work.used)
    return None
