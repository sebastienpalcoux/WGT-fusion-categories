"""Exact simultaneous CRT constraints on rational local Weyl entries.

Input parts are the pairwise coprime parts of one nonunit denominator over
the twist components. This checks necessary entry constraints only; it is
not a modular-data construction. Optional parity applies to odd parts.
"""
from fractions import Fraction as F
from math import gcd, prod


def _setup(parts, parities):
    parts = tuple(parts)
    if not parts or any(type(u) is not int or u < 1 for u in parts):
        raise ValueError('positive integral component parts are required')
    if any(gcd(u, v) != 1 for j, u in enumerate(parts) for v in parts[j+1:]):
        raise ValueError('component parts must be pairwise coprime')
    b = prod(parts)
    if b < 2:
        raise ValueError('nonunit denominators must exceed one')
    if parities is None:
        parities = (None,)*len(parts)
    parities = tuple(parities)
    if len(parities) != len(parts) or any(x not in (None, 0, 1) for x in parities):
        raise ValueError('one parity in {None,0,1} per component is required')
    if any(u % 2 == 0 and parity is not None for u, parity in zip(parts, parities)):
        raise ValueError('2-adic parity may only be prescribed at odd parts')
    residues = tuple(pow(b//u, -1, u) if u > 1 else 0 for u in parts)
    qs = tuple(F(r, u) for r, u in zip(residues, parts))
    K = sum(qs, F()) - F(1, b)
    assert K.denominator == 1
    return parts, b, residues, qs, int(K), parities


def coupled_weyl_minimum(parts, parities=None, *, forced_negative=(), forced_positive=()):
    """Optimal total squared norm and exact feasibility of entry constraints.

    For each u>1 the two lifts are r/u and r/u-1. Exactly K are negative.
    Their positive sum must be <1, equivalently every union entry has
    absolute value <1. With forced choices (including 2-adic parity), both
    the smallest norm and smallest positive sum select the largest free
    residues for the remaining negative positions.

    Return None if these necessary constraints are inconsistent. Otherwise
    return a certificate, with exact Fraction objects, of the minimizing
    necessary-condition solution. No claim of realization is made.
    """
    parts, b, residues, qs, K, parities = _setup(parts, parities)
    active = {j for j, u in enumerate(parts) if u > 1}
    negative, positive = set(forced_negative), set(forced_positive)
    if not (negative | positive).issubset(range(len(parts))):
        raise ValueError('invalid component index')
    if negative - active:
        return None
    for j, (u, r, parity) in enumerate(zip(parts, residues, parities)):
        if parity is None:
            continue
        if u == 1:
            if parity:
                return None
            continue
        # For odd u, the negative lift r-u reverses the numerator parity.
        (negative if (r & 1) != parity else positive).add(j)
    if negative & positive:
        return None
    free = active - negative - positive
    needed = K-len(negative)
    if needed < 0 or needed > len(free):
        return None
    negative.update(sorted(free, key=lambda j: qs[j], reverse=True)[:needed])
    lifts = tuple(q-int(j in negative) for j, q in enumerate(qs))
    positive_mass = sum((x for x in lifts if x > 0), F())
    if positive_mass >= 1:
        return None
    assert sum(lifts, F()) == F(1, b)
    return {'parts': parts, 'denominator': b, 'negative_count': K,
            'negative_indices': tuple(sorted(negative)), 'lifts': lifts,
            'positive_mass': positive_mass,
            'squared_norm': sum((x*x for x in lifts), F())}


def coupled_weyl_entry_floors(parts, parities=None):
    """Minimum squared magnitude of each individual entry, under coupling.

    The minimizing choice may differ between components. These independent
    lower bounds are safe for separate column budgets; adding them only
    weakens the stronger jointly minimized total norm.
    """
    parts = tuple(parts)
    best = coupled_weyl_minimum(parts, parities)
    if best is None:
        return None
    floors, alternatives = [], []
    for j, u in enumerate(parts):
        if u == 1:
            floors.append(F())
            alternatives.append((F(),))
            continue
        lifts = set()
        for negative in (False, True):
            answer = coupled_weyl_minimum(parts, parities,
                forced_negative=(j,) if negative else (),
                forced_positive=() if negative else (j,))
            if answer is not None:
                lifts.add(answer['lifts'][j])
        assert lifts
        floors.append(min(x*x for x in lifts))
        alternatives.append(tuple(sorted(lifts)))
    return {'parts': parts, 'entry_squared_floors': tuple(floors),
            'entry_alternatives': tuple(alternatives),
            'joint_squared_floor': best['squared_norm']}
