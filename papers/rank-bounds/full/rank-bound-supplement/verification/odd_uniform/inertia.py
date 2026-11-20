"""Short exact extension of the saved odd neutral-orbit screen.

No profile regeneration, final Gauss search, or rank exclusion is performed.
All files read from the authoritative manuscript remain unchanged.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUPPLEMENT = Path(__file__).resolve().parents[2]
BASELINE = SUPPLEMENT / 'verification/neutral_orbits/check.py'
spec = importlib.util.spec_from_file_location('baseline_neutral_screen', BASELINE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def extra_screen(row, rank):
    """Apply the new congruence only after all prior screens have retained."""
    verdict = old.screen(row, rank)
    if not verdict['keep'] or verdict['reason'] != 'retained':
        return verdict
    support, classes, edges, nines = row
    marked = {p: tuple(v) for p, v in zip((p for p in support if p >= 5), classes)
              if v[0]}
    isolated = [p for p in support if all(p not in edge for edge in edges)]
    base = [(1, prod(support)), (nines, 9)]
    base += [(v[1], p) for p, v in marked.items()]
    active = {p: old.twopart(p-1)*(p-1)
              for p in isolated if p != 3 and p not in marked}
    base += [(n, p*p) for p, n in active.items()]
    if edges:
        base.append((80, 15))

    def mass(records):
        slots = rank-sum(n for n, _ in records)
        return slots, sum((Q(n, b*b) for n, b in records), Q(0))+Q(slots, 225)

    def fits(records):
        slots, upper = mass(records)
        return slots >= 0 and upper >= 1

    fields, strengthened = [], []
    for p in isolated:
        if p == 3:
            continue
        d = old.twopart(p-1)
        if p in marked:
            val, total, act = marked[p]
            neutral = total-act
            if val != 1 or neutral <= 0:
                continue
            quota = max(2*d, ((neutral+2*d-1)//(2*d))*(2*d))
            if p % 4 == 3:
                previous = max(quota, 8) if (act//(p-1)) % 4 == 0 else quota
                quota += (-act-quota) % 8
                if quota > previous:
                    strengthened.append({'prime': p, 'active': act,
                                         'present_neutral': neutral,
                                         'old_quota': previous, 'new_quota': quota})
            fields.append((quota-neutral, 3*p))
        else:
            dense = base + [(max(0, p*(p-1)-active[p]), p*p)]
            if not fits(dense):
                fields.append((2*d, p*p))
    if not fits(base+fields):
        return {'keep': False, 'reason': 'inertia_eight',
                'upper': str(mass(base+fields)[1]), 'strengthened': strengthened}
    # The complete pure population is at least dp at every absent isolated p.
    # These extra pure labels are disjoint from all mixed-component labels.
    pure = [(old.twopart(p-1), p*p) for p in active]
    components = ([{3, 5}]+[{p} for p in support if p not in (3, 5)]
                  if edges else [{p} for p in support])
    c = len(components)
    if c >= 2:
        mixed_floor = min(p*q for j, C in enumerate(components)
                          for E in components[j+1:] for p in C for q in E)
        mixed = [(2*max(c-1, 3), max(15, mixed_floor))]
        if not fits(base+pure+mixed):
            return {'keep': False, 'reason': 'mixed_cut_and_pure_orbits',
                    'upper': str(mass(base+pure+mixed)[1])}
    return {**verdict, 'strengthened': strengthened}
