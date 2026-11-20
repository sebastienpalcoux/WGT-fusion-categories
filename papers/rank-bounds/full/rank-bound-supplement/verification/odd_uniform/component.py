"""Exact linear mass envelopes on saved odd-square reservations only.

This is a necessary-condition check, not a Gauss or denominator census.
All rank-derived assertions are explicit; the saved inputs are left unchanged.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod, ceil
from residue import height_lower_bounds, reservation_failure, support_budget_failure
from inertia import old
from global_weyl import global_constraints
from coupled_weyl import coupled_weyl_entry_floors


def rounded(n, modulus):
    return ((n+modulus-1)//modulus)*modulus


def local_norm_upper(p, active_count, mass_ceiling):
    """Maximize the proved concave bound with E<=active_count/p^4."""
    x = min(Q(active_count, p**4), mass_ceiling, Q(p*p-1, 2*p*p))
    return Q(p+1, p-1)*x-Q(p*p, (p-1)**2)*x*x


FLOOR_CAP = 10001


@lru_cache(None)
def nine_witness_quota(p, multiplicity):
    """Paired unit-weight form of the existing weighted prime-free bound."""
    numerator = 0
    for a in range(1, (p+1)//2):
        u = next(k for k in range(3, 3*p, 3) if (k+a) % p == 0)
        l = u-3*p
        best = min((multiplicity-j)*l*l+j*u*u
                   +Q(((multiplicity-j)*l+j*u)**2, 81-multiplicity)
                   for j in range(0, multiplicity+1, 2))
        numerator += 2*ceil(best/162)
    return 2*ceil(Q(numerator, 2*sum(a*a for a in range(1, (p+1)//2))))


@lru_cache(None)
def supported_denominators(support):
    result = []
    for b in range(15, FLOOR_CAP, 2):
        if b in support:
            continue
        n = b
        for p in support:
            while n % p == 0:
                n //= p
        if n == 1:
            result.append(b)
    return tuple(result)


@lru_cache(None)
def field_entry_constraints(b, components, minimal, field):
    parts, parities = [], []
    for component in components:
        part = 1
        for p in component:
            remaining = b
            while remaining % p == 0:
                part *= p
                remaining //= p
        parts.append(part)
        p = component[0]
        parities.append(int(field == component) if len(component) == 1 and p in minimal
                        else None)
    return coupled_weyl_entry_floors(tuple(parts), tuple(parities))


def component_envelope(row, rank):
    support, classes, edges, nines = row
    marked = {p: tuple(v) for p, v in zip((p for p in support if p >= 5), classes)
              if v[0]}
    if (3 not in support or nines <= 0 or edges not in ([], [[3, 5]])
            or any(v[0] != 1 for v in marked.values())):
        return {'keep': True, 'reason': 'outside_scope'}
    if len(support) == 3 and marked:
        return {'keep': False, 'reason': 'prime_denominator_with_three_primes'}
    if any(p in edge for p in marked for edge in edges):
        return {'keep': False, 'reason': 'prime_denominator_not_isolated'}
    failure = reservation_failure(row) or support_budget_failure(row)
    if failure:
        return {'keep': False, **failure}
    isolated = [p for p in support if all(p not in edge for edge in edges)]
    absent = [p for p in isolated if p != 3 and p not in marked]
    minimal = {p for p, (v, m, act) in marked.items()
               if p in isolated and v == 1 and m == 2*p and act == 2*(p-1)}
    heights = height_lower_bounds(row)
    global_failure, global_heights = global_constraints(row)
    if global_failure:
        return {'keep': False, **global_failure}
    for p, h in global_heights.items():
        heights[p] = max(h, heights.get(p, 1))
    residual_mass = 1-Q(nines, 81)-sum((Q(v[1], p*p) for p, v in marked.items()), Q(0))
    if residual_mass <= 0:
        return {'keep': False, 'reason': 'component_mass', 'upper': str(residual_mass),
                'maximizing_relaxation': None}
    best, best_case = None, None
    # Height>=3 is relaxed to3. Exact height2 forbids every denominator27|b.
    lower3 = max(2, heights.get(3, 1))
    three_heights = ((max(3, lower3),) if lower3 >= 3 or
                     (len(support) == 3 and nines >= 9) else (2, 3))
    for dense_bits in product((False, True), repeat=len(absent)):
        dense = {p for p, bit in zip(absent, dense_bits) if bit}
        sparse = set(absent)-dense
        for height3 in three_heights:
            lower_heights = {p: max(heights.get(p, 1),
                                  height3 if p == 3 else 3 if p in sparse else
                                  2 if p in absent or len(support) == 3 else 1)
                             for p in support}
            lower_M = prod(p**lower_heights[p] for p in support)
            components = [{3, 5}] + [{p} for p in support if p not in (3, 5)] if edges else [{p} for p in support]
            component_key = tuple(tuple(sorted(C)) for C in components)
            component_index = {tuple(sorted(C)): j for j, C in enumerate(components)}
            specifications = []
            for component in components:
                p = min(component)
                modulus = max(8, 2*old.twopart(p-1)) if len(component) == 1 else 8
                if 3 in component:
                    free = (6 if height3 == 2 else 8) if nines >= 42 else 0
                    overlap = min(8, free) if edges else 0
                    count = nines + (80 if edges else 0) + free-overlap
                elif p in marked:
                    count = marked[p][1]
                elif p in sparse:
                    count = old.twopart(p-1)*(p+1)
                else:
                    count = 2*p*(p-1)
                specifications.append((component, p, count, rounded(count, modulus)))
            used = sum(t for _, _, _, t in specifications)
            slots = rank-1-used
            if slots < 0:
                continue
            limits = {p: local_norm_upper(p, slots+total, residual_mass)
                      for _, p, _, total in specifications if p in absent}

            def allowed(b):
                if b < 15 or b in support or b == 9:
                    return False
                if height3 == 2 and b % 27 == 0:
                    return False
                if any(b % (p*p) == 0 for p in marked):
                    return False
                if any(b % p == 0 and b % (p*p) != 0 for p in sparse):
                    return False
                for p, (_, m, act) in marked.items():
                    if b % p == 0:
                        k = b//p
                        a = pow(k, -1, p)
                        a = min(a, p-a)
                        j = act//(p-1)
                        H = 2*p*j-j*j-m
                        if 2*a*a+4 > H:
                            return False
                for p, bound in limits.items():
                    if b % p:
                        continue
                    k, u = b, 1
                    while k % p == 0:
                        k //= p
                        u *= p
                    if k == 1:
                        continue
                    a = pow(k, -1, u)
                    a = min(a, u-a)
                    if 2*a*a*bound.denominator > u*u*bound.numerator:
                        return False
                return True

            possible = [b for b in supported_denominators(tuple(support)) if allowed(b)]
            # If no candidate below FLOOR_CAP exists, that cap itself remains
            # a valid lower bound. This weakens the mass test; it never rejects
            # a branch merely because a denominator was above the cap.
            def field_allowed(b, component):
                for p in minimal:
                    if b % p:
                        continue
                    residue = pow(b//p, -1, p)
                    odd = component == {p}
                    coefficient = residue if residue % 2 == odd else residue-p
                    if coefficient*coefficient > p-(10 if odd else 7):
                        return False
                coupled = field_entry_constraints(
                    b, component_key, tuple(sorted(minimal)), tuple(sorted(component)))
                if coupled is None:
                    return False
                for p, bound in limits.items():
                    j = component_index[(p,)]
                    if coupled['parts'][j] == b:
                        continue  # The parabolic entry bound is used only on mixed denominators.
                    if 2*coupled['entry_squared_floors'][j] > bound:
                        return False
                for p, (_, m, act) in marked.items():
                    if b % p:
                        continue
                    j = act//(p-1)
                    budget = 2*p*j-j*j-m
                    coordinate = component_index[(p,)]
                    if 2*p*p*coupled['entry_squared_floors'][coordinate]+4 > budget:
                        return False
                return True

            def floor(test, component):
                return next((b for b in possible
                             if test(b) and field_allowed(b, component)), FLOOR_CAP)

            records = []
            for component, p, count, total in specifications:
                b = floor(lambda b: any(b % p == 0 for p in component), component)
                if 3 in component:
                    mass = Q(nines, 81) + (Q(80, 225) if edges else 0)
                    if nines >= 42:
                        free = 6 if height3 == 2 else 8
                        bf = floor(lambda b: b % (3**height3) == 0, component)
                        if edges:
                            # Of the eighty edge labels, at least seventy-two
                            # are five-active, hence not three-free. Allow the
                            # remaining eight to overlap the free reservation.
                            overlap = min(8, free)
                            mass += Q(free, bf*bf)-Q(overlap, 225)
                        else:
                            mass += Q(free, bf*bf)
                elif p in marked:
                    mass = Q(count, p*p)
                elif p in sparse:
                    # Every p-field denominator has p²|b in this regime.
                    mass = Q(count, p**4)
                else:
                    mass = Q(count, p**4)
                if p in absent:
                    # The largest active layer contains >=2(p-1) labels,
                    # at B >= p^ceil((v_p(M)+1)/2), in both density cases.
                    high = p**((lower_heights[p]+2)//2)
                    quota = 2*(p-1)
                    assert quota <= count and high >= p*p
                    mass -= Q(quota, p**4)-Q(quota, high*high)
                mass += Q(total-count, b*b)
                records.append((total, mass, b, p not in minimal))
            upper = Q(1, lower_M*lower_M) + sum((m for _, m, _, _ in records), Q(0))
            # Minimal valuation-one components have exactly 2p+2 labels.
            # Every unassigned label therefore belongs to another component.
            upper += Q(slots, min(b for _, _, b, open_ in records if open_)**2)
            omission_upper = None
            if len(support) == 3 and nines >= 9 and not marked:
                # Different omitted primes give disjoint nonunit degrees:
                # each such degree uses both of the other support primes.
                separate = [(nines, Q(nines, 81))]
                for p in absent:
                    count = (2*p*(p-1) if p in dense else old.twopart(p-1)*p)
                    high = p**((lower_heights[p]+2)//2)
                    quota = 2*(p-1)
                    mass = Q(count, p**4)-Q(quota, p**4)+Q(quota, high*high)
                    separate.append((count, mass))
                if edges:
                    # These are five-active; all selected omission witnesses
                    # and all three-free labels have three-primary twists.
                    # Discard the remaining eight edge positions altogether.
                    separate.append((72, Q(72, 225)))
                if nines >= 42:
                    comp3 = {3, 5} if edges else {3}
                    bf = floor(lambda b: b % (3**height3) == 0, comp3)
                    separate.append((8, Q(8, bf*bf)))
                for p in support:
                    if p == 3:
                        continue
                    count = nine_witness_quota(p, nines)
                    factor = 9*p**lower_heights[p]
                    bf = next((b for b in possible if b % factor == 0), FLOOR_CAP)
                    separate.append((count, Q(count, bf*bf)))
                left = rank-1-sum(n for n, _ in separate)
                if left < 0:
                    continue
                least = possible[0] if possible else FLOOR_CAP
                omission_upper = (Q(1, lower_M*lower_M)
                                  +sum((m for _, m in separate), Q(0))
                                  +Q(left, least*least))
                upper = min(upper, omission_upper)
            if best is None or upper > best:
                best, best_case = upper, {'dense': sorted(dense), 'height3': height3,
                                         'component_counts': [n for n, _, _, _ in records],
                                         'component_floors': [b for _, _, b, _ in records],
                                         'height_lowers': lower_heights,
                                         'omission_upper': str(omission_upper)
                                            if omission_upper is not None else None}
    return {'keep': best is not None and best >= 1, 'reason': 'component_mass',
            'upper': None if best is None else str(best), 'maximizing_relaxation': best_case}

