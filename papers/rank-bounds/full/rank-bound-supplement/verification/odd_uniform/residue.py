"""Small, uniform residue relaxations of minimal prime-Weyl support.

No degree list or Gauss sum is enumerated. Unknown resource-limited cases
are retained. A retained reservation is not asserted to be realizable.
"""
from functools import lru_cache
from itertools import combinations


SUPPORT_LIMIT = 7
PRIME_LIMIT = 257


def least_residue(a, p):
    a %= p
    return min(a, p-a)


@lru_cache(None)
def positive_power_residues(q, p, maximum_exponent=None):
    """All q^e mod p for e>=1, or 1<=e<=the known valuation bound."""
    if q == p:
        raise ValueError('a proper p-free degree cannot contain p')
    out = set()
    value = q % p
    e = 1
    while maximum_exponent is None or e <= maximum_exponent:
        if value in out:
            break
        out.add(value)
        value = value*q % p
        e += 1
    return frozenset(out)


def possible_degree_residues(p, support, upper_valuations, omitted):
    """Conservative residues of nonunit degrees with prescribed omissions.

    Nonunit degrees use at least two distinct primes.  A supplied valuation
    is only an upper bound; absent valuations permit all positive powers.
    Consequently every actual eligible degree is represented, even though
    the output generally includes degrees that do not divide the actual M.
    """
    allowed = tuple(q for q in support if q != p and q not in omitted)
    residues = set()
    supports = []
    for size in range(2, len(allowed)+1):
        for subset in combinations(allowed, size):
            current = {1}
            for q in subset:
                powers = positive_power_residues(q, p, upper_valuations.get(q))
                current = {a*b % p for a in current for b in powers}
            residues.update(current)
            supports.append({'support': list(subset), 'residues': sorted(current)})
    return residues, supports


def fixed_pair_omission_failure(p, support, upper_valuations, omitted):
    """Exclude if a required proper, p-inertia-fixed omission has no residue.

    Hypotheses supplied by the caller: p is isolated, v_p(M)=1,
    primitive multiplicity j=2, and some proper p-free degree avoids every
    prime in `omitted` and has row field in another twist component.
    """
    if p > PRIME_LIMIT or len(support) > SUPPORT_LIMIT:
        return None
    residues, supports = possible_degree_residues(
        p, tuple(support), upper_valuations, frozenset(omitted))
    admissible = sorted(a for a in residues
                        if 4*least_residue(a, p)**2 <= p-7)
    if admissible:
        return None
    return {'prime': p, 'omitted': sorted(omitted),
            'upper_valuations': dict(sorted(upper_valuations.items())),
            'pair_square_budget': p-7,
            'possible_support_residues': supports}


@lru_cache(None)
def minimum_support_budget(p, support, upper_items, forced_odd_primes):
    """Relaxed least square cost of a pairwise-noncoprime gcd-one support.

    There is exactly one odd coefficient.  Every other coefficient is even.
    Equal prime-support masks are allowed to hide different actual degrees;
    retaining only a cheapest representative can only decrease the cost.
    We require selected supports to meet pairwise and have empty intersection.
    """
    if p > PRIME_LIMIT or len(support) > SUPPORT_LIMIT:
        return None
    upper = dict(upper_items)
    primes = tuple(q for q in support if q != p)
    costs = []
    for size in range(2, len(primes)+1):
        for subset in combinations(primes, size):
            residues = {1}
            for q in subset:
                residues = {a*b % p for a in residues for b in
                            positive_power_residues(q, p, upper.get(q))}
            mask = sum(1 << primes.index(q) for q in subset)
            even = min((4*least_residue(a, p)**2 for a in residues), default=p*p)
            odd = min((least_residue(2*a, p)**2 for a in residues
                       if least_residue(2*a, p) % 2), default=p*p)
            if not set(forced_odd_primes).issubset(subset):
                odd = p*p
            costs.append((mask, even, odd))
    even_candidates = [(mask, cost) for mask, cost, _ in costs if cost <= p-7]
    compatibility = []
    for mask, _ in even_candidates:
        compatibility.append(sum(1 << j for j, (other, _) in
                                 enumerate(even_candidates) if mask & other))

    @lru_cache(None)
    def finish(intersection, allowed):
        if not intersection:
            return 0
        best = p*p
        bits = allowed
        while bits:
            bit = bits & -bits
            bits -= bit
            j = bit.bit_length()-1
            mask, cost = even_candidates[j]
            smaller = intersection & mask
            # A redundant support can be omitted from any putative cover.
            if smaller == intersection:
                continue
            rest = finish(smaller, allowed & compatibility[j])
            best = min(best, cost+rest)
        return best

    best = p*p
    for mask, _, odd_cost in costs:
        if odd_cost > p-10:
            continue
        allowed = sum(1 << j for j, (other, _) in
                      enumerate(even_candidates) if mask & other)
        best = min(best, odd_cost+finish(mask, allowed))
    return best


def support_budget_failure(row):
    support, classes, edges, nines = row
    class_by_prime = dict(zip((q for q in support if q >= 5), classes))
    # Marker 2 denotes valuation >=2, not a finite upper bound.
    upper = tuple(sorted((q, 1) for q, c in class_by_prime.items() if c[0] == 1))
    touching = {q for edge in edges for q in edge}
    for p, (valuation, total, active) in class_by_prime.items():
        if (p in touching or valuation != 1 or total != 2*p
                or active != 2*(p-1)):
            continue
        forced = (3,) if 3 in support and nines >= 42 else ()
        lower = minimum_support_budget(p, tuple(support), upper, forced)
        if lower is not None and lower > p-2:
            return {'reason': 'proper_support_residue_budget', 'prime': p,
                    'available_pair_square_budget': p-2,
                    'relaxed_minimum_pair_square_budget': lower,
                    'forced_odd_degree_primes': list(forced)}
    return None


@lru_cache(None)
def _height_lowers(support, classes, edges, nines, cap):
    class_by_prime = dict(zip((q for q in support if q >= 5), classes))
    upper = {q: 1 for q, c in class_by_prime.items() if c[0] == 1}
    touching = {q for edge in edges for q in edge}
    result = {}
    for p, (valuation, total, active) in class_by_prime.items():
        if (p in touching or valuation != 1 or total != 2*p
                or active != 2*(p-1)):
            continue
        forced = (3,) if 3 in support and nines >= 42 else ()
        for q in support:
            if q == p or q in upper:
                continue
            lower = 1
            for h in range(1, cap+1):
                trial = dict(upper)
                trial[q] = h
                cost = minimum_support_budget(
                    p, support, tuple(sorted(trial.items())), forced)
                if cost is None or cost <= p-2:
                    break
                lower = h+1
            if lower > result.get(q, 1):
                result[q] = lower
    return tuple(sorted(result.items()))


def height_lower_bounds(row, cap=6):
    """Proved lower bounds only; hitting the cap is never a rejection."""
    support, classes, edges, nines = row
    return dict(_height_lowers(tuple(support), tuple(map(tuple, classes)),
                              tuple(map(tuple, edges)), nines, cap))


def reservation_failure(row):
    support, classes, edges, nines = row
    # The large nine-class forces every three-free field into the component
    # containing three.  No such deduction is made at smaller multiplicity.
    if 3 not in support or nines < 42:
        return None
    class_by_prime = dict(zip((q for q in support if q >= 5), classes))
    # Marker 2 denotes valuation >=2, not a finite upper bound.
    upper = {q: 1 for q, c in class_by_prime.items() if c[0] == 1}
    touching = {q for edge in edges for q in edge}
    for p, (valuation, total, active) in class_by_prime.items():
        if (p in touching or valuation != 1 or total != 2*p
                or active != 2*(p-1)):
            continue
        failure = fixed_pair_omission_failure(p, support, upper, {3})
        if failure:
            return {'reason': 'fixed_proper_omission_residue', **failure}
    return None

