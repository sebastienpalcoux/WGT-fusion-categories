"""Safe prime-cut mass bounds for ordered nonunit denominator prefixes.

No factorization, smooth-number enumeration, or completion search is done.
The allowed support may be a superset of the final support. Each tested cut
must meet required_primes on both sides, so it remains nontrivial in every
completion. Supplied support entries are distinct primes by caller contract.
"""
from fractions import Fraction
from heapq import heappop, heappush


def _crossing_candidates(lower, excluded, left, right, count):
    """Least crossing integers, permitting primes outside the allowed support.

    This enlarged candidate set gives lower bounds on actual denominators.
    A number crosses exactly when it is divisible by some p*q with p on
    one side and q on the other. Merge these arithmetic progressions.
    """
    steps = sorted({p*q for p in left for q in right})
    heap = [(((lower+step-1)//step)*step, step) for step in steps]
    # The steps are sorted, but their first multiples need not be.
    from heapq import heapify
    heapify(heap)
    found, previous = [], None
    while len(found) < count:
        value, step = heappop(heap)
        heappush(heap, (value+step, step))
        if value != previous and value not in excluded:
            found.append(value)
        previous = value
    return found


def prime_cut_prefix_bound(prefix, remaining, primes, next_lower, *,
                           delta=1, required_primes=None, max_primes=12):
    """Return an exact upper bound on the future NONUNIT reciprocal mass.

    prefix contains nonunit denominators only; remaining is their future
    label count; every future denominator is at least next_lower. delta=2
    requires the odd-dimensional dual-pair theorem. required_primes are
    proved present in the final M; by default the whole support is fixed.

    The vacuum is excluded from the count and mass. A caller whose residual
    includes it must add a separately certified vacuum mass upper bound.
    Resource limits weaken to the ordinary mass bound, never to exclusion.
    """
    primes = tuple(primes)
    required = set(primes if required_primes is None else required_primes)
    if (primes != tuple(sorted(set(primes))) or any(p < 2 for p in primes)
            or not required.issubset(primes)):
        raise ValueError('invalid supplied prime support')
    if (type(remaining) is not int or remaining < 0 or
            type(next_lower) is not int or next_lower < 1 or delta not in (1, 2)):
        raise ValueError('invalid count, lower bound, or multiplicity')
    if any(type(b) is not int or b < 1 for b in prefix):
        raise ValueError('prefix denominators must be positive integers')
    ordinary = Fraction(remaining, next_lower*next_lower)
    answer = {'status': 'checked', 'upper_bound': ordinary,
              'certificate': None, 'cuts_tested': 0}
    if len(primes) > max_primes:
        answer['status'] = 'skipped'
        return answer
    whole = (1 << len(primes))-1
    observed = {b: sum(1 << j for j, p in enumerate(primes) if b % p == 0)
                for b in prefix}
    for mask in range(1, whole, 2):
        left = tuple(p for j, p in enumerate(primes) if mask >> j & 1)
        right = tuple(p for j, p in enumerate(primes) if not (mask >> j & 1))
        if not (required.intersection(left) and required.intersection(right)):
            continue
        answer['cuts_tested'] += 1
        crossing = {b for b, seen in observed.items()
                    if seen & mask and seen & (whole ^ mask)}
        deficit = max(0, 3-len(crossing))
        if not deficit:
            continue
        cert = {'left': list(left), 'right': list(right),
                'observed_crossing': sorted(crossing), 'deficit': deficit,
                'reserved_labels': delta*deficit}
        if remaining < delta*deficit:
            answer.update(status='impossible', upper_bound=None, certificate=cert)
            return answer
        candidates = _crossing_candidates(next_lower, crossing, left, right, deficit)
        bound = Fraction(remaining-delta*deficit, next_lower*next_lower)
        bound += delta*sum((Fraction(1, b*b) for b in candidates), Fraction())
        cert['new_value_lower_bounds'] = candidates
        if bound < answer['upper_bound']:
            answer.update(upper_bound=bound, certificate=cert)
    return answer


def prime_cut_prefix_failure(prefix, remaining, primes, next_lower, residual, *,
                             vacuum_upper=Fraction(), delta=1,
                             required_primes=None, max_primes=12):
    """An exclusion certificate, or None. residual and vacuum_upper are exact.

    residual is the required future mass, including the vacuum if and only
    if vacuum_upper is supplied. No equality endpoint is excluded here.
    """
    if not isinstance(residual, (int, Fraction)) or not isinstance(vacuum_upper, (int, Fraction)):
        raise TypeError('exact integers or fractions are required')
    if residual < 0 or vacuum_upper < 0:
        raise ValueError('masses must be nonnegative')
    result = prime_cut_prefix_bound(prefix, remaining, primes, next_lower,
        delta=delta, required_primes=required_primes, max_primes=max_primes)
    if result['status'] == 'impossible':
        return result
    if residual > result['upper_bound'] + vacuum_upper:
        return result
    return None
