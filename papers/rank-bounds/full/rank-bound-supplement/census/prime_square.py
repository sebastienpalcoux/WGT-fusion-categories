"""Uniform actual/effective prime-square trace obstructions, p-free version.

The mathematical scope is the non-pointed based-simple rotation-admissible
class of the manuscript (prime-square and new:prime-free-square-trace). No rank bound or saved degree list is
used. A returned certificate is a necessary-condition failure, not a claim
that a retained degree list admits fusion rules or modular data.

prime_square_failure(degrees, q, M) validates the dimension/denominator
identities and returns a JSON-serializable certificate or None.
inspect_prime_class is the explicit local mathematical helper: unlike the
public checker, it does not assert a complete dimension identity, permitting
small synthetic tests of the support and block-capacity calculations.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, isqrt


@lru_cache(maxsize=2048)
def _factor(n):
    result = []
    p = 2
    while p*p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            result.append((p, exponent))
        p = 3 if p == 2 else p+2
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def _prime(n):
    return n >= 2 and _factor(n) == ((n, 1),)


@lru_cache(maxsize=128)
def unequal_patterns(ell):
    """All degree-count patterns (a,b,k) for unequal defect-one ell-blocks.

Such a block has a labels of degree u and b labels of degree k*u.
For ell=2 no unequal block exists. Equal blocks are handled separately.
"""
    if not _prime(ell):
        raise ValueError('the block prime must be prime')
    if ell == 2:
        return ()
    divisors = []
    for a in range(1, isqrt(ell-1)+1):
        if (ell-1) % a == 0:
            divisors.append(a)
            if a*a != ell-1:
                divisors.append((ell-1)//a)
    patterns = set()
    for a in divisors:
        for b in ((ell-1)//2, ell-1):
            square = (ell-1)*a//b
            k = isqrt(square)
            if k >= 2 and k*k == square:
                patterns.add((a, b, k))
    return tuple(sorted(patterns))


def partner_capacity(available, options):
    """Maximum target-label yield using at most available partner labels.

Each option is (partner_labels_consumed, target_labels_covered). This is a
small exact knapsack calculation, not a search for compatible full blocks.
"""
    if available < 0 or any(need <= 0 or gain <= 0 for need, gain in options):
        raise ValueError('invalid partner capacity problem')
    best = [0]*(available+1)
    for capacity in range(1, available+1):
        best[capacity] = best[capacity-1]
        for need, gain in options:
            if need <= capacity:
                best[capacity] = max(best[capacity], best[capacity-need]+gain)
    return best[available]


def unequal_capacity(counts, ell, degree):
    """An upper bound on labels of this degree coverable by unequal blocks.

Partner classes are optimized independently. Ignoring their obligations to
other classes enlarges the possibilities, hence cannot cause false failure.
"""
    partners = defaultdict(set)
    for a, b, k in unequal_patterns(ell):
        if counts.get(k*degree, 0):
            partners[k*degree].add((b, a))
        if degree % k == 0 and counts.get(degree//k, 0):
            partners[degree//k].add((a, b))
    details = []
    total = 0
    for other, options in sorted(partners.items()):
        choices = sorted(options)
        capacity = partner_capacity(counts[other], choices)
        total += capacity
        details.append({'degree': other, 'available': counts[other],
                        'options': [list(t) for t in choices],
                        'capacity': capacity})
    return total, details


def inspect_prime_class(counts, q, M, p):
    """Inspect one denominator-p class; no complete dimension check here.

Counts are nonnegative multiplicities of positive integral degrees dividing
M. This helper is useful for mathematical unit tests and independently
validated complete profiles. It returns None conservatively when q is
composite or when the proposed removed prime equals p.
"""
    if not _prime(p) or M % p:
        return None
    if q != 1 and (not _prime(q) or q == p):
        return None
    n = M//p
    multiplicity = counts.get(n, 0)
    if n <= 1 or not multiplicity:
        return None
    # Principality gives v(Y)=0 modulo the unique p-prime whenever
    # p divides d_Y. Its absolute trace is then a multiple of p, and
    # the strict norm bound |trace|<p forces zero. This holds at EVERY
    # p-height, so the target is the full prime-to-p part of M.
    target = M
    while target % p == 0:
        target //= p
    eligible = sorted(d for d, copies in counts.items()
                      if copies and d > 1 and target % d == 0)
    common = gcd(*eligible)  # gcd of the empty set is zero.
    support = sorted({s for d in eligible for s, _ in _factor(d)})
    violations = []
    if common != 1:
        violations.append('trace_support_gcd')
    if len(support) < 3:
        violations.append('trace_support_prime_count')
    if p % 2 and n % p == 0:
        violations.append('prime_square_height')
    if not violations:
        return None

    forced = None
    if q != 1:
        capacity, partners = unequal_capacity(counts, q, n)
        if multiplicity <= capacity:
            return None
        forced = {'prime': q, 'unequal_capacity': capacity,
                  'partner_capacities': partners,
                  'reason': 'multiplicity_exceeds_unequal_capacity'}
    result = {'criterion': 'prime_square_trace',
              'criterion_version': 'prime-to-p-support-v2',
              'kind': 'actual' if q == 1 else 'effective',
              'q': q, 'p': p, 'M': M, 'degree': n,
              'multiplicity': multiplicity,
              'trace_degree_bound': target,
              'trace_support_prime_to': p,
              'eligible_degrees': eligible, 'eligible_gcd': common,
              'eligible_primes': support, 'violations': violations}
    if forced is not None:
        result['forced_equal_block'] = forced
    return result


def possible_prime_packet(p, multiplicity, odd):
    """Exact necessary size and parity of a prime denominator's full class."""
    j, remainder = divmod(multiplicity, p)
    if remainder or not 0 < j < p-1 or (p-1) % j:
        return False
    if p == 5 and j != 1 or p == 7 and j > 2:
        return False
    if odd:
        two_part = (p-1) & -(p-1)
        return j % two_part == 0 and (p-1)//two_part >= 5
    return True


def prime_square_failure(degrees, q, M):
    """Return a structured actual/effective prime-square failure, or None.

Invalid dimension/denominator input raises ValueError and is never counted
as a mathematical exclusion. The support obstruction uses the based-simple
coprime-vanishing theorem; this API is not a test for arbitrary fusion rings.
"""
    if type(q) is not int or q < 1 or type(M) is not int or M < 1:
        raise ValueError('positive integral q and M are required')
    if any(exponent != 1 for _, exponent in _factor(q)):
        raise ValueError('q must be squarefree')
    if not isinstance(degrees, (list, tuple)) or not degrees:
        raise ValueError('a nonempty degree list is required')
    if any(type(d) is not int or d < 1 or M % d for d in degrees):
        raise ValueError('every degree must be a positive integral divisor of M')
    counts = Counter(degrees)
    if counts[1] != 1 or sum(d*d for d in degrees) != q*M*M:
        raise ValueError('the unique-unit or dimension identity fails')
    if q != 1 and not _prime(q):
        return None
    for p, _ in _factor(M):
        if q == 1 and p >= 5 and M > p and counts[M // p]:
            multiplicity = counts[M // p]
            if not possible_prime_packet(p, multiplicity, odd=M % 2):
                return {"kind": "actual", "p": p, "degree": M // p,
                        "multiplicity": multiplicity, "violations": ["prime_packet_multiplicity"]}
        failure = inspect_prime_class(counts, q, M, p)
        if failure is not None:
            return failure
    return None
