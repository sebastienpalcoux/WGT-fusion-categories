"""Reproduce the support-by-support dimension comparison in Theorem 17.1."""

from heapq import heappop, heappush
from math import prod


def covered(primes, exponents):
    """Test the sufficient exponent conditions used in this comparison."""
    ordered = sorted(exponents)
    if len(primes) == 3:
        # The first three bounds use increasing prime order. The box (7,5,4)
        # can be applied after any permutation of the primes.
        a, b, c = exponents
        in_box = all(e <= limit for e, limit in zip(ordered, (4, 5, 7)))
        return a <= 5 or b <= 4 or c <= 3 or in_box

    if len(primes) == 4:
        d, c, b, a = ordered
        if b <= 3:
            return True
        if d != 2 or c > 3:
            return False
        dyadic = primes[0] == 2 and exponents[0] >= 4
        return a // 2 + b // 2 <= 5 or (dyadic and a <= 9 and b <= 5)

    if len(primes) == 5:
        least, _, _, second_largest, largest = ordered
        return largest <= 3 or (largest == 4 and second_largest <= 3 and least == 2)

    # At six primes the product of the initial powers already gives the bound.
    return False


def first(primes, center):
    """Return the least product outside the covered families on this support."""
    minimum_exponent = 1 if center else 2
    initial = (minimum_exponent,) * len(primes)
    queue = [(prod(primes) ** minimum_exponent, initial)]
    seen = {initial}
    while queue:
        dimension, exponents = heappop(queue)
        # For a center, enumerate N and test the exponents of N squared.
        tested = tuple(2 * e for e in exponents) if center else exponents
        if not covered(primes, tested):
            return dimension
        for index, prime in enumerate(primes):
            neighbor = list(exponents)
            neighbor[index] += 1
            successor = tuple(neighbor)
            if successor not in seen:
                seen.add(successor)
                heappush(queue, (dimension * prime, successor))


def main():
    primes = (2, 3, 5, 7, 11, 13, 17)
    # Columns in the proof: Q(D), Q(D) odd, N, N odd.
    expected_rows = {
        3: (2**6 * 3**6 * 5**4, 3**6 * 5**6 * 7**4, 5400, 165375),
        4: (6300**2, 121275**2, 6300, 121275),
        5: (2**3 * 2310**2, 3**3 * 15015**2, 9240, 135135),
        6: (30030**2, 255255**2, 30030, 255255),
    }
    rows = []
    for size, expected in expected_rows.items():
        row = []
        for center in (False, True):
            even = first(primes[:size], center)
            odd = first(primes[1:size + 1], center)
            row.extend((min(even, odd), odd))
        if tuple(row) != expected:
            raise SystemExit(f"{size} repeated primes: expected {expected}, got {tuple(row)}")
        rows.append(row)

    braided, braided_odd, general, general_odd = map(min, zip(*rows))
    print(braided, braided_odd)
    print(general, general_odd)


if __name__ == "__main__":
    main()
