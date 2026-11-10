"""Reciprocal-tuple counts for ranks 6, 9 and 11 in Corollary 18.2 of the 229-page coauthor edition."""
from fractions import Fraction
from math import lcm


def dimensions(remaining, weight, minimum, odd, remainder=Fraction(1), divisor=1):
    """Yield the last denominator D of each admissible nondecreasing tuple."""
    if remaining == 0:
        dim = remainder.denominator
        if (remainder.numerator == 1 and dim >= minimum
                and dim % divisor == 0 and (not odd or dim % 2 == 1)):
            yield dim
        return

    low = max(minimum, weight // remainder + 1)
    high = (remaining * weight + 1) // remainder
    if odd and low % 2 == 0:
        low += 1
    step = 2 if odd else 1
    for x in range(low, high + 1, step):
        yield from dimensions(remaining - 1, weight, x, odd,
                              remainder - Fraction(weight, x), lcm(divisor, x))


if __name__ == "__main__":
    high = sorted({d for d in dimensions(5, 1, 3, False) if d >= 5400})
    assert high == [6240, 7140, 7320, 8268, 8364, 12324, 14520, 24492]
    print("Rank 6 dimensions at least 5400:", *high)
    for pairs, expected in ((4, (3, 105)), (5, (30, 22155))):
        dims = list(dimensions(pairs, 2, 7, True))
        assert (len(dims), max(dims)) == expected
        print(f"Rank {2 * pairs + 1}: {len(dims)} tuples, maximum {max(dims)}")
