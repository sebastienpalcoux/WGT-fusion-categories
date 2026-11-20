"""Small mathematical tests; no census or saved survivor lists are imported."""
from collections import Counter
from itertools import product
from math import isqrt
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prime_square import (inspect_prime_class, partner_capacity,
                          prime_square_failure, unequal_capacity,
                          unequal_patterns)


def naive_patterns(ell):
    if ell == 2:
        return ()
    answer = []
    for a in range(1, ell):
        for b in range(1, ell):
            if (ell-1) % a or (ell-1) % b or b % ((ell-1)//2):
                continue
            for k in range(2, ell):
                if b*k*k == (ell-1)*a:
                    answer.append((a, b, k))
    return tuple(sorted(answer))


def naive_capacity(available, options):
    """Exhaustively list option multiplicities, independently of the DP."""
    return max(sum(t*gain for t, (_, gain) in zip(mults, options))
               for mults in product(*(range(available//need+1)
                                      for need, _ in options))
               if sum(t*need for t, (need, _) in zip(mults, options)) <= available)


class PatternTests(unittest.TestCase):
    def test_known_patterns(self):
        self.assertEqual(unequal_patterns(2), ())
        self.assertEqual(unequal_patterns(3), ((2, 1, 2),))
        self.assertEqual(unequal_patterns(5), ((2, 2, 2), (4, 4, 2)))
        self.assertEqual(unequal_patterns(7), ((2, 3, 2),))

    def test_independent_equation_enumeration(self):
        for ell in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            self.assertEqual(unequal_patterns(ell), naive_patterns(ell))

    def test_capacity_against_naive_enumeration(self):
        for options in ([(2, 1)], [(2, 2), (4, 4)], [(2, 3)],
                        [(2, 1), (3, 2), (5, 4)], [(2, 8), (4, 16)]):
            for available in range(13):
                self.assertEqual(partner_capacity(available, options),
                                 naive_capacity(available, options))

    def test_small_prime_capacity_formulas(self):
        # The target is 30, with potential lower/higher partners 15 and60.
        # These are local block-count problems, not full dimension profiles.
        for lower in range(9):
            for upper in range(9):
                counts = {15: lower, 30: 100, 60: upper}
                self.assertEqual(unequal_capacity(counts, 3, 30)[0],
                                 lower//2 + 2*upper)
                self.assertEqual(unequal_capacity(counts, 5, 30)[0],
                                 2*(lower//2) + 2*(upper//2))
                self.assertEqual(unequal_capacity(counts, 7, 30)[0],
                                 3*(lower//2) + 2*(upper//3))


class LocalCriterionTests(unittest.TestCase):
    """Synthetic local data are intentional; use the documented helper."""
    def test_actual_gcd_failure(self):
        failure = inspect_prime_class({1: 1, 6: 1, 12: 1, 18: 1, 180: 1},
                                      1, 900, 5)
        self.assertEqual(failure['kind'], 'actual')
        self.assertEqual(failure['trace_degree_bound'], 36)
        self.assertEqual(failure['eligible_gcd'], 6)

    def test_empty_support_is_failure(self):
        failure = inspect_prime_class({1: 1, 25: 2}, 2, 125, 5)
        self.assertEqual(failure['eligible_degrees'], [])
        self.assertEqual(failure['eligible_gcd'], 0)

    def test_prime_count_can_fail_with_gcd_one(self):
        failure = inspect_prime_class({1: 1, 2: 1, 3: 1, 6: 1, 120: 3},
                                      3, 600, 5)
        self.assertEqual(failure['eligible_gcd'], 1)
        self.assertEqual(failure['violations'], ['trace_support_prime_count', 'prime_square_height'])

    def test_possible_unequal_blocks_are_retained(self):
        base = {1: 1, 6: 1, 12: 1, 18: 1, 90: 6, 180: 3}
        self.assertIsNone(inspect_prime_class(base, 3, 900, 5))
        # Six lower labels can cover at most three target labels.
        more = dict(base, **{})
        more[180] = 4
        failure = inspect_prime_class(more, 3, 900, 5)
        self.assertEqual(failure['forced_equal_block']['unequal_capacity'], 3)

    def test_q5_and_q7_capacity_boundary(self):
        for ell, cover in ((5, 2), (7, 3)):
            counts = {1: 1, 6: 1, 12: 1, 18: 1, 90: 2, 180: cover}
            self.assertIsNone(inspect_prime_class(counts, ell, 900, 5))
            # At ell=5 the proposed removed prime equals p and must be skipped.
            if ell != 5:
                counts[180] += 1
                failure = inspect_prime_class(counts, ell, 900, 5)
                self.assertEqual(failure['forced_equal_block']['unequal_capacity'], cover)
        # A distinct target prime exercises q=5 rather than its skip branch.
        counts = {1: 1, 6: 1, 12: 1, 18: 1, 126: 2, 252: 2}
        self.assertIsNone(inspect_prime_class(counts, 5, 1764, 7))
        counts[252] = 3
        self.assertEqual(inspect_prime_class(counts, 5, 1764, 7)
                         ['forced_equal_block']['unequal_capacity'], 2)

    def test_prime_two_actual_and_effective(self):
        counts = {1: 1, 6: 1, 18: 1, 36: 3}
        for q in (1, 3):
            failure = inspect_prime_class(counts, q, 72, 2)
            self.assertEqual(failure['p'], 2)
            self.assertEqual(failure['trace_degree_bound'], 9)
            self.assertEqual(failure['eligible_degrees'], [])

    def test_valuation_one_uses_n_not_n_over_p(self):
        counts = {1: 1, 6: 1, 14: 1, 21: 1, 42: 3}
        self.assertIsNone(inspect_prime_class(counts, 3, 210, 5))
        bad = {1: 1, 6: 1, 42: 3}
        failure = inspect_prime_class(bad, 3, 210, 5)
        self.assertEqual(failure['trace_degree_bound'], 42)

    def test_high_valuation_removes_p_divisible_trace_degrees(self):
        # The older M/p^2 bound admitted 6,10,15. Principality removes
        # 10 and15 at every height, leaving a genuine gcd obstruction.
        counts = {1: 1, 6: 1, 10: 1, 15: 1, 600: 3}
        failure = inspect_prime_class(counts, 3, 3000, 5)
        self.assertEqual(failure['trace_degree_bound'], 24)
        self.assertEqual(failure['eligible_degrees'], [6])
        self.assertEqual(failure['eligible_gcd'], 6)

    def test_prime_free_triangle_still_requires_valuation_one(self):
        for height in (1, 2, 3, 7):
            M = 252 * 5**height
            for q, copies in ((1, 1), (2, 2), (3, 3)):
                counts = {1: 1, 6: 1, 14: 1, 21: 1, M//5: copies}
                failure = inspect_prime_class(counts, q, M, 5)
                if height == 1:
                    self.assertIsNone(failure)
                else:
                    self.assertEqual(failure['violations'], ['prime_square_height'])

    def test_composite_kernel_and_equal_primes_are_retained(self):
        counts = {1: 1, 6: 1, 12: 1, 18: 1, 180: 5}
        self.assertIsNone(inspect_prime_class(counts, 6, 900, 5))
        self.assertIsNone(inspect_prime_class(counts, 5, 900, 5))

    def test_nonoccurring_or_unit_class_is_retained(self):
        self.assertIsNone(inspect_prime_class({1: 1, 6: 1}, 1, 900, 5))
        self.assertIsNone(inspect_prime_class({1: 1}, 1, 5, 5))


class PacketSizeTests(unittest.TestCase):
    def test_ordinary_and_odd_packets(self):
        from prime_square import possible_prime_packet
        for p, expected, expected_odd in [(5, [1], []), (7, [1, 2], []),
                                          (11, [1, 2, 5], [2]), (13, [1, 2, 3, 4, 6], [])]:
            for odd, j_values in ((False, expected), (True, expected_odd)):
                valid = [m for m in range(1, p*p) if possible_prime_packet(p, m, odd)]
                self.assertEqual(valid, [p*j for j in j_values])


class PublicApiTests(unittest.TestCase):
    def test_complete_arithmetic_profiles(self):
        square = [1, 2, 2, 3, 3, 3]
        self.assertEqual(sum(d*d for d in square), 6**2)
        self.assertEqual(prime_square_failure(square, 1, 6)['kind'], 'actual')
        twice = square + [6]
        self.assertEqual(sum(d*d for d in twice), 2*6**2)
        self.assertEqual(prime_square_failure(twice, 2, 6)['kind'], 'effective')
        composite = square + [6]*5
        self.assertEqual(sum(d*d for d in composite), 6*6**2)
        self.assertIsNone(prime_square_failure(composite, 6, 6))
        self.assertIsNone(prime_square_failure([1], 1, 1))

    def test_invalid_data_raise_instead_of_excluding(self):
        for degrees, q, M in (([1, 2], 1, 6), ([1, 1], 1, 1),
                               ([1], 4, 1), ([1, 5], 1, 6),
                               ([1, True], 1, 1)):
            with self.assertRaises(ValueError):
                prime_square_failure(degrees, q, M)


if __name__ == '__main__':
    unittest.main()
