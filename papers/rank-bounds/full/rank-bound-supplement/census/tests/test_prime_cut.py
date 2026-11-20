"""Exact small-completion audits of terminal and prefix prime-cut tests."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prime_cut import prime_cut_diagnostics, prime_cut_failure, prefix_count_failure
from prime_cut_prefix import prime_cut_prefix_bound, prime_cut_prefix_failure


def degrees(M, denominators):
    return sorted([1]+[M//b for b in denominators])


def every_cut_has_three(values, primes):
    for size in range(1, len(primes)):
        for left in combinations(primes, size):
            right = set(primes)-set(left)
            crossing = {b for b in values if any(b % p == 0 for p in left)
                        and any(b % p == 0 for p in right)}
            if len(crossing) < 3:
                return False
    return True


class TerminalCutTests(unittest.TestCase):
    def test_distinct_classes_gcd_and_vacuum(self):
        support = (2, 3, 5, 7)
        good = (6, 10, 14, 15, 21, 35)
        self.assertIsNone(prime_cut_failure(degrees(210, good), 210, support))
        bad_height = prime_cut_diagnostics(degrees(420, good), 420, support)
        self.assertTrue(bad_height['failures'])
        self.assertTrue(all(not f['fails_count'] for f in bad_height['failures']))
        bad_count = prime_cut_diagnostics(degrees(30, (6, 10, 15)), 30, (2, 3, 5))
        self.assertTrue(all(f['fails_count'] and not f['fails_gcd'] for f in bad_count['failures']))
        repeated = prime_cut_diagnostics(degrees(30, (6, 6, 10, 10, 15, 15)), 30, (2, 3, 5))
        self.assertEqual(repeated, bad_count)
        self.assertTrue(all(f['distinct_count'] == 2 for f in bad_count['failures']))
        skipped = prime_cut_diagnostics(degrees(210, good), 210, support, max_primes=3)
        self.assertEqual(skipped['status'], 'skipped')
        self.assertEqual(skipped['failures'], [])

    def test_count_prefixes_of_all_small_admissible_completions(self):
        for length in range(6):
            for full in combinations_with_replacement((2, 3, 5, 6, 10, 15, 30), length):
                if prefix_count_failure(full, 0, (2, 3, 5)) is None:
                    for end in range(length+1):
                        self.assertIsNone(prefix_count_failure(full[:end], length-end, (2, 3, 5)))


class PrefixMassTests(unittest.TestCase):
    def test_all_small_completions_including_open_odd_classes(self):
        values = (6, 10, 12, 14, 15, 20, 21, 28, 30, 35, 42)
        actual = (2, 3, 5, 7)
        projections = 0
        for size in (6, 7):
            for base in combinations(values, size):
                if not every_cut_has_three(base, actual):
                    continue
                for delta in (1, 2):
                    full = tuple(b for b in base for _ in range(delta))
                    for support in (actual, actual+(11,)):
                        for end in range(len(full)+1):
                            prefix, suffix = full[:end], full[end:]
                            lower = suffix[0] if suffix else full[-1]
                            answer = prime_cut_prefix_bound(prefix, len(suffix), support,
                                lower, delta=delta, required_primes=actual)
                            mass = sum((F(1, b*b) for b in suffix), F())
                            self.assertNotEqual(answer['status'], 'impossible')
                            self.assertLessEqual(mass, answer['upper_bound'])
                            self.assertIsNone(prime_cut_prefix_failure(prefix, len(suffix),
                                support, lower, mass+F(1, 210**2), vacuum_upper=F(1, 210**2),
                                delta=delta, required_primes=actual))
                            projections += 1
        self.assertEqual(projections, 12110)

    def test_unknown_support_and_safe_skips(self):
        answer = prime_cut_prefix_bound([], 0, (2, 3), 6, required_primes=(2,))
        self.assertEqual(answer['cuts_tested'], 0)
        self.assertEqual(answer['status'], 'checked')
        self.assertEqual(prime_cut_prefix_bound([], 5, (3, 5), 15, delta=2)['status'], 'impossible')
        self.assertEqual(prime_cut_prefix_bound([], 1, (2, 3, 5, 7), 6, max_primes=3)['status'], 'skipped')

    def test_exact_mass_input_and_strict_endpoint(self):
        upper = prime_cut_prefix_bound([], 3, (3, 5), 15)['upper_bound']
        self.assertIsNone(prime_cut_prefix_failure([], 3, (3, 5), 15, upper))
        self.assertIsNotNone(prime_cut_prefix_failure([], 3, (3, 5), 15, upper+F(1, 10**9)))
        with self.assertRaises(TypeError):
            prime_cut_prefix_failure([], 3, (3, 5), 15, 0.5)


if __name__ == '__main__':
    unittest.main()
