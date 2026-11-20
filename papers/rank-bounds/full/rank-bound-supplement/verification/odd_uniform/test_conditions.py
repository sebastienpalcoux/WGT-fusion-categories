"""Small independent checks of the uniform odd reservation conditions."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd, lcm, prod
import json
import unittest
from unittest.mock import patch

import component
import coupled_weyl
import global_weyl
import residue
import check


def brute_entries(parts, parities):
    """Enumerate integer numerators and test every component union directly."""
    b = prod(parts)
    options = []
    for u, parity in zip(parts, parities):
        options.append([Q(a, u) for a in range(1-u, u)
                        if (u == 1 and a == 0 or u > 1 and (a*(b//u)-1) % u == 0)
                        and (parity is None or a % 2 == parity)])
    answers = []
    for entries in product(*options):
        if sum(entries, Q()) != Q(1, b):
            continue
        if any(abs(sum((entries[j] for j in range(len(parts)) if mask >> j & 1), Q())) >= 1
               for mask in range(1, 1 << len(parts))):
            continue
        answers.append(entries)
    if not answers:
        return None
    return (min(sum(x*x for x in e) for e in answers),
            tuple(min(e[j]*e[j] for e in answers) for j in range(len(parts))))


def brute_support_feasible(p, support, upper, forced):
    """Enumerate actual bounded divisors, then all support subfamilies."""
    others = tuple(q for q in support if q != p)
    candidates = {}
    for exponents in product(*(range(upper[q]+1) for q in others)):
        mask = sum(1 << j for j, e in enumerate(exponents) if e)
        if mask.bit_count() < 2:
            continue
        d = prod(q**e for q, e in zip(others, exponents))
        least = min(d % p, (-d) % p)
        twice = min((2*d) % p, (-2*d) % p)
        even, odd = candidates.setdefault(mask, (p*p, p*p))
        even = min(even, 4*least*least)
        if twice % 2 and all(mask >> others.index(q) & 1 for q in forced):
            odd = min(odd, twice*twice)
        candidates[mask] = even, odd
    for odd_mask, (_, odd_cost) in candidates.items():
        for count in range(len(candidates)+1):
            for even_masks in combinations(candidates, count):
                masks = (odd_mask,)+even_masks
                if any(not (a & b) for a, b in combinations(masks, 2)):
                    continue
                intersection = (1 << len(others))-1
                for mask in masks:
                    intersection &= mask
                if intersection:
                    continue
                if odd_cost+sum(candidates[m][0] for m in even_masks) <= p-2:
                    return True
    return False


class UniformOddConditions(unittest.TestCase):
    def test_receipt_structured_outcomes_roundtrip_exactly(self):
        outcome = {'maximizing_relaxation': {'height_lowers': {3: 6, 5: 3}},
                   'marked': ((11, 2), (19, 2))}
        canonical = check.canonical_result(outcome)
        self.assertEqual(canonical['maximizing_relaxation']['height_lowers'],
                         {'3': 6, '5': 3})
        self.assertEqual(canonical['marked'], [[11, 2], [19, 2]])
        self.assertEqual(canonical, json.loads(json.dumps(canonical)))

    def test_coupled_entries_against_all_lifts_and_unions(self):
        for parts in ((1, 3, 5), (3, 5, 7), (9, 5, 7), (3, 25, 11), (27, 5, 7)):
            for parities in product((None, 0, 1), repeat=len(parts)):
                with self.subTest(parts=parts, parities=parities):
                    expected = brute_entries(parts, parities)
                    got = coupled_weyl.coupled_weyl_entry_floors(parts, parities)
                    self.assertEqual(got is None, expected is None)
                    if expected:
                        self.assertEqual(got['joint_squared_floor'], expected[0])
                        self.assertEqual(got['entry_squared_floors'], expected[1])

    def test_support_budget_against_bounded_divisors(self):
        cases = ((11, (3, 5, 7, 11)), (19, (3, 5, 11, 19)),
                 (23, (3, 5, 11, 23)), (31, (3, 5, 7, 31)))
        for p, support in cases:
            for height in (1, 2):
                upper = {q: height if q in (3, 5) else 1 for q in support if q != p}
                for forced in ((), (3,)):
                    with self.subTest(p=p, upper=upper, forced=forced):
                        expected = brute_support_feasible(p, support, upper, forced)
                        got = residue.minimum_support_budget(p, support,
                                tuple(sorted(upper.items())), forced)
                        self.assertEqual(got <= p-2, expected)

    def test_positive_powers_include_full_period(self):
        for q, p in ((3, 11), (5, 19), (7, 23)):
            self.assertEqual(residue.positive_power_residues(q, p),
                             frozenset(pow(q, e, p) for e in range(1, p)))
            for h in (1, 2, 4):
                self.assertEqual(residue.positive_power_residues(q, p, h),
                                 frozenset(pow(q, e, p) for e in range(1, h+1)))

    def test_higher_marker_is_not_an_upper_height(self):
        row = [[3, 5, 11, 19], [[0, 0, 0], [1, 22, 20], [1, 38, 36]], [], 44]
        higher = [row[0], [[2, 4, 4], *row[1][1:]], row[2], row[3]]
        self.assertEqual(residue.height_lower_bounds(row), residue.height_lower_bounds(higher))
        self.assertEqual(residue.support_budget_failure(row), residue.support_budget_failure(higher))

    def test_edge_prime_height_one_is_admitted(self):
        support, marked = (3, 5, 11, 13), ((11, 4), (13, 4))
        def brute(height5):
            for a in range(2, 2+lcm(10, 12)):
                M = 3**a*5**height5*11*13
                values = tuple(j*(M//p) % p for p, j in marked)
                if values in ((1, 1), (10, 12)):
                    return True
            return False
        self.assertTrue(brute(1))
        self.assertFalse(brute(2))
        self.assertIs(global_weyl.possible_targets(support, marked, 5, 1, 1), True)
        self.assertIs(global_weyl.possible_targets(support, marked, 5, 2, 2), False)
        row = [list(support), [[0, 0, 0], [1, 44, 40], [1, 52, 48]], [[3, 5]], 42]
        failure, heights = global_weyl.global_constraints(row)
        self.assertIsNone(failure)
        self.assertNotIn(5, heights)

    def test_resource_caps_retain_unknown(self):
        row = [[3, 5, 11, 19], [[0, 0, 0], [1, 22, 20], [1, 38, 36]], [], 44]
        global_weyl.possible_targets.cache_clear()
        try:
            with patch.object(global_weyl, 'CAP', 0):
                failure, heights = global_weyl.global_constraints(row)
                self.assertIsNone(failure)
                self.assertEqual(heights, {})
        finally:
            global_weyl.possible_targets.cache_clear()
        self.assertIsNone(residue.minimum_support_budget(263, (3, 5, 7, 263), (), ()))
        self.assertIsNone(residue.fixed_pair_omission_failure(263, (3, 5, 7, 263), {}, {3}))

    def test_invalid_primitive_multiplicity_has_an_exact_rejection(self):
        row = [[3, 5, 7, 11], [[0, 0, 0], [0, 0, 0], [1, 132, 110]], [], 42]
        failure, heights = global_weyl.global_constraints(row)
        self.assertEqual(failure, {'reason': 'valuation_one_primitive_multiplicity',
                                  'prime': 11, 'primitive_multiplicity': 11})
        self.assertEqual(heights, {})
        self.assertIs(global_weyl.possible_targets((3, 5, 7, 11), ((11, 11),)), False)

    def test_outside_scope_is_retained(self):
        rows = (
            [[5, 7, 11], [[0, 0, 0]]*3, [], 0],
            [[3, 5, 7], [[0, 0, 0]]*2, [], 0],
            [[3, 5, 7], [[0, 0, 0]]*2, [[3, 7]], 42],
            [[3, 5, 7], [[2, 4, 4], [0, 0, 0]], [], 42],
        )
        for row in rows:
            with self.subTest(row=row):
                self.assertEqual(component.component_envelope(row, 209),
                                 {'keep': True, 'reason': 'outside_scope'})

    def test_larger_rank_does_not_inherit_a_fixed_rank_cutoff(self):
        row = [[3, 5, 7], [[0, 0, 0], [0, 0, 0]], [], 66]
        verdict = component.component_envelope(row, 1009)
        self.assertTrue(verdict['keep'])
        self.assertGreaterEqual(Q(verdict['upper']), 1)

    def test_three_prime_prime_class_is_explicitly_impossible(self):
        row = [[3, 5, 11], [[0, 0, 0], [1, 22, 20]], [], 42]
        self.assertEqual(component.component_envelope(row, 209)['reason'],
                         'prime_denominator_with_three_primes')


if __name__ == '__main__':
    unittest.main(verbosity=2)
