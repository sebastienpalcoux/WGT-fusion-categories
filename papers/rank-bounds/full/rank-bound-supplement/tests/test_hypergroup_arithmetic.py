"""Exact checks of the probability example and the three small-degree arguments.

These short checks use only rational arithmetic and nonnegative integral
multiplicities. They can be run independently of the native generators.
"""
from fractions import Fraction
import unittest


def solutions(target, weights):
    """Nonnegative coefficients of all weighted partitions of target."""
    if len(weights) == 1:
        if target % weights[0] == 0:
            yield (target // weights[0],)
        return
    first, *rest = weights
    for n in range(target // first + 1):
        for tail in solutions(target - n * first, rest):
            yield (n, *tail)


class AlgebraicArithmetic(unittest.TestCase):
    def test_rank_two_probability_group(self):
        for degree in range(2, 17):
            k = Fraction(1, degree**2)
            nontrivial = -k
            self.assertEqual(nontrivial**2, k + (1-k)*nontrivial)
            self.assertEqual(degree**2*k, 1)
            self.assertEqual(degree*(1-k), degree-Fraction(1, degree))
            self.assertNotEqual((degree*(1-k)).denominator, 1)
            self.assertEqual(Fraction(1+degree**2, degree**2), 1+k)
            self.assertLess(-(degree**2-1), -2)  # The required theta + theta^{-1}.

    def test_first_small_degree_profile(self):
        for c66, c610, c612, c615 in solutions(35, (6, 10, 12, 15)):
            self.assertEqual((c66, c610, c612, c615), (0, 2, 0, 1))
            for _, c1012, _ in solutions(60 - 6*c610, (10, 12, 15)):
                self.assertEqual(c1012, 4)
                self.assertEqual(list(solutions(144 - 10*c1012, (12, 15))), [])
        self.assertEqual(len(list(solutions(35, (6, 10, 12, 15)))), 1)
        self.assertTrue(list(solutions(48, (10, 12, 15))))

    def test_second_small_degree_profile(self):
        self.assertEqual(list(solutions(99, (10, 18, 40, 45))), [(0, 3, 0, 1)])
        second = list(solutions(330, (18, 40, 45)))
        self.assertTrue(second)
        for _, c1840, _ in second:
            remainder = 400 - 18*c1840
            self.assertEqual(remainder, 292)
            self.assertLess(remainder, 324)
            self.assertFalse(list(solutions(remainder, (40, 45))))

    def test_third_small_degree_profile(self):
        first = list(solutions(287, (12, 15, 20, 144)))
        second = list(solutions(360, (12, 15, 20, 144)))
        third = list(solutions(240, (12, 15, 20, 144)))
        # Impose symmetry of the class totals, then enforce the individual
        # 12*12 and 12*15 decompositions. No constituent can have degree 250.
        match = [(A,B,C,E,V) for A,B,C,_ in first for B2,_,E,V in second
                 if B == B2 for C2,E2,_,_ in third if (C,E) == (C2,E2)]
        self.assertTrue(match)
        remain = [(A,B,C,E,V) for A,B,C,E,V in match
                  if any(a <= A and b <= B and c <= C
                         for a,b,c in solutions(143, (12, 15, 20)))]
        self.assertTrue(remain)
        for _,B,_,E,V in remain:
            self.assertEqual((B,E,V), (1,0,2))
            self.assertFalse(any(a <= B for a,_ in solutions(180-144, (12,15))))


if __name__ == '__main__':
    unittest.main()
