"""Run with python test_three.py from any directory."""
import unittest
from fractions import Fraction
from itertools import permutations, product
from math import comb
from three import SPECIAL, covered_three, g, is_prime, quotient_possible, support, word


class FamilyTests(unittest.TestCase):
    def test_scalar(self):
        for s in range(3, 500, 2):
            if not is_prime(s):
                continue
            for k in range(2, 11):
                e = comb(k,2)
                c = 1+max(Fraction(h*(h-1),g(h)) for h in range(1,e+1))
                for odd in (True, False):
                    lhs = Fraction(2*s,3) if odd else Fraction(s,2)
                    self.assertEqual(support(s,k,odd),1+2*e<lhs and (s-1)//2>c)

    def test_bounds_and_gaps(self):
        for a in permutations((7,5,4)):
            self.assertTrue(covered_three((2,3,5),a))
        for p,q,r,a,b,c in [(2,3,7,6,97,13),(3,5,7,6,99,4),
                               (5,7,29,6,99,61),(7,11,53,6,99,61),
                               (11,19,127,99,99,4),(13,17,83,99,99,4),
                               (5,7,1013,57,6,62)]:
            self.assertTrue(covered_three((p,q,r),(a,b,c)))
        for triple, exponents in [((2,3,5),(6,6,4)),((5,7,23),(6,6,4)),
                                   ((3,5,7),(6,8,8))]:
            self.assertIsNone(covered_three(triple,exponents))

    def test_quotient(self):
        for b,cap in [(3,2),(4,2),(3,3),(5,1)]:
            H = b*cap
            for h in range(1,H+1):
                j = g(h)
                feasible = any(sum(x)==h and
                               all(sum(x[i]*x[(i+t)%b] for i in range(b)) >=
                                   cap*j+(h-j if t==0 else 0) for t in range(b))
                               for x in product(range(cap+1),repeat=b))
                self.assertEqual(quotient_possible(H,b,h,j,[10000]),feasible)
                self.assertTrue(quotient_possible(H,b,h,j,[0]))

    def test_odd_word(self):
        self.assertFalse(word(23,4))
        self.assertTrue(word(29,4))
        self.assertFalse(word(7,3))
        for s,k in SPECIAL.items():
            self.assertTrue(word(s,k))


if __name__ == '__main__':
    unittest.main()
