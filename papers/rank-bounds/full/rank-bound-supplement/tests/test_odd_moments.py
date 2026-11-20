"""Regression tests for the odd-fiber moment and height refinements."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'arithmetic'))
from edges import axis
from dense_boundary import check
from small_phases import small_counts
from joint_gauss import next_level_cap,zero_moments


class OddMomentBounds(unittest.TestCase):
    def test_exact_primitive_moment_models_are_preserved(self):
        # One coefficient -p and (p+1) joint pairs at each primitive phase.
        # Both moment identities are exact, as is the p^2 height congruence.
        for p in (3,5,7):
            c=p+1
            self.assertEqual(-p,1-c)
            self.assertEqual(p*p,1+(p-1)*c)
            self.assertTrue(axis(p,(p-1)*c,(p-1)*c,1,1,c,c))
        # Three coefficient -5 pairs at relative denominator3, C=2,J=8.
        self.assertEqual(3*F(-5,9),F(1,3)-2)
        self.assertEqual(3*F(25,9),F(1,3)+8)
        self.assertTrue(axis(5,8,8,3,0,2,2))
        # A five-pair ternary fiber: coefficients3,0,0,0,0.
        self.assertTrue(axis(3,4,4,5,5,2,2))

    def test_exact_fiber_models_survive_small_occupancy_screen(self):
        # Two ternary phases each have three -5 pairs at relative denominator3;
        # the fifth-root phases each have five pairs at the smallest denominator.
        three,five=small_counts(32,0,84)
        self.assertIn(12,three[0])
        self.assertIn(40,five[40])
        # In the doubled joint family the -5 and -3 models give 24 and8 labels.
        three,five=small_counts(64,0,64)
        self.assertIn(24,three[0])
        self.assertIn(8,five[8])

    def test_half_split_five_fiber_needs_another_pair(self):
        self.assertFalse(axis(5,8,4,4,0,2,1))

    def test_ternary_height_restricts_the_next_denominator(self):
        # Two larger pairs cannot increase3-height, so their relative
        # denominators are at least5, contradicting the required mass1/5.
        self.assertLess(2*F(1,25),F(1,5))
        self.assertFalse(axis(3,4,4,5,3,2,2))

    def test_coprime_cap_preserves_exact_gauss_tails(self):
        # These are actual positive packet sums with either vacuum-error sign.
        # Future denominators need not divide the earlier scale.
        for b,k,charge in ((343,12,F(14,343**2)),(10125,4,F(7,2*10125**2))):
            root=27*25*7**3*11
            from math import lcm
            M=lcm(root,b*b//7 if b%7==0 else b)
            factor=min(7**3*11,7*27*25*11)
            for sign in (-1,1):
                residual=charge+F(sign,M)
                cap=next_level_cap(residual.numerator,residual.denominator,k,factor)
                self.assertGreaterEqual(cap,b)
                # The next integer violates the exact necessary quadratic.
                self.assertGreater(2*factor*abs(residual)*(cap+1)**2,
                                   2*(cap+1)+3*factor*k)

    def test_reflection_retains_exact_models_and_removes_fake_small_counts(self):
        self.assertFalse(axis(5,8,4,4,2,2,1))
        self.assertFalse(axis(3,4,2,6,2,2,1))
        self.assertTrue(zero_moments(5,3,0,72,72))
        self.assertTrue(zero_moments(3,5,4,4,4))
        self.assertTrue(zero_moments(3,5,4,3,5))
        self.assertFalse(zero_moments(3,4,4,4,4))

    def test_dense_boundary_certificate(self):
        self.assertEqual(len(check()),4)


if __name__=='__main__':
    unittest.main()
