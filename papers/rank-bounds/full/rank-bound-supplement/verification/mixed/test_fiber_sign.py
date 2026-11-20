"""Arithmetic regression for the oriented first moment and dyadic multiplicities."""
from fractions import Fraction
from dense_phase import groups, intervals

# At the least denominator the four zero coefficients and one coefficient
# three have first moment sum(1-a)=2 and squared norm sum((1-a)^2)=8.
assert (-2,8) in groups(5,1,3,8,2)
assert (Fraction(5),Fraction(5)) in intervals(5,5,15)

# One four-element 5-inertia packet cannot have different coefficients.
assert groups(1,1,5,1,4)==()
assert (-4,4) in groups(4,1,5,4,4)
