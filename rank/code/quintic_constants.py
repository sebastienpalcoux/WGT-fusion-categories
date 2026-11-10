"""Exact checks of the margins and finite phase cases in the 14/2 and 13/4 prefix proofs."""
from fractions import Fraction as F


def verify(M0=5400):
    for j in (2,3):
        delta = 1-F(14,16)-F(j,25)
        width=delta+F(1,M0)
        assert 256*width**2<1
        assert 16*(F(j,25)+width)<3
        assert F(1,8)-F(j,100)>width
        assert 25*(F(1,8)-width)>j-1
        if j==2:assert (25*width)**2<F(3,2)
        else:assert 25*width<1
        assert delta*M0*M0>1
    delta=F(11,400)
    width=delta+F(1,M0)
    assert 256*width**2<1 and 625*width**2<F(5,4)
    possible=[(A,B) for A in range(-13,14,2) for B in (4,-1)
              if abs(F(A,16)+F(B,25))<=width]
    assert possible==[(-3,4),(1,-1)]
    defect=delta-(F(1,16)-F(1,25))+F(1,M0)
    pure_five_mass=F(4,3)*defect
    assert 25*pure_five_mass<1 and 625*pure_five_mass<5
    assert 625*defect/2<2
    return possible


if __name__=='__main__':
    print('proved admissible quintic phase pairs before the local step',verify())
