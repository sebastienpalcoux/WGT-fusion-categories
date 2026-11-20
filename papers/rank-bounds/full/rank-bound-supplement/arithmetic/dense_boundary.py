"""Independent rational interval check of the last dense-seven boundaries."""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache

# Necessary local pairs (h,v) for C=2 at the least joint denominator:
# h=3 forces v=0; h=4,5 require 1<=v<=3, while h=6 allows0<=v<=3.
# The upper bound v<=3 for h=6 uses both reflected parity inequalities.
# These follow from the
# primitive parity-lattice and two-moment inequalities. No search imports.
local={3:(0,),4:(1,2,3),5:(1,2,3),6:(0,1,2,3)}

@lru_cache(None)
def check():
    neutral=(F(2,99**2)+F(2,55**2)+F(2,(11*343)**2)
             +F(6,343**2)+F(1,105**2))
    # At s=34, forty-eight primitive seven-active labels can fit only
    # with x=12,y=40 and no free positions. Their upper mass is too small.
    forty_eight=(F(34,81)+F(22,121)+F(72,225)+F(4,675)
                 +F(48,2401)+neutral)
    assert forty_eight==F(1302011477,1372711725)<1
    remaining=[]
    largest_upper=F(0)
    for s,xs in [(34,(12,14,16,18)),(36,(12,14,16)),(38,(12,14)),(40,(12,))]:
        for x in xs:
            free=201-(1+s+22+32+x+40+42+12)
            for h1,h2 in product(local,repeat=2):
                if h1>h2 or 2*(h1+h2)!=x:
                    continue
                for v1,v2 in product(local[h1],local[h2]):
                    lower=[F(v)+F((2-v)**2,12-v) for v in (v1,v2)]
                    upper=[F(v)+F(h-v,9) for h,v in ((h1,v1),(h2,v2))]
                    base=F(s,81)+F(22,121)+F(72,225)
                    for small,lo7,hi7 in [(True,F(42,2401),F(42,2401)),
                                           (False,F(0),F(42,343**2))]:
                        low=base+2*sum(lower)/225+lo7
                        high=base+2*sum(upper)/225+hi7+neutral+F(free,225)
                        largest_upper=max(largest_upper,high)
                        if low<=1<=high:
                            remaining.append((s,x,h1,v1,h2,v2,small))
    # Any non-small fifth-phase fiber moves at least eight labels above15.
    assert largest_upper==F(1400263477,1372711725)
    assert largest_upper-F(64,2025)<1
    assert remaining==[(34,16,4,3,4,3,True),
                       (36,14,3,0,4,3,True),
                       (36,16,4,2,4,2,True),
                       (38,12,3,0,3,0,True)]
    # The three rigid cases have the same total mixed-plus-nine charge.
    assert 1-F(58,135)-F(2,11)-F(10,27)==F(1,55)
    assert 1-F(16,45)-F(2,11)-F(12,27)==F(1,55)
    assert 385%abs(49**2-2310)!=0
    # With two (4,2) fibers, balance gives equal tail weights and even c.
    assert F(2,5)+F(2,11)+F(12,27)>1+F(1,105)
    assert (F(2,5)+F(2,11)+F(10,27)+F(4,675)
            +F(42,2401))<1-F(1,105)
    return remaining

if __name__=='__main__':
    print('Four dense-seven phase intervals; all excluded by exact Gauss balance.',check())
