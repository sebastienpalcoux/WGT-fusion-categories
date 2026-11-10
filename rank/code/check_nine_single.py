"""Exact arithmetic cross-check of the analytic multiplicity-one proof."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm


def check_thresholds():
    gap = F(1,1800)
    assert F(1,20) > gap
    assert F(3,2)+F(5,4) < 3-gap
    assert F(1,2)+F(4,5) < F(3,2)-gap
    assert F(1,4)+F(1,7)<F(1,2)-gap
    assert F(1,5)+F(2,8)<F(1,2)-gap
    assert min(F(1,8)-F(1,10),F(1,10)-F(1,11))==F(1,110)>gap
    assert F(3,7)<F(1,2)-gap
    assert F(1,7)+F(1,10)<F(1,4)-gap
    assert min(F(1,26)-F(1,28),F(1,28)-F(1,29))==F(1,812)>gap


def remaining_lists():
    """Optional enumeration of the four positive terms, unused by the proof."""
    out=[]
    def visit(bs,rem,n):
        count=4-len(bs)
        if count==0:
            if n<1800 or abs(rem)!=F(1,n):
                return
            ds=[n//j for j in (1,2,2,2,*bs)]
            if (min(ds)<6 or reduce(gcd,ds)!=1
                    or any(gcd(x,y)==1 for x,y in combinations(ds,2))):
                return
            out.append((n,tuple(bs)))
            return
        if rem<=0:
            # Remaining terms are positive and each is a nonunit degree.
            return
        # rem <= count/b + 1/n <= (count+1)/b.
        for b in range(bs[-1] if bs else 4,int((count+1)/rem)+1):
            if b%3==1 and not (b==4 and 4 in bs):
                visit(bs+[b],rem-F(1,b),lcm(n,b))
    visit([],F(1,2),2)
    return out


if __name__=='__main__':
    check_thresholds()
    rows=remaining_lists()
    assert rows==[(1820,(4,7,10,130)),(2660,(4,7,10,133)),
                  (4760,(4,7,10,136)),(19460,(4,7,10,139)),
                  (9940,(4,7,10,142)),(4060,(4,7,10,145))]
    for n,bs in rows:
        deficit=1+sum(n//b for b in (1,2,2,2,*bs))-3*n
        assert deficit in (0,2)
        print(n,bs,deficit)
