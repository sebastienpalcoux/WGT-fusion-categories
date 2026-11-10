"""Check the rational endpoints in the two-square rank proposition."""
from fractions import Fraction as F


def lower(x,y,u,v):
    n = x*y*u*v
    return (17-F(16,n*n) + (u-1)*(1-F(16,u*u)) + (v-1)*(1-F(16,v*v))
            + 4-16*(F(1,x*x)+F(1,y*y))*(F(1,u*u)+F(1,v*v))
            + 1-F(16,(u*v)**2))


def endpoint(x,y):
    return lower(x,y,5,7)+2-max(F(16,x*x)+F(16,y*y),F(16,(x*y)**2)+F(16,36))


if __name__ == '__main__':
    values = (lower(8,9,7,11),endpoint(8,27),endpoint(16,9))
    assert values == (F(352721,10164),F(5180713,178605),F(11512811,396900))
    assert min(values)>29
    print(*values,sep='\n')
