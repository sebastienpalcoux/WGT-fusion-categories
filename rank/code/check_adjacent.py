"""Exact rational-twist checks for the near-saturated denominator 4,5 classes."""
def possibilities(i,j):
    t=400-25*i-16*j
    if not 1<=t<=17:
        raise ValueError('outside the variance range')
    return [(a,b) for a in range(-i,i+1,2) for b in range(-j//4-1,j+1)
            if 4*b>=-j and b%5==j%5 and abs(25*a+16*b)<t]

if __name__=='__main__':
    excluded=[(i,j) for i in range(1,16) for j in range(1,25)
              if 1<=400-25*i-16*j<=17 and not possibilities(i,j)]
    assert excluded==[(1,23),(2,21),(3,20),(5,17),(7,14),(10,9),(12,6),(14,3),(15,1)]
    assert possibilities(11,7)==[(-1,2)]
    assert possibilities(13,4)==[(1,-1)]
    print(excluded)
