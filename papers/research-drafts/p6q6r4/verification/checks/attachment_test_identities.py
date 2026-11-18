"""Exact small fusion-ring identity tests; standard library only."""
from fractions import Fraction as Q
from itertools import permutations,product
import json


def test(name,elements,unit,mul,inv,linear,s,F,symplectic=False):
 n=len(elements);assert len(linear(unit))+4==n
 f=n # All linear formal codegrees in a finite group ring.
 def A(i,j):return sum(Q(a*b,f) for a,b in zip(linear(i),linear(j)))
 def t(i):return Q(s(i),2)
 def det(i):return Q(s(i)**2-F*(int(i==inv(i))-A(i,i)),2)
 def G(i,j):return int(i==inv(j))-A(i,j)-2*t(i)*t(j)/F
 square=0;jordan=0;norm=0
 for i,j in product(elements,repeat=2):
  z=sum(Q((a*a-s(i)*a+det(i))*b,f) for a,b in zip(linear(i),linear(inv(j))))+s(i)*(i==j)-det(i)*(j==unit)
  assert z==int(mul(i,i)==j),(name,'square',i,j,z)
  square+=1
  if symplectic:
   w=sum(Q((a*b-det(i))*c,f) for a,b,c in zip(linear(i),linear(inv(i)),linear(inv(j))))+det(i)*(j==unit)
   assert w==int(mul(i,inv(i))==j)
   norm+=1
 for i,j,k in product(elements,repeat=3):
  z=inv(k)
  L=sum(Q(a*b*c,f) for a,b,c in zip(linear(i),linear(j),linear(z)))
  result=L+2*t(i)*t(j)*t(k)/F+t(i)*G(j,z)+t(j)*G(i,z)+t(k)*G(i,j)
  assert result==Q(int(mul(i,j)==k)+int(mul(j,i)==k),2),(name,'Jordan',i,j,k,result)
  jordan+=1
 return dict(group=name,square_coefficients=square,jordan_coefficients=jordan,norm_coefficients=norm,status='PASS')

S3=tuple(permutations(range(3)));identity=(0,1,2)
def pmul(p,q):return tuple(p[q[i]] for i in range(3))
def pinv(p):return tuple(p.index(i) for i in range(3))
def sign(p):return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
results=[test('S3',S3,identity,pmul,pinv,lambda p:(1,sign(p)),lambda p:sum(p[i]==i for i in range(3))-1,3)]
# Q8 element (sign, axis), axes 0=1, 1=i, 2=j, 3=k.
Q8=tuple(product((1,-1),range(4)))
def qmul(x,y):
 a,i=x;b,j=y
 if i==0:return(a*b,j)
 if j==0:return(a*b,i)
 if i==j:return(-a*b,0)
 return(a*b*(1 if (i,j) in ((1,2),(2,3),(3,1)) else -1),6-i-j)
def qinv(x):return x if x[1]==0 else(-x[0],x[1])
def qlinear(x):
 _,axis=x
 return tuple(1 if axis==0 else a if axis==1 else b if axis==2 else a*b for a,b in product((1,-1),repeat=2))
results.append(test('Q8',Q8,(1,0),qmul,qinv,qlinear,lambda x:2*x[0] if x[1]==0 else 0,4,True))
print(json.dumps(results,indent=2))
