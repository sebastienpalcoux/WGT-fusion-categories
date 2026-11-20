"""Small rational-interval test for the pure-five three-prime rank44 relaxation.
Requires the three-prime reservations and inverse-residue lemma in the manuscript.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import isqrt
from row_quota import row_quota
from full_census import census
from three_primes import bounds

def booked(groups, rank=44, other=6):
 remaining=rank-sum(n for n,b in groups)
 return sum(Q(n,b*b) for n,b in groups)+Q(remaining,other*other) if remaining>=0 else Q(0)

def bmin(a):
 return next(b for b in range(2,52) if pow(2,a,125)*pow(3,b,125)%125 in (1,124))

def initial(m,a,b):
 k3=row_quota(3,m);n3=k3+(2-k3)%3;k5=row_quota(5,m);n5=census(5,m)
 h=44-(1+m+n3+3+5+n5)
 A=2**a
 upper=Q(m,16)+Q(h,36)+Q(k3,16*3**(2*b))+Q(n3-k3,3**(2*b))+Q(2,A*A)+Q(1,(5*A)**2)+Q(5,625)+Q(k5,500**2)+Q(n5-k5,125**2)+Q(1,(A*3**b*125)**2)
 return h,n3,k3,k5,n5,upper

def intervals(m,a,b):
 A=2**a;h,n3,k3,k5,n5,upper=initial(m,a,b)
 if upper<1:return []
 assert upper-1<4*(Q(1,36)-Q(1,A*A)),(m,a,b,'odd')
 assert upper-1<3*(Q(1,36)-Q(1,3**(2*b))),(m,a,b,'3-free')
 evenfloor=9 if a==3 else 8
 extra_max=int((upper-1)/(Q(1,36)-Q(1,evenfloor**2)))
 odd=[(A*3**j) for j in range(min(b,8)) if A*3**j<200]+[0]
 ordinary=sorted(set(2**i*3**j for i in range(min(a,8)) for j in range(min(b,7)) if 6<2**i*3**j<200))+[0]
 results=[]
 for c in range(2,6):
  B=5**c;Mp=B*B//5
  odd5=A*25
  if c==2:
   odd5=next(A*25*3**j for j in range(100) if min(pow(A*3**j,-1,25),-pow(A*3**j,-1,25)%25)<=2)
  tail=Q(k3,16*3**(2*b))+Q(n3-k3,3**(2*b))+Q(1,(A*3**b*Mp)**2)+Q(1,odd5**2)
  # At B25 split exactly whether the sole possible non-hit is denominator125.
  packets=[(0,Q(k5,(4*Mp)**2)+Q(n5-k5,Mp**2))]
  if c==2 and n5-k5==1:
   packets=[(1,Q(k5,27000**2)),(0,Q(k5,27000**2)+Q(1,2250**2))]
  for at125,free_tail in packets:
   fixed=tail+free_tail+(Q(5,B*B) if c==5 else 0)
   target=1-Q(m,16)-(Q(5,B*B) if c<5 else 0)-Q(at125,125**2)
   for extra in range(extra_max+1):
    for d,e in combinations_with_replacement(odd,2):
     need=target-Q(h-extra,36)-(Q(1,d*d) if d else 0)-(Q(1,e*e) if e else 0)
     margin=fixed+Q((not d)+(not e),max(A,200)**2)
     for xs in combinations_with_replacement(ordinary,extra):
      rem=need-sum(Q(1,x*x) for x in xs if x)
      slack=margin+Q(xs.count(0),50**2)
      # Vacuum is strictly positive, so equality at zero is impossible.
      if 0<rem<=slack:results.append((m,a,b,c,at125,d,e,xs,rem,slack))
 return results


def finish_intervals(rows):
 expected={(9,29,3,2,1),(12,36,2,2,1),(12,36,2,2,0),(12,36,2,3,0),(12,36,2,4,0),(12,36,2,5,0)}
 assert {row[:5] for row in rows}==expected and len(rows)==6
 residual=Q(79,162000)
 tiny=Q(3,2**58)+Q(1,(2**29*27*125)**2)
 assert Q(4,108**2)+Q(1,216**2)+Q(1,125**2)+Q(4,27000**2)+tiny<residual
 assert Q(5,108**2)+Q(1,125**2)>residual
 values=[Q(1,(36*2**i)**2) for i in range(7)]+[Q(0)]
 six={}
 for terms in combinations_with_replacement(values,6):
  value=sum(terms);six[value]=max(six.get(value,0),terms.count(0))
 values=[Q(1,(9*2**i)**2) for i in range(9)]+[Q(0)]
 two=[(sum(terms),terms.count(0)) for terms in combinations_with_replacement(values,2)]
 base=Q(12,16)+Q(8,36)
 for c in range(2,6):
  B=5**c;P=B*B//5
  target=1-base-(Q(5,B*B) if c<5 else 0)
  tiny=Q(Q(6,16)+1,P**2)+Q(3,2**72)+Q(1,(2**36*9*125)**2)
  if c==5:tiny+=Q(5,B*B)
  for a,za in six.items():
   for b,zb in two:
    assert not 0<target-a-b<=tiny+Q(za,900**2)+Q(zb,225**2)

def preliminary44():
 for (p,t,mode,m),base in bounds(44):
  if base<1:continue
  if p>3:
   g=[(1,900),(m,4),(t,t*t),(2,8),(1,8*t),(row_quota(p,m),4*p*p),(row_quota(t,m),4*t*t)]
   assert mode=='pure' and booked(g,44,8)<1
  elif t==7:
   assert mode=='pure'
   if m==3:g=[(1,900),(m,4),(7,49),(1,9),(1,441),(1,343),(1,1029),(1,196)]
   else:
    k3=row_quota(3,m);k7=row_quota(7,m);n7=census(7,m)
    g=[(1,900),(m,4),(7,49),(k3,36),((2-k3)%3,9),(2,8),(1,56),(k7,1372),(n7-k7,343)]
   assert booked(g,44)<1,(p,t,mode,m,booked(g,44))
  elif mode=='dyadic':
   assert (p,t,m)==(3,5,3)
   # These odd-witness reservations even suffice at rank45, where m2
   # is the only additional dyadic case.
   assert booked([(1,900),(3,4),(16,10),(2,8),(1,40)],45)<1
   assert booked([(1,900),(2,4),(16,10),(1,20)],45)<1
   assert booked([(1,900),(3,4),(16,10),(3,12),(1,20)],45)<1
  else:assert (p,t,mode)==(3,5,'pure')

def main():
 preliminary44()
 smooth=sorted({2**i*3**j for i in range(8) for j in range(6) if 2**i*3**j<=216})
 permitted=[k for k in smooth if min(pow(k,-1,125),-pow(k,-1,125)%125)<=11]
 assert permitted[0]==1 and next(k for k in permitted if k>1)==18
 assert next(k for k in permitted if k%4==0)==216
 assert next(j for j in range(10) if min(pow(8*3**j,-1,25),-pow(8*3**j,-1,25)%25)<=2)==1
 assert next(j for j in range(10) if min(pow(16*3**j,-1,25),-pow(16*3**j,-1,25)%25)<=2)==5
 for m in range(1,4):
  assert booked([(1,900),(m,4),(20,25)])<1
  old=[(1,900),(m,4),(5,25),(1,45),(1,9),(1,375),(1,20),(2,125)]
  if m<3:assert booked(old)<1
  else:
   assert booked(old)-1<3*(Q(1,36)-Q(1,625))
   assert booked([(1,900),(3,4),(5,25),(1,45),(1,9),(1,375),(2,125),(2,8),(1,200)])<1
   assert booked([(1,900),(3,4),(5,25),(1,45),(1,9),(1,375),(2,125),(3,12),(1,100)])<1
 for m in range(4,14):
  k=row_quota(3,m)
  shared=[(1,900),(m,4),(k,36),((2-k)%3,9),(2,8),(1,40)]
  assert booked(shared+[(20,25)])<1
  k5=row_quota(5,m);n5=census(5,m)
  baseline=booked(shared+[(5,25),(k5,500),(n5-k5,125)])
  if baseline>=1:assert baseline-1<3*(Q(1,36)-Q(1,625))

 out=[]
 for m in range(4,14):
  for a in range(3,58):
   b=bmin(a)
   for bb in (b,b+50):out.extend(intervals(m,a,bb))
 finish_intervals(out)
 print('Three-prime support through rank44: six final intervals, all excluded.')

if __name__=='__main__':main()
