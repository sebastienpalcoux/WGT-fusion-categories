"""Inspect candidate types of the positive fermionic quotient (dimension Delta)."""
import json
import sys
from collections import Counter
from math import gcd,lcm

def factors(n):
    answer={}
    p=2
    while p*p<=n:
        while n%p==0:
            answer[p]=1+answer.get(p,0)
            n//=p
        p+=1 if p==2 else 2
    if n!=1:answer[n]=1
    return answer

def covered(P):
    E=sorted(P.values())
    if len(E)<=2:return True
    if len(E)==3:
        a,b,c=P.values()
        return a<=5 or b<=4 or c<=3 or all(x<=k for x,k in zip(E,(4,5,7)))
    if len(E)==4:
        d,c,b,a=E
        return b<=3 or d==2 and c<=3 and (a//2+b//2<=5 or P.get(2,0)>=4 and a<=9 and b<=5)
    if len(E)==5:
        return E[-1]<=3 or (E[-1]==4 and E[-2]<=3 and E[0]==2)
    return False

def four_supports(d,x):
    """Look for m mutually orthogonal integral rows of codegree four."""
    from itertools import combinations
    m=d.count(x)
    candidates=[]
    for i,j,k in combinations(range(1,len(d)),3):
        a,b,c=d[i],d[j],d[k]
        if len({a,b,c})<3 or x%a or x%b or x%c:continue
        if gcd(a,b,c)!=1 or min(gcd(a,b),gcd(a,c),gcd(b,c))==1:continue
        if x%2==0 and any((x//e)%2 for e in (a,b,c)):continue
        signs=a+b-c
        if abs(signs)!=1:continue
        candidates.append({i:-signs,j:-signs,k:signs})
    for selected in combinations(candidates,m):
        if all(sum(z.get(i,0)*v for i,v in w.items())==-1
               for z,w in combinations(selected,2)):
            return True
    return False

def check(d):
    D=sum(z*z for z in d)
    if 2*D<29160000 or covered(factors(2*D)):return 'BDP'
    if not all(gcd(*[y for y in set(d) if y>1 and gcd(x,y)>1 and x*x+lcm(x,y)**2<=D],0)==1 for x in set(d) if x>1): return 'row gcd'
    if any(not four_supports(d,x) for x in set(d) if D==4*x*x):return 'four'
    if any(D==3*x*x for x in d):return 'three'
    return 'survived'

if __name__=='__main__':
    for name in sys.argv[1:]:
        c=Counter(); good=[]
        for row in open(name):
            d=json.loads(row);s=check(d);c[s]+=1
            if s=='survived':good.append(d)
        print(name,c,'remaining first 5',good[:5])
