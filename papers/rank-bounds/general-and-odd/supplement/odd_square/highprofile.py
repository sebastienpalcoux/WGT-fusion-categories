from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import isqrt,prod
from dense_phase import groups
from full_blocks import count as free_count, neutral_pair_floors

def merge(intervals):
    acc=[]
    for lo,hi in sorted(set(intervals)):
        if not acc or lo>acc[-1][1]:acc.append([lo,hi])
        else:acc[-1][1]=max(acc[-1][1],hi)
    return tuple(map(tuple,acc))

@lru_cache(None)
def phase(p,k,c,h,v):
    if not h:return ()
    d=(p-1)&-(p-1); eta=k+2*c
    if k==(p-1)*c:
        allowed=[e for e in range(0,h+1,4)
          if 2*d*k >= (p-1)*(2*d+e*(d-1))]
        if not allowed:return ()
        if allowed==[0]:
            if not c or (p+1)*v>c:return ()
            V=F(c,p+1)
            if V.denominator%2==0 or h%8!=(V.numerator*pow(V.denominator,-1,8))%8:return ()
            if v>V or v+F(h-v,9)<V or (h>v and V==v):return ()
            return ((V,V),)
    levels=(1,3)
    out=[]
    def visit(ix,left,tot,first,spent):
        if not left:
            if spent==eta and first==-c:out.append((tot,tot))
            return
        l=levels[ix] if ix<len(levels) else 5
        if spent>eta or (c+first)**2>(eta-spent)*F(left,l*l):return
        if ix==len(levels):
            out.append((tot,tot+F(left,l*l)));return
        for n in ((v,) if ix==0 else range(left+1)):
            if n>left:continue
            for m,e in groups(n,l,p,eta,d):
                visit(ix+1,left-n,tot+F(n,l*l),first+F(m,l*l),spent+F(e,l*l))
    visit(0,h,F(),F(),F())
    return merge(out)

@lru_cache(None)
def sums(p,types,total,high):
    # types = (Galois phase count, joint count/phase, primitive count/phase, order)
    out={}
    def rec(i,left,large,low,hi,low9,hi9):
        if i==len(types):
            if left==0 and large>=high:
                out.setdefault((low9,hi9),[]).append((low,hi))
            return
        m,k,c,order=types[i]
        if not k:
            rec(i+1,left,large,low,hi,low9,hi9);return
        for h in range(1,left//(2*m)+1):
            if (len(types)==i+1 and 2*m*h!=left):continue
            for v in range(h+1):
                for a,b in phase(p,k,c,h,v):
                    a,b=F(2*m,225)*a,F(2*m,225)*b
                    # Track the ordinary contribution of primitive-third
                    # sectors separately, all order-nine charge rates are one.
                    rec(i+1,left-2*m*h,large+2*m*(h-v),
                        low+a,hi+b,low9+(a if order==9 else 0),hi9+(b if order==9 else 0))
    rec(0,total,0,F(),F(),F(),F())
    return [(lo,hi,l9,h9) for (l9,h9),vals in out.items() for lo,hi in merge(vals)]

def possibilities(P,marked,s,shape,witnesses):
    z,x,y,zh,xh,yh=shape
    R=209
    if zh or z not in (96,128) or s<9 or set(P) not in ({3,5,7},{3,5,7,11}):return True
    fixed=1+s+z+x+y+sum(m for _,m,a in marked.values())
    if R-fixed>=42:return True
    hs={3:9,5:5,7:343,**{p:p for p in marked}}
    root=prod(hs.values());three=len(P)==3
    if three:
        root=27*25*343
        q5=witnesses(5,s);q7=witnesses(7,s);fq=max(q7,free_count(7,s))
        oms=[(o3,o5,8-o3+q5-o5+fq,F(8-o3,27**2)+F(q5-o5,225**2)+F(q7,(9*343)**2)+F(fq-q7,343**2),F(3*(8-o3),2*27**2))
             for o3 in range(0,9,2) for o5 in range(0,q5+1,2) if o3+o5<=x]
    else:
        ns=neutral_pair_floors(P,marked,hs,s)
        oms=[(0,0,2*len(ns),sum(F(2,b*b) for b in ns),F())]
    M=F(s,81)+F(z,225)+sum(F(m,p*p) for p,(_,m,a) in marked.items())
    GP=sum(F(a,p*(p-1)) for p,(_,m,a) in marked.items())
    lam=F(64,45)
    result=[]
    for a,b,c in product(range(0,10,2),repeat=3):
        if b>c or 16*a+24*(b+c)!=z:continue
        GJ=GP+F(14*a+24*(b+c),225)
        X=sums(5,tuple((g,4*t,t,order) for g,t,order in [(1,a,3),(1,a,3),(3,b,9),(3,c,9)]),x,xh)
        Y=sums(3,((2,2*a+3*(b+c),a,5),)*2,y,yh)
        for o3,o5,n,neutral,c_wit in oms:
            # Add an inactive pure 7-pair disjoint from every omission witness.
            available=R-fixed-n-2
            if available<12:continue
            xcap=F(x-o3-o5-max(0,xh-o3-o5),225)+F(max(0,xh-o3-o5),45**2)+F(o3,135**2)+F(o5,225**2)
            for n7 in range(12,available+1,12):
                for t7 in range(0,n7+1,12):
                    m7=F(t7,49**2);g7=F(7*t7,6*49**2)
                    mh=F(n7-t7,343**2);gh=F(7*(n7-t7),6*343**2)
                    spare=available-n7
                    for c9 in range(s//2+1):
                        g=GJ+g7+F(c9,27)
                        for xl,xu,x9l,x9u in X:
                            xu=min(xu,xcap)
                            if xl>xu:continue
                            for yl,yu,_,_ in Y:
                                maxmass=M+m7+mh+neutral+F(2,49**2)+F(1,root**2)+F(spare,225)+xu+yu
                                if maxmass<1 or M+m7+xl+yl>1:continue
                                gl=g+F(3,2)*xl-x9u/2+F(5,4)*yl
                                gu=g+F(3,2)*xu-x9l/2+F(5,4)*yu
                                if gl>1+F(1,root) or gu+gh+c_wit+F(spare,486)<1-F(1,root):continue
                                if maxmass+lam*(gu+gh+c_wit)<1+lam-lam/root:continue
                                result.append((a,b,c,n7,t7,spare,c9,xl,xu,yl,yu))
                                break
                            else:continue
                            break
    return result
