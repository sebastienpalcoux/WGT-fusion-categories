from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

@lru_cache(None)
def groups(n,k,p,energy,dyadic):
    """Integer first/squared moments in k^{-2} units at one level."""
    max_a=isqrt(energy*k*k)
    choices=range(-max_a,max_a+1)
    states={0:{(0,0)}}
    for a in choices:
        # a = (coefficient - 1), coefficient divisible by p.
        if (a+1)%p: continue
        tile=dyadic if (a+1)%2==0 else 1
        for j in range(tile,n+1):
            states.setdefault(j,set()).update((m+tile*a,e+tile*a*a)
                for m,e in states.get(j-tile,()) if e+tile*a*a<=energy*k*k)
    return tuple(states.get(n,()))

@lru_cache(None)
def intervals(h,v,order):
    # Four joint labels per primitive full phase; each 5-phase has
    # J=4,C=2 in pairs at order15; J=6,C=0 at order45.
    p=5 if order==64 else 3
    C=4 if order==64 else int(order==15)*2
    energy=24 if order==64 else (8 if C else 6)
    C=-C  # summands tracked below are a_i-1, so the first moment is -C
    levels=(1,3,5)
    out=[]
    def visit(idx,left,total,first,used):
        if left==0:
            if used==energy and first==C:out.append((total,total))
            return
        k=levels[idx] if idx<len(levels) else 9
        if used>energy or (C-first)**2>(energy-used)*F(left,k*k):return
        if idx==len(levels):
            out.append((total,total+F(left,k*k)))
            return
        for n in ((v,) if idx==0 else range(left+1)):
            if n>left:continue
            for m,e in groups(n,k,p,energy,4 if p==5 else 2):
                visit(idx+1,left-n,total+F(n,k*k),first+F(m,k*k),used+F(e,k*k))
    visit(0,h,F(),F(),F())
    if not out:return ()
    # Merge overlapping rational intervals.
    acc=[]
    for lo,hi in sorted(set(out)):
        if not acc or lo>acc[-1][1]:acc.append([lo,hi])
        else:acc[-1][1]=max(acc[-1][1],hi)
    return tuple(map(tuple,acc))

@lru_cache(None)
def yintervals(y,high,z):
    out=[]
    # two square classes, each has two phases, each counted in pairs
    for h in range(y//4+1):
        l=y//4-h
        if y%4 or (z==32 and min(h,l)<5):continue
        for v in range(h+1):
            for w in range(l+1):
                if 4*((h-v)+(l-w))<high:continue
                for a,b in intervals(h,v,15 if z==32 else 45):
                    for c,d in intervals(l,w,15 if z==32 else 45):
                        out.append((4*(a+c)/225,4*(b+d)/225))
    acc=[]
    for lo,hi in sorted(set(out)):
        if not acc or lo>acc[-1][1]:acc.append([lo,hi])
        else:acc[-1][1]=max(hi,acc[-1][1])
    return tuple(map(tuple,acc))

@lru_cache(None)
def x64(x,hi):
    out=[]
    for h in range(x//2+1):
        l=x//2-h
        for v in range(h+1):
            for w in range(l+1):
                if 2*(h+l-v-w)<hi: continue
                for a,b in intervals(h,v,64):
                    for c,d in intervals(l,w,64):
                        out.append((2*(a+c)/225,2*(b+d)/225))
    acc=[]
    for lo,up in sorted(set(out)):
        if not acc or lo>acc[-1][1]:acc.append([lo,up])
        else:acc[-1][1]=max(up,acc[-1][1])
    return tuple(map(tuple,acc))

@lru_cache(None)
def y_five_pair_refinement():
    # Minimal primitive joint family, exactly five solely-five pairs
    # above each of the four primitive fifth phases. At b=15 a nonzero
    # coefficient forces saturation. Otherwise every small coefficient
    # is zero and its positive multiplicity at each phase is even.
    # Four even coefficients in the fiber follow from E=4: E=0 would
    # require the non-2-integral equation V=C/4=1/2.
    # Pair these four equal-degree/equal-coefficient entries in twos.
    ks=[1,3,5,9,15,25,27,45,75]
    one=[(F(5),F(5))]
    for a in ks:
        for b in ks:
            for c in ks[1:]: # one odd coefficient, outside 15 in this branch
                if a>b:continue
                lo=sum((F(m,k*k) for m,k in ((2,a),(2,b),(1,c)) if k!=75),F())
                hi=sum(F(m,k*k) for m,k in ((2,a),(2,b),(1,c)))
                if hi>=F(1,2):one.append((lo,hi))
    both=[]
    for lo,hi in set(one):
        for a,b in set(one):
            both.append((F(4,225)*(lo+a),F(4,225)*(hi+b)))
    return tuple(set(both))
