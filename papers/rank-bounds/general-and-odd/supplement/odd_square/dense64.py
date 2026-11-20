from fractions import Fraction as F
from dense_phase import x64
from full_blocks import count as free_count
def try64(s,x,xh,witnesses):
    root=27*25*343;lam=F(64,45)
    J=F(64,225);Y=F(8,225)
    q7=witnesses(7,s);fq=max(q7,free_count(7,s))
    q5=witnesses(5,s)
    m0=F(s,81)+J+Y
    c0=J*F(7,8)+Y*F(5,4)
    # Let O3/O5 be possible overlaps of the omission groups with x.
    for o3 in range(0,9,2):
        for o5 in range(0,q5+1,2):
            if o3+o5>x:continue
            neutral=(F(q7,(9*343)**2)+F(fq-q7,343**2)+F(q5-o5,225**2)
                +F(8-o3,27**2)+F(1,root**2))
            cap=209-(1+s+64+8+x+fq+q5+8-o3-o5)
            if cap<42:continue
            # budget upper for X: overlap b>=135 resp 225. Nonoverlap high>=45.
            eh=max(0,xh-o3-o5)
            ceiling=F(x-o3-o5-eh,225)+F(eh,45**2)+F(o3,135**2)+F(o5,225**2)
            for u in range(cap//12+1):
                for v in range(2):
                    n7=12*u+42*v
                    if not 42<=n7<=cap:continue
                    for uh in range(u+1):
                        for vh in range(v+1):
                            spare=cap-n7
                            m7=F(12*uh+42*vh,49**2)
                            g7=F(14*uh+42*vh,49**2)
                            mh=F(n7-12*uh-42*vh,343**2)
                            gh=F(14*(u-uh)+42*(v-vh),343**2)+F(3*(8-o3),2*27**2)
                            for c in range(s//2+1):
                                g=c0+g7+F(c,27)
                                for low,high in x64(x,xh):
                                    upper=min(high,ceiling,1-m0-m7,(1+F(1,root)-g)*F(2,3))
                                    lo=max(low,1-m0-m7-mh-neutral-F(spare,225),
                                        (1-F(1,root)-g-gh-F(spare,486))*F(2,3),
                                        (1+lam-lam/root-m0-m7-lam*g-mh-lam*gh-neutral-F(spare,225))/(1+lam*F(3,2)))
                                    if lo<=upper:return (o3,o5,n7,uh,vh,c,spare,lo,upper)
    return False
