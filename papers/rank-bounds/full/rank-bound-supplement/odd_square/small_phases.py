"""Exact denominator-fifteen occupancies in the small joint phase families."""
from functools import lru_cache
from edges import axis

@lru_cache(None)
def small_counts(z,big_z,cap):
    """Maps the number at fifteen on either side to all possible full multiplicities."""
    if z not in (32,48,64) or big_z not in (0,z//2):
        return None
    if z==64 and big_z:
        return None
    # (p complementary to the fiber, multiplicity of its orbit, full joint
    # count, joint count in one primitive p-phase), all counts per dual pair.
    if z in (32,64):
        c=z//16
        three=[(5,1,4*c,c)]*2
        five=[(3,2,2*c,c)]*2
    else:
        three=[(5,3,8,2)]
        five=[(3,2,6,0)]*2
    def combine(axes):
        costs={0:{0}}
        for p,w,k,c in axes:
            j=k*(z-big_z)//z
            uc=c*(z-big_z)//z
            opts=[]
            for h in range(1,cap//(2*w)+1):
                for v in range(h+1):
                    if axis(p,k,j,h,v,c,uc):
                        opts.append((v,h))
            new={}
            for small,old_counts in costs.items():
                for v,h in opts:
                    key=small+2*w*v
                    for old in old_counts:
                        total=old+2*w*h
                        if total<=cap:
                            new.setdefault(key,set()).add(total)
            costs=new
        return {small:frozenset(counts) for small,counts in costs.items()}
    return combine(three),combine(five)
