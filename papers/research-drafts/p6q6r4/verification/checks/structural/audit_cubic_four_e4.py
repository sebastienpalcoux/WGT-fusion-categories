#!/usr/bin/env python3
"""Independent direct integer inventory and cone check for the e4 branch."""
from itertools import product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    data=json.loads((HERE/'cubic_four_e4_certificates.json').read_text())
    A={(a,q+3,q) for a,q in product(range(3,42),range(3,36)) if 3*a*a+2*(q+3)**2+2*q*q==5053}
    assert A=={(r['a'],r['p'],r['q']) for r in data['A']}
    for a,p,q in A:assert 2*p>23
    # Solve the square equation by a table rather than the producer's isqrt.
    sq={e*e:e for e in (0,)+tuple(range(3,72))};B=set()
    for p,z in product(range(3,42),range(72)):
        if (p+z)%2!=1:continue
        x=tuple(v//2 for v in (z+3+p,z+3-p,3-z+p,3-z-p))
        if any(abs(a)<3 for a in x):continue
        if 5044-z*z-3*p*p not in sq:continue
        B.add((p,z,sq[5044-z*z-3*p*p],x))
    assert B=={(r['p'],r['z'],r['extra'],tuple(r['signed_dimensions'])) for r in data['B']}
    direct_cones=0;cores=0
    for p,z,e,x in B:
        # Directly enumerate every selfdual active coefficient vector of
        # mass <=23 and test whether its matrix image is I.
        s=tuple(1 if a>0 else -1 for a in x);d=tuple(map(abs,x));hit=False
        for n in product(*(range(23//a+1) for a in d)):
            f=tuple(a*b for a,b in zip(s,n));k=f[0]-f[1]
            if k<0 or f[2]-f[3]!=k or sum(f)!=2:continue
            if sum(a*b for a,b in zip(d,n))+2*p*k<=23:hit=True;break
        if not hit:direct_cones+=1
        else:
            assert (p,z,e,x)==(3,16,69,(11,8,-5,-8))
            assert 4*max(12,*d,p)<e
            cores+=1
    result={'passed':True,'A_types':len(A),'B_types':len(B),'B_direct_scalar_cone_exclusions':direct_cones,'B_generator_core_exclusions':cores,
            'independent_methods':['direct dimension equations','square lookup','four-variable scalar-I cone enumeration']}
    assert direct_cones==20 and cores==1
    (HERE/'cubic_four_e4_independent_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
