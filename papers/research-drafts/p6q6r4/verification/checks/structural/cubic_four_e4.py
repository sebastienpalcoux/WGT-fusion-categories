#!/usr/bin/env python3
"""Complete tiny dimension inventories after the e4 local product identities."""
from math import isqrt
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    A=[]
    for a in range(3,42):
        for q in range(3,36):
            p=q+3
            if 3*a*a+2*p*p+2*q*q+64!=5117:continue
            assert 2*p>23
            A.append({'a':a,'p':p,'q':q,'scalar_I_minimum':2*p,'available_mass':23})
    B=[]
    for p in range(3,42):
        for z in range(72):
            if (z+p+3)%2:continue
            x=((z+3+p)//2,(z+3-p)//2,(3-z+p)//2,(3-z-p)//2)
            if min(map(abs,x))<3:continue
            if z>=p+9:minimum=z+3
            else:
                assert p>=z+9
                minimum=3*p+3
            rem=5044-z*z-3*p*p
            if rem<0:continue
            extra=isqrt(rem)
            if extra*extra!=rem or (extra and extra<3):continue
            if minimum>23:cert={'reason':'scalar_I_mass','minimum':minimum,'available':23}
            else:
                assert (p,z,extra)==(3,16,69)
                small=(1,4,4,4,4,5,7,12,8)+tuple(map(abs,x))+(p,p)
                assert len(small)==15 and max(small)==12 and 4*max(small)<extra
                cert={'reason':'proper_generator_core','generator_dimension':4,'core_maximum_dimension':12,'outside_dimension':69}
            B.append({'p':p,'z':z,'signed_dimensions':x,'extra':extra,'certificate':cert})
    assert len(A)==2 and len(B)==21
    out={'scope':'positive cubic d4 e4 T7, unique M2 f6, ranks <=16','all_excluded':True,'A':A,'B':B}
    (HERE/'cubic_four_e4_certificates.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'all_excluded':True,'A_types':len(A),'B_types':len(B),'B_scalar_mass_exclusions':sum(r['certificate']['reason']=='scalar_I_mass' for r in B),'B_generator_exclusions':sum(r['certificate']['reason']=='proper_generator_core' for r in B)},indent=2))

if __name__=='__main__':main()
