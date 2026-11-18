#!/usr/bin/env python3
"""The complete small B inventory in the e3 branch with W8 distinct from T8."""
from math import isqrt
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    records=[]
    for p in range(3,42):
        for z in range(73):
            if (z+p)%2:continue
            x=((z+2+p)//2,(z+2-p)//2,(2-z+p)//2,(2-z-p)//2)
            if min(map(abs,x))<3:continue
            if z>=p+8:cost=z+2
            else:
                assert p>=z+8
                cost=3*p+2
            if cost>17:continue
            mass=5194-z*z-3*p*p
            for d in (0,)+tuple(range(3,isqrt(mass//2)+1)):
                e=isqrt(mass-d*d)
                if e<d or d*d+e*e!=mass or e<3:continue
                assert (p,z)==(4,12) and x==(9,5,-3,-7)
                assert (d,e) in ((39,59),(49,51))
                records.append({'p':p,'z':z,'extras':(d,e),'signed_dimensions':x,
                                'scalar_base_mass':14,'remaining_mass':3,
                                'least_zero_image_dimension':8,'least_extra_active_cone_mass':24})
    assert len(records)==2
    out={'scope':'positive cubic d4 e3 T8, W8 distinct from T8, B support, ranks <=16',
         'all_excluded':True,'descriptions':records}
    (HERE/'cubic_four_e3_certificates.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'descriptions':len(records),'all_excluded':True},indent=2))

if __name__=='__main__':main()
