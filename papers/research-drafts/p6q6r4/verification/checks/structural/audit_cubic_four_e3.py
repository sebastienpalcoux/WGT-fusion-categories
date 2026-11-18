#!/usr/bin/env python3
"""Independent direct square and matrix-cone audit of the e3 final cases."""
from itertools import combinations_with_replacement,product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    data=json.loads((HERE/'cubic_four_e3_certificates.json').read_text())
    sq={}
    for d,e in combinations_with_replacement((0,)+tuple(range(3,73)),2):sq.setdefault(d*d+e*e,[]).append((d,e))
    regenerated=set()
    for p,z in product(range(3,8),range(16)):
        if z<p+8 or(p-z)%2:continue
        for d,e in sq.get(5194-z*z-3*p*p,[]):regenerated.add((p,z,d,e))
    assert regenerated=={(r['p'],r['z'],*r['extras']) for r in data['descriptions']}
    for p,z,d,e in regenerated:
        signed=(9,5,-3,-7);dims=tuple(map(abs,signed));signs=(1,1,-1,-1)
        zero_dimensions=(8,d,e);reachable={0}
        for mass in range(1,18):
            if any(mass-v in reachable for v in zero_dimensions):reachable.add(mass)
        for n in product(*(range(17//v+1) for v in dims)):
            f=tuple(a*b for a,b in zip(n,signs));k=f[0]-f[1]
            if k<0 or f[2]-f[3]!=k or sum(f)!=2:continue
            mass=sum(a*b for a,b in zip(n,dims))+2*p*k
            assert 17-mass not in reachable
    result={'passed':True,'descriptions':len(regenerated),'independent_methods':['raw two-square table','four-variable matrix equations','zero-image dimension semigroup']}
    (HERE/'cubic_four_e3_independent_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
