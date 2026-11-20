"""Exact scalar checks for the rank-44 odd-edge isolation extension.
No degree-type enumeration is performed.
"""
from fractions import Fraction as F
import json
from pathlib import Path

vac=F(1,210**2)
checks={
 'two_odd_edges_15_21':F(12,225)+F(18,441)+F(14,16),
 'edge21_and_dyadic5':F(16,100)+F(18,441)+F(10,16),
 'four_odd_primes_with_edge':F(11,16)+F(1,15)+F(1,7)+F(1,11),
 '21_absent5':F(13,16)+F(22,441)+F(5,625)+F(3,36)+vac,
 '21_present5_at_m11':F(11,16)+F(22,441)+F(1,5)+F(4,400)+F(1,36)+vac,
 '21_present5_at_m12_rigid':F(12,16)+F(1,5)+F(73,1764)+F(2,400)+F(1,900)+F(1,2025)+F(1,420**2),
 '15_u11_m12':F(12,16)+F(18,225)+F(1,11)+F(2,36)+vac,
 '15_u11_absent':F(13,16)+F(18,225)+F(11,121**2)+F(1,36)+vac,
 '15_u7_m6':F(6,16)+F(18,225)+F(1,7)+F(12,36)+vac,
 '15_u7_absent_m13':F(13,16)+F(18,225)+F(7,49**2)+F(3,64)+F(2,36)+vac,
 '15_u7_present_m9':F(9,16)+F(18,225)+F(1,7)+F(2,64)+F(1,56**2)+F(6,36)+vac,
 '15_u7_present_m10_even_witness':F(10,16)+F(18,225)+F(1,7)+F(2,64)+F(1,56**2)+F(4,36)+F(1,21**2)+vac,
 '15_u7_present_m11':F(11,16)+F(18,225)+F(1,7)+F(6,42**2)+F(1,36)+vac,
 '15_u7_present_m12':F(12,16)+F(18,225)+F(1,7)+F(6,42**2)+vac,
 'dyadic7_with5':F(8,16)+F(24,196)+F(1,5)+F(4,400)+F(2,36)+vac,
 'dyadic7_without5_m10':F(10,16)+F(24,196)+F(9,36)+vac,
 'twenty_active5_high_four':F(20,625)+F(10,16)+F(3,64)+F(7,49)+F(3,36)+vac,
 'eighteen_active5_m11':F(18,625)+F(7,49)+F(3,36)+F(11,16)+F(3,125**2)+vac,
}
assert all(v<1 for v in checks.values()), checks
assert 1+13+18+11+3==46
assert 1+24+6+11+3==45
out={k:str(v) for k,v in checks.items()}
Path(__file__).with_name('isolation44-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
