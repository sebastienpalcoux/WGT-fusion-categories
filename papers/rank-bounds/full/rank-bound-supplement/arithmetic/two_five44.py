"""Exact short completion of all rank44 high-four two--five modes.
Imports only the existing omission-reservation calculation; no census.
"""
import sys,json
from pathlib import Path
from fractions import Fraction as F
from large_four import preliminary,profiles,MIXED

survivors=[]
for m,five,support,modes in preliminary(44):
    if 5 not in support or modes[support.index(5)]!=MIXED:
        continue
    out=list(profiles(44,m,five,support,modes))
    if out:
        assert support==(3,5,7) and modes==(0,1,2) and m in (10,11)
        survivors.append((m,out))
assert [(m,len(out)) for m,out in survivors]==[(10,5),(11,1)]
maxima={m:max(x[3] for x in out) for m,out in survivors}
assert maxima=={10:F(77237,75600),11:F(233293,226800)}
for m,out in survivors:
    for odd,even,other,mass in out:
        assert odd==(1,2,4)
        assert 5 in even and 6 in even
        if m==11:assert even==(3,5,6) and not other
    assert maxima[m]-1 < 2*(F(1,36)-F(1,100))
    assert maxima[m]-1 < 4*(F(1,100)-F(1,400))
    assert maxima[m]-1 < F(3,81)-F(3,729)
assert maxima[10]-1 < 2*(F(1,36)-F(1,81))
assert maxima[10]-1 < F(1,36)-F(1,729)

base=lambda m,n:F(m,16)+F(16,100)+F(n,81)+F(1,7)
assert base(11,3)>1
old_neutral=sum((F(1,b*b) for b in (72,40,56,252,140)),F(1,2520**2))
assert old_neutral==F(13,10800)
four_nine_residual=1-base(10,4)
assert F(1,36)>four_nine_residual>F(1,64)+old_neutral

residual=1-base(10,3)
assert residual==F(1327,37800)
gap=residual-F(1,36)-F(1,144)
assert gap==F(29,75600)
large_five=sum((F(1,b*b) for b in (216,200,56,756,700)),F(1,37800**2))
small_five=sum((F(1,b*b) for b in (864,160,224,756,140)),F(1,30240**2))
assert large_five==F(43949,119070000)<gap
assert small_five==F(2873,25401600)<gap
assert F(2,64)+max(large_five,small_five)<residual
assert F(2,36)>residual
assert F(1,36)+F(1,121)>residual
assert [a for a in range(3,9) if (135*pow(2,a,7))%7 in (1,6)]==[5,8]

out={'rank':44,'modes':[(m,len(v)) for m,v in survivors],
     'maximum_masses':{m:str(v) for m,v in maxima.items()},
     'three_nine_residual':str(residual),'four_nine_residual':str(four_nine_residual),
     'required_gap':str(gap),'neutral_upper_c_ge2':str(large_five),
     'neutral_upper_c1':str(small_five),'surviving_modes':0}
Path(__file__).with_name('two-five44-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
