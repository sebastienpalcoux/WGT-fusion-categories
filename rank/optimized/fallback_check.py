"""Exercise the production subtree handoff with an artificially small native range."""
from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run
import json

here = Path(__file__).resolve().parent
key = '4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,9,9,24,40,120,150'
with TemporaryDirectory() as temp:
    tmp = Path(temp)
    source = (here/'fast.cpp').read_text()
    old = 'if (L>(U(1)<<58)) {delegate(q); return;}'
    assert source.count(old) == 1
    (tmp/'small.cpp').write_text(source.replace(old,'if (L>U(64)) {delegate(q); return;}'))
    run(['c++','-O2','-std=c++17',str(tmp/'small.cpp'),'-o',str(tmp/'small')],check=True)
    run(['c++','-O2','-std=c++17',str(here/'exact.cpp'),'-lgmpxx','-lgmp','-o',str(tmp/'exact')],check=True)
    first = run([tmp/'small','28','8','1',key],capture_output=True,text=True,check=True)
    prefixes=[line.split()[2] for line in first.stderr.splitlines() if line.startswith('GMP ')]
    assert prefixes and all((p+',').startswith(key+',') for p in prefixes)
    merged=first.stdout.splitlines()
    for prefix in prefixes:
        part=run([tmp/'exact','--shard','28','1',prefix],capture_output=True,text=True,check=True)
        assert part.stderr.startswith('COMPLETE search rank=28 ')
        merged.extend(part.stdout.splitlines())
    whole=run([tmp/'exact','--shard','28','1',key],capture_output=True,text=True,check=True)
    assert merged and sorted(merged)==sorted(whole.stdout.splitlines()) and len(merged)==len(set(merged))
    print(f'{len(prefixes)} checked-arithmetic handoffs, {len(merged)} exact candidates')
