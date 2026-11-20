"""Re-enter the lattice bounds independently and compare every archived decision."""
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
from zipfile import ZipFile
import argparse
import gzip
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'census'))
from degree_slices import slice_failure

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)


def independently_fails(row, dimension, primes):
    by_degree = Counter(row[1:])
    for p in primes:
        eligible = []
        def pure(n):
            while n % p == 0:
                n //= p
            return n == 1
        with_pure = any(pure(dimension//(y*y)) for y in by_degree)
        for y,n in by_degree.items():
            f = dimension//(y*y)
            all_required = pure(f) or (with_pure and y%p != 0)
            admissible = []
            if not all_required:
                admissible.append(0)
            if f % (p*p) == 0:
                for a in ([n] if all_required else range(1,n+1)):
                    energy = a*(f-a)
                    cells=[]
                    for x,nx in by_degree.items():
                        g=gcd(x,y)
                        if g==1 or (pure(f) and dimension//(x*x)%(p*p)):
                            if x%p:
                                energy=-1
                                break
                            continue
                        u=x//g
                        spacing=p if u%p else 1
                        least=(a*g+spacing//2)%spacing-spacing//2
                        if abs(least)>=a*g:
                            energy=-1
                            break
                        energy-=nx*(u*least)**2
                        cells.append((x,nx,g,u,spacing,least))
                    if energy<0:continue
                    anchor=0
                    div=0
                    for x,nx,g,u,s,z in cells:
                        # All lifts in the per-coordinate ball are computed
                        # from its endpoints, rather than adjacent-lift costs.
                        cap=min(a*g-1,isqrt(energy+(u*z)**2)//u)
                        kmin=-((cap+z)//s)
                        kmax=(cap-z)//s
                        if kmin<kmax:div=gcd(div,x*u*s)
                        anchor+=x*u*nx*z
                    if ((anchor+a)%div==0 if div else anchor+a==0):
                        admissible.append(a)
            if not admissible:return True
            if max(admissible):eligible.append((y,admissible))
        if gcd(*(y for y,a in eligible))!=1:return True
        h=1
        while dimension%(h*p)==0:h*=p
        if h>200000:continue
        # Shift each cardinality from the unrotated source, using a
        # doubled periodic bit vector, and start at the vacuum residue.
        mask=(1<<h)-1;bits=2 % (1<<h)
        for y,choices in eligible:
            doubled=(bits<<h)|bits
            nbits=0
            for a in choices:
                nbits|=(doubled>>((-a*y*y)%h))&mask
            bits=nbits
        if not bits&1:return True
    return False


def replay(archive):
    counts={}
    with ZipFile(archive) as z:
        for r in range(1,45):
            section=f'ranks/{r:02}/'
            tasks=json.loads(z.read(section+'manifest.json'))['tasks']
            stats=Counter()
            for record in gzip.decompress(z.read(section+'raw-census.jsonl.gz')).splitlines():
                rank,tid,origin,reason,ds=json.loads(record)
                if r != rank:raise ValueError('rank mismatch')
                task=tasks[tid];M=task['B']*ds[-1];D=task['q']*M*M
                ps=[p for p in PRIMES if M%p==0]
                rem=M
                for p in ps:
                    while rem%p==0:rem//=p
                if rem!=1:raise ValueError('unknown prime support')
                outcome=slice_failure(ds,D,ps)
                if bool(outcome)!=independently_fails(ds,D,ps):
                    raise ValueError(('independent decision differs',r,tid,ds))
                stats['profiles']+=1;stats['pass']+=outcome is None
                if outcome:stats[outcome[1]]+=1
            counts[str(r)]=dict(stats)
    return dict(archive_sha256=sha256(Path(archive).read_bytes()).hexdigest(),ranks=counts)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    result = replay(args.archive)
    table_path = Path(__file__).with_name('table1.json')
    table = json.loads(table_path.read_text())
    if (table['archive_sha256'] != result['archive_sha256'] or
            any(row['degree_slices'] != result['ranks'][str(row['rank'])].get('pass', 0)
                for row in table['table'])):
        raise ValueError('Independent degree-slice totals differ from the seven-test comparison')
    print(json.dumps({rank: counts for rank, counts in result['ranks'].items() if counts}))
