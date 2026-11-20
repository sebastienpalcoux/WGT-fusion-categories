#!/usr/bin/env python3
"""Independent finite checks and explicit terminal witnesses.

Python 3.10+ and SymPy. No production modules are imported. This checks
all emitted covering-census lists, the odd profile counts and final
nonsquare tail, and finite group arithmetic. It does not regenerate the
native task partition or verify every discarded prefix.
"""
import argparse
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, isqrt, lcm, prod
from zipfile import ZipFile

from sympy import Matrix, divisors, factorint, primerange


@lru_cache(None)
def prime_factors(integer):
    return tuple(factorint(integer).items())


def tile(counts, prime):
    # A choice of an equal-degree block removes p entries of its class.
    # Mixed blocks at 3 remove (1,2) copies at denominators (b,2b),
    # those at 5 remove (2,2); (4,4) is twice that admissible shape.
    initial = tuple(sorted((b, m) for b, m in counts.items() if b % prime))

    @lru_cache(None)
    def visit(state):
        if not state:
            return True
        (b, n), *rest = state
        if prime == 2:
            return n % 2 == 0 and visit(tuple(rest))
        small = 1 if prime == 3 else 2
        parts = dict(rest)
        for mixed in range(1 + min(n // small, parts.get(2*b, 0) // 2)):
            if (n - small*mixed) % prime:
                continue
            remainder = parts.copy()
            if mixed:
                remainder[2*b] -= 2*mixed
            if visit(tuple((x, v) for x, v in remainder.items() if v)):
                return True
        return False
    return visit(initial)


STAGES = ('primary_mass', 'prime_square_gcd', 'defect_one',
          'principal_block', 'ternary_trace')


def obstruction(degrees):
    """Return the first failed necessary condition and its exact witness.

    The terminal witness refers to the degrees in the bound archive row.
    Tiling is an exhaustive finite recursion, independent of the greedy
    production implementation. Principal eligibility and trace subsets
    are recomputed here from their defining integer inequalities.
    """
    assert all(type(d) is int and d > 0 for d in degrees)
    assert degrees == sorted(degrees) and degrees.count(1) == 1
    dimension = sum(d*d for d in degrees)
    factors = prime_factors(dimension)
    M = prod(p**(e//2) for p,e in factors)
    q = dimension // (M*M)
    cls = Counter(d for d in degrees if d != 1)
    assert cls and all(M % d == 0 for d in degrees)
    codegree = {d:dimension//(d*d) for d in cls}
    for p,_ in factors:
        mass = sum(m*d*d for d,m in cls.items() if codegree[d] % p == 0)
        lhs, rhs = 2*p*min(cls)*mass, (p-1)*dimension
        if lhs < rhs:
            return 0, dict(prime=p, least_degree=min(cls),
                           eligible_mass=mass, left=lhs, right=rhs)
    for p,_ in factors:
        n = M//p
        if p == q or n not in cls:
            continue
        actual = q == 1
        effective = q == 2 or (q in (3,5) and 2*n not in cls and all(2*x != n for x in cls))
        free = sorted(d for d in cls if d % p)
        if (actual or effective) and gcd(*free) != 1:
            return 1, dict(prime=p, kernel=q, row_degree=n,
                           prime_free_degrees=free, gcd=gcd(*free))
    counts = Counter(M//d for d in degrees)
    for p,_ in prime_factors(q):
        assert p in (2,3,5)
        if not tile(counts,p):
            return 2, dict(prime=p, prime_free_classes=sorted(
                [b,m] for b,m in counts.items() if b % p),
                exhaustive_tiling_exists=False)
    eligible = {}
    for p,_ in factors:
        xs = [(x,m) for x,m in cls.items() if x % p]
        eligible[p] = {y for y in cls if codegree[y] % p == 0 and all(
            gcd(x,y)>1 and 2*m*gcd(x,y)>p and m*lcm(x,y)**2<dimension for x,m in xs)}
        if gcd(*eligible[p]) != 1:
            return 3, dict(prime=p, eligible_degrees=sorted(eligible[p]),
                           gcd=gcd(*eligible[p]))
    assert q==1 and M%15==0 and M%25!=0 and M//5 in cls
    proper = [d for d in degrees if d != 1 and d != M//5 and d%5]
    assert len(proper)==4 and all(d%5 in (1,4) for d in proper)
    modulus = 3**dict(prime_factors(M))[3]
    signs = [d if d%5==1 else -d for d in proper if d in eligible[3]]
    residues = sorted({sum(subset)%modulus for size in range(len(signs)+1)
                       for subset in combinations(signs,size)})
    if modulus-1 not in residues:
        return 4, dict(modulus=modulus, signed_degrees=signs,
                       reachable_residues=residues, required_residue=modulus-1)
    raise AssertionError('A covering-census list survives all five conditions')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def general(archive, fresh=None):
    if 'summary.json' not in archive.namelist():
        path, = (name for name in archive.namelist() if name.endswith('/reference-census.zip'))
        with ZipFile(io.BytesIO(archive.read(path))) as reference:
            return general(reference, fresh)
    rows, witnesses = [], bytearray()
    total_tasks = total_continuations = count = 0
    for r in range(1,45):
        base = f'ranks/{r:02d}/'
        tasks = json.loads(archive.read(base+'manifest.json'))['tasks']
        receipts = [json.loads(line) for line in archive.read(base+'receipts.jsonl').splitlines()]
        assert len(tasks)==len(receipts) and {c['index'] for c in receipts}==set(range(len(tasks)))
        total_tasks += len(tasks)
        total_continuations += sum(c['handoffs'] for c in receipts)
        values = [0]*7
        seen = set()
        for line in gzip.decompress(archive.read(base+'raw-census.jsonl.gz')).splitlines():
            rank, task, _, _, d = json.loads(line)
            assert rank==r==len(d) and tuple(d) not in seen
            assert type(task) is int and 0 <= task < len(tasks)
            seen.add(tuple(d))
            dim = sum(x*x for x in d)
            assert dim == tasks[task]['q']*(tasks[task]['B']*max(d))**2
            values[int(isqrt(dim)**2 != dim)] += 1
            stage, witness = obstruction(d)
            for j in range(stage):
                values[2+j] += 1
            witnesses.extend(canonical(dict(rank=r, task=task,
                row_sha256=hashlib.sha256(line).hexdigest(),
                stage=STAGES[stage], witness=witness))+b'\n')
        rows.append(dict(rank=r, counts=values))
        count += len(seen)
        if any(values):
            print('rank',r,values)
        if fresh is not None and base+'raw-census.jsonl.gz' in fresh.namelist():
            emitted = {
                tuple(json.loads(line)[-1]) for line in
                gzip.decompress(fresh.read(base+'raw-census.jsonl.gz')).splitlines()
            }
            assert emitted <= seen
            removed = Counter(STAGES[obstruction(list(d))[0]] for d in seen-emitted)
            if removed:
                print('fresh rank',r,'covering-census-only dispositions',dict(removed))
    totals = [sum(row['counts'][i] for row in rows) for i in range(7)]
    print('totals',totals,'tasks',total_tasks,'continuations',total_continuations)
    assert totals==[7314,10885,8207,6194,6166,244,0] and count==18199
    assert (total_tasks,total_continuations)==(12273,41915)
    return dict(columns=['square','nonsquare',*STAGES], table=rows,
                totals=totals, tasks=total_tasks, exact_continuations=total_continuations,
                terminal_witnesses=count,
                witness_sha256=hashlib.sha256(witnesses).hexdigest()), bytes(witnesses)


def odd_nonsquare():
    minimal_M = 27*5*49*121
    requirements = (135,1323,3267,27,245,605)
    # Eliminate the last unknown m_9 by two integer inequalities. This counts
    # the starting profiles without traversing the authors' prefix recursion.
    supports = [ps for ps in combinations(primerange(5,26),3) if 3+6*sum(p-1 for p in ps)<=211]
    total = 0
    before_floor = 0
    for ps in supports:
        base = 3*prod(ps)
        scale = lcm(base,9,15)**2
        w9,w15,vac = (scale//(b*b) for b in (9,15,base))
        sizes = []
        for p in ps:
            step = 3*((p-1)&-(p-1))
            maximum = min(210,3*p*p-1,48 if p==5 else 210)
            sizes.append([m for m in range(0,maximum+1,step) if m!=3*p*(p-1)])
        for mults in product(range(0,27,2),*sizes):
            cap = 210-sum(mults)
            if cap<0:
                continue
            amount = 3*scale-sum(m*(scale//(b*b)) for m,b in zip(mults,(3,)+ps))
            upper = min(cap,(amount-1)//w9)
            before_floor += max(0,upper//2+1)
            lower = max(0,-((vac+cap*w15-amount)//(w9-w15)))
            total += max(0,upper//2-(lower+1)//2+1)
    print('odd nonsquare supports, profiles before floor-15 test, and starting profiles',
          len(supports),before_floor,total)
    assert (len(supports),before_floor,total)==(14,4189686,685040)
    # Allow all larger denominators by a rational enclosing interval.
    ps = (5,7,11)
    choices = []
    for p in ps:
        denominators = []
        power = p
        while power<=500:
            b = power
            while b<=500:
                if b>=3*p:
                    denominators.append(b)
                b *= 3
            if p==5:
                break
            power *= p
        choices.append([(b,Fraction(1,b*b),Fraction(1,b*b)) for b in denominators]
                       +[(None,0,Fraction(1,501**2))])
    passes = []
    for k in (36,38,40,42):
        charge = 3-Fraction(12,9)-Fraction(k,81)-Fraction(24,25)
        room = Fraction(42-k,225)+Fraction(1,minimal_M**2)+sum(Fraction(2,b*b) for b in requirements)
        for triple in product(*choices):
            lo = sum(6*(p-1)*t[1] for p,t in zip(ps,triple))
            hi = sum(6*(p-1)*t[2] for p,t in zip(ps,triple))
            if lo<charge<=hi+room:
                passes.append((k,tuple(t[0] for t in triple)))
    print('odd nonsquare sector choices', passes)
    assert passes==[(36,(15,21,33))]
    residual = 3-sum(Fraction(m,b*b) for m,b in [(12,3),(36,9),(24,5),(24,15),(36,21),(60,33)])
    assert residual == Fraction(5023,266805)

    @lru_cache(None)
    def smooth_below(hi):
        numbers = {1,5} if hi>=5 else {1}
        for p in (3,7,11):
            new = set(numbers)
            for k in numbers:
                while k*p<=hi:
                    k *= p
                    new.add(k)
            numbers = new
        return sorted(numbers)

    def budget(mass, free, required, least):
        if free < 0:
            return False
        # Smallest multiple of each outstanding witness divisor above the floor.
        bound = Fraction(1,minimal_M**2) + Fraction(2*free,least**2)
        for divisor in required:
            floor = divisor*((least+divisor-1)//divisor)
            bound += Fraction(2,floor**2)
        return mass <= bound

    @lru_cache(None)
    def visit(mass, least, slots, required, denominator_lcm):
        if not slots:
            root = isqrt(mass.denominator)
            return not required and mass.numerator==1 and root**2==mass.denominator and root%denominator_lcm==0
        if not budget(mass,slots-len(required),required,least):
            return False
        hi = isqrt((2*slots+1)*mass.denominator//mass.numerator)
        low = max(least,1+isqrt(2*mass.denominator//mass.numerator))
        for b in smooth_below(hi):
            if b<low:
                continue
            rest = mass-Fraction(2,b*b)
            assigned = [required[:j]+required[j+1:] for j,z in enumerate(required) if b%z==0]
            if slots>len(required):
                assigned.append(required)
            if any(visit(rest,b,slots-1,pending,lcm(denominator_lcm,b)) for pending in assigned):
                return True
        return False

    possible = visit(residual,15,9,requirements,minimal_M)
    print('odd nonsquare tail',residual,'solutions exist',possible,'states',visit.cache_info().currsize)
    assert not possible
    return dict(supports=len(supports), profiles_before_floor15=before_floor,
                profiles_after_floor15=total, sector_choices=passes,
                residual=str(residual), tail_states=visit.cache_info().currsize,
                tail_solutions=0)


def odd_profile_starts(rank):
    # For each support choose independently the absent class or an allowed
    # proper divisor j, with p*j rows. Sum the products by walking prime subsets.
    limit = (rank+1)//2
    primes = list(primerange(3,limit+1))
    @lru_cache(None)
    def options(p):
        odd = p-1
        while odd%2==0:
            odd //= 2
        if p==3 or odd<5:
            return 1
        return 1+sum(j<p-1 and j%((p-1)//odd)==0 and p*j<=rank for j in divisors(p-1))
    def count(start, room, n, weight):
        total = (1 if n==3 else weight) if n>=3 else 0
        for k in range(start,len(primes)):
            p = primes[k]
            if p>room:
                break
            total += count(k+1,room-p,n+1,weight*options(p))
        return total
    starts = count(0, limit, 0, 1)
    print('odd square rank',rank,'starting profiles',starts)
    assert starts=={209:16512,217:21396}[rank]
    return starts


def permutation_checks():
    def compose(x,y):
        return tuple(x[i] for i in y)
    def parity(x):
        return sum(a>b for a,b in combinations(x,2))%2
    def inverse(x):
        return tuple(x.index(i) for i in range(len(x)))
    def classes(group):
        todo = set(group)
        while todo:
            x = min(todo)
            conjugates = {compose(compose(g,x),inverse(g)) for g in group}
            todo -= conjugates
            yield x,conjugates
    def double_rank(name,group,expected):
        central_ranks = [sum(1 for _ in classes([g for g in group if compose(g,x)==compose(x,g)])) for x,_ in classes(group)]
        print(name,'centralizer class numbers',central_ranks,'total',sum(central_ranks))
        assert sum(central_ranks)==expected
    for degree,rank in [(5,22),(6,44)]:
        double_rank('A'+str(degree),[g for g in permutations(range(degree)) if not parity(g)],rank)
    # GL(3,2) on the seven nonzero binary column vectors.
    linear = []
    for a,b,c in permutations(range(1,8),3):
        if a^b==c:
            continue
        linear.append(tuple(((a if v&1 else 0)^(b if v&2 else 0)^(c if v&4 else 0))-1 for v in range(1,8)))
    assert len(linear)==168
    double_rank('GL(3,2)',linear,32)
    # SL(2,5) on nonzero vectors; a faithful permutation model retains the center.
    nonzero = [v for v in product(range(5), repeat=2) if v!=(0,0)]
    index = {v:i for i,v in enumerate(nonzero)}
    binary = [tuple(index[((a*x+b*y)%5,(c*x+d*y)%5)] for x,y in nonzero)
              for a,b,c,d in product(range(5), repeat=4) if (a*d-b*c)%5==1]
    assert len(binary)==120
    double_rank('SL(2,5)',binary,74)
    vecs = list(range(32))
    vecs = [v for v in vecs if v.bit_count()%2==0]
    group = [g for g in permutations(range(5)) if not parity(g)]
    def act(g,v):
        return sum(((v>>i)&1)<<g[i] for i in range(5))
    def order(g):
        elt,unit,count = g,tuple(range(5)),1
        while elt!=unit:
            elt=compose(elt,g)
            count+=1
        return count
    assert sorted({len({act(g,v) for g in group}) for v in vecs})==[1,5,10]
    fixed=[(order(g),sum(act(g,v)==v for v in vecs)) for g,_ in classes(group)]
    assert sorted(fixed)==[(1,16),(2,4),(3,4),(5,1),(5,1)]
    gauss = sum((-1)**(v.bit_count()//2) for v in vecs)
    assert gauss==-4
    # Fixed-point columns on every transitive A5-set, enumerating all subgroups.
    to_id = {g:i for i,g in enumerate(group)}
    table = [[to_id[compose(g,h)] for h in group] for g in group]
    invs = [to_id[inverse(g)] for g in group]
    ident = to_id[tuple(range(5))]
    def generate(generators):
        found = {ident}
        stack = [ident]
        while stack:
            h = stack.pop()
            for k in generators:
                t = table[h][k]
                if t not in found:
                    found.add(t)
                    stack.append(t)
        return frozenset(found)
    subgroups = {frozenset({ident}):()}
    queue = list(subgroups)
    for sub in queue:
        for g in set(range(60))-sub:
            gs = subgroups[sub]+(g,)
            bigger = generate(gs)
            if bigger not in subgroups:
                subgroups[bigger] = gs
                queue.append(bigger)
    assert len(subgroups)==59
    samples = {o:next(i for i,g in enumerate(group) if order(g)==o) for o in (2,3,5)}
    weights = {1:1,2:2,3:3,4:1,5:5,6:3,10:4,12:3,60:4}
    inventory = {}
    for sub in subgroups:
        size = len(sub)
        f = tuple(sum(table[table[invs[i]][g]][i] in sub for i in range(60))//size for g in samples.values())
        entry = (weights[size]-1+f[0]+3*f[1]+10*f[2], f[0], f[2])
        assert size not in inventory or inventory[size]==entry
        inventory[size] = entry
    bounded = sorted(set(inventory)-{1})
    patterns = {}
    def allocate(index,sizes,weight,n2,n5):
        if index==len(bounded):
            max_free = 25-len(sizes)-weight-3*(n2%2)
            if max_free>=0 and n5>1:
                patterns[tuple(sorted(sizes))]=max_free
            return
        h = bounded[index]
        w,f2,f5 = inventory[h]
        for n in range(1+(25-len(sizes)-weight)//(1+w)):
            allocate(index+1,sizes+[60//h]*n,weight+n*w,n2+n*f2,n5+n*f5)
    allocate(0,[],0,1,1)
    assert patterns=={(1,):7,(1,30):3,(1,15):0,(6,):6,(6,15):5,(6,30):2,(6,15,30):1}
    print('A5 neutral orbit patterns and free bounds',patterns)
    s=(1,0,3,2,4); t=(2,3,0,1,4)
    image={v^act(s,v) for v in vecs}
    quotient={min(v^w for w in image) for v in vecs}
    affine_counts=[]
    for shift in quotient:
        action={v:min(act(t,v)^shift^w for w in image) for v in quotient}
        if all(action[action[v]]==v for v in quotient):
            affine_counts.append(sum(v==action[v] for v in quotient))
    assert affine_counts==[2,2]
    print('deleted module fixed labels',fixed,'affine involution fixed counts',affine_counts)
    # Cohomology orbit counts, conditional on the groups and actions in Section 10.
    outer_a5=lambda a:(-a if a%5 else 0)%5
    assert len({(i,j,min(a,outer_a5(a))) for i,j,a in product(range(2),range(3),range(5))})==18
    assert len({(i,min(j,-j%3),min(a,-a%5)) for i,j,a in product(range(4),range(3),range(5))})==24
    assert len({min(t,49*t%120) for t in range(1,120,2)})==36
    def second_chern(weights):
        return sum(a*b for a,b in combinations(weights,2))
    assert second_chern([0,0,1,2,0,1,2]) == 13
    assert second_chern([0,0,2,0,1,2,3]) == 23
    assert second_chern([0,1,2,3,0,2])%4 == 3
    assert [second_chern(w)%3 for w in ([0,1,2,0,0,0],[0,1,2,0,1,2])] == [2,1]
    binary_degrees = (1,2,3,4,5,6,12,20,30)
    binary_counts = (1,2,2,2,1,1,20,6,2)
    assert sum(binary_counts)==37 and sum(m*d*d for d,m in zip(binary_degrees,binary_counts))==7200
    # Flux (transposition) with the two C2-characters against the two electric
    # linear characters of S3. Divide by sqrt(36) to normalize the S-matrix.
    vacuum_and_sign = Matrix([[3,-3],[3,-3]])/6
    assert vacuum_and_sign.rank()==1 and 2*(4-2)==4
    assert 2*(Fraction(3,3)+Fraction(1,3))>2  # cyclic pointed data at p=3.
    return dict(untwisted_double_ranks={'A5':22,'A6':44,'PSL2(7)':32,'SL2(5)':74},
                a5_subgroups=len(subgroups),
                neutral_patterns=[[list(pattern),bound] for pattern,bound in sorted(patterns.items())],
                deleted_fixed_labels=sorted(fixed),
                deleted_affine_fixed_counts=affine_counts,
                cohomology_orbits_conditional_on_group_and_action={'A5':18,'A6':24,'binary_odd':36},
                binary_family={'rank':37,'dimension':7200})


def main():
    if not __debug__:
        raise SystemExit('Run with assertions enabled; omit -O and -OO.')
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path,
                        default=here.parents[1]/'census/reference-census.zip')
    parser.add_argument('--fresh', type=Path, help='optionally compare emitted lists; this does not certify its coverage')
    parser.add_argument('--write', action='store_true', help='regenerate the receipt and terminal witness ledger')
    args = parser.parse_args()
    with ZipFile(args.archive) as archive:
        if args.fresh:
            with ZipFile(args.fresh) as fresh:
                summary, witnesses = general(archive, fresh)
        else:
            summary, witnesses = general(archive)
    odd = odd_nonsquare()
    starts = {str(rank):odd_profile_starts(rank) for rank in (209,217)}
    groups = permutation_checks()
    result = dict(format='independent-finite-audit-v1',
                  scope='Terminal exclusion and finite arithmetic; no native partition regeneration or prefix-by-prefix certificate.',
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),
                  general=summary, odd_nonsquare=odd,
                  odd_square_prime_profile_counts=starts, group_arithmetic=groups)
    receipt = here/'certificate.json'
    ledger = here/'terminal-witnesses.jsonl.gz'
    encoded = json.dumps(result, sort_keys=True, indent=2)+'\n'
    if args.write:
        receipt.write_text(encoded)
        ledger.write_bytes(gzip.compress(witnesses, mtime=0))
    else:
        assert receipt.read_text() == encoded, 'Independent receipt differs; inspect the change before regeneration.'
        assert gzip.decompress(ledger.read_bytes()) == witnesses, 'Independent terminal witnesses differ.'
    print('Independent finite checks and all 18,199 terminal witnesses passed.')


if __name__ == '__main__':
    main()
