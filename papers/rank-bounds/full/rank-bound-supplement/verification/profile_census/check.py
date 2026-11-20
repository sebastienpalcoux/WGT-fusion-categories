"""Count odd-square starts and reindex the archived reservations after field-population tests."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt, lcm, prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PROVENANCE = 'verification/neutral_orbits/provenance.json'
PROVENANCE_SHA256 = '734cb5cf623e7b1b1bebfefbce4ff2285c5a1b256c86de3b8987cffb9c6bb9d9'


def read_verified(path, digest):
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ValueError(f'Digest mismatch at {path.name}')
    return raw


def possible_class(p, choice):
    valuation, multiplicity, active = choice
    if not valuation:
        return choice == (0,0,0)
    j, remainder = divmod(active, p-1)
    return (p >= 11 and valuation == 1 and remainder == 0 and 0 < j < p-1
            and (p-1) % j == 0 and j % ((p-1) & -(p-1)) == 0 and multiplicity == p*j)


def fits(rank, support, choices, mixed, nines=None):
    marked = {p:c for p,c in zip((q for q in support if q >= 5),choices) if c[0]}
    if mixed and not {3,5}.issubset(support):
        return False
    if len(support)==3 and marked:
        return False
    # A free slot has a composite denominator. The complete nine class is closed
    # only at the reservation stage. No other odd denominator between 9 and 15 occurs.
    floor = min(support)**2
    if nines is not None:
        floor = max(15, floor)
    minimum_m = prod(p**(2 if len(support)==3 or (p not in marked and not (mixed and p in (3,5))) else 1)
                     for p in support)
    if nines:
        minimum_m = lcm(minimum_m, 9)
    groups = [(1,minimum_m)]
    open_floors = []
    def book_field(total, modulus, exact, charges, extra_floor):
        rounded = ((total+modulus-1)//modulus)*modulus
        groups.extend(charges)
        charged = sum(n for n,_ in charges)
        if charged > rounded:
            raise ArithmeticError((charges,rounded))
        groups.append((rounded-charged, extra_floor))
        if not exact:
            open_floors.append(extra_floor)
    if mixed:
        s = nines or 0
        book_field(80+s,8,False,[(80,15),(s,9)],floor)
    for p in support:
        if mixed and p in (3,5):
            continue
        d = (p-1) & -(p-1)
        modulus = max(8,2*d)
        extra_floor = max(floor,p*min(support))
        if p in marked:
            _,m,a = marked[p]
            minimal = a == 2*(p-1)
            book_field(2*p+2 if minimal else m,modulus,minimal,[(m,p)],extra_floor)
            continue
        # Sparse absent fields have p^2 | every denominator. Dense activity
        # provides at least 2p(p-1) pure labels, more than the core d(p+1).
        if rank-1 < 2*p*(p-1):
            extra_floor = max(extra_floor,p*p)
        core = d*(p+1)
        if p == 3 and nines is not None:
            book_field(max(core,nines),modulus,False,[(nines,9)],extra_floor)
        else:
            book_field(core,modulus,False,[(core,p*p)],extra_floor)
    free = rank-sum(n for n,_ in groups)
    if free < 0:
        return False
    # In the all-exact case the actual rank could be smaller than the ceiling.
    if open_floors:
        groups.append((free,min(open_floors)))
    common = lcm(*(b for n,b in groups if n))
    return sum(n*(common//b)**2 for n,b in groups) >= common**2


def run(root, regenerate=False):
    manifest = json.loads(read_verified(root/PROVENANCE,PROVENANCE_SHA256))
    for path,digest in manifest['baseline_source_sha256'].items():
        read_verified(root/path,digest)
    sys.path.insert(0,str(root/'arithmetic'))
    import odd_square as historical
    reports = []
    for rank in (201,209):
        historical.RANK = rank
        primes = tuple(p for p in range(3,(rank+1)//2+1,2) if all(p%q for q in range(3,isqrt(p)+1,2)))
        profile_counts = Counter()
        starts = {}
        for size in range(3,len(primes)+1):
            if 2*sum(primes[:size])-1 > rank:
                break
            for support in combinations(primes,size):
                if 2*sum(support)-1 > rank:
                    continue
                allowed = []
                hist_count = 1
                for p in support:
                    if p<5:
                        continue
                    original = tuple(historical.profiles(p,support))
                    hist_count *= len(original)
                    allowed.append([c for c in original if possible_class(p,c)])
                profile_counts['archived'] += hist_count
                for choices in product(*allowed):
                    profile_counts['packet'] += 1
                    graphs = tuple(mixed for mixed in (False,True) if fits(rank,support,choices,mixed))
                    if not graphs:
                        continue
                    profile_counts['open_nines_components'] += 1
                    ns = range(0,77,2) if 3 in support else (0,)
                    by_graph = {g:tuple(n for n in ns if fits(rank,support,choices,g,n)) for g in graphs}
                    if any(by_graph.values()):
                        starts[support,choices] = {g:ns for g,ns in by_graph.items() if ns}
                        profile_counts['component_mass'] += 1
        entry = manifest['inputs'][str(rank)]
        raw = json.loads(read_verified(root/entry['path'],entry['sha256']))
        rows = raw if rank==201 else raw['rows']
        if profile_counts['archived'] != entry['prime_profiles'] or len(rows) != entry['reservations']:
            raise ValueError('Historical count mismatch')
        old = []; new = []; packets = []; rejected = Counter()
        for i,(p,c,g,s) in enumerate(rows):
            key = tuple(p),tuple(map(tuple,c))
            if all(possible_class(p,choice) for p,choice in zip((t for t in key[0] if t>=5),key[1])):
                packets.append(i)
            mixed = bool(g)
            if g and g != [[3,5]]:
                raise ValueError('Unexpected archived graph')
            if mixed not in starts.get(key,()):
                rejected['front-end'] += 1
            else:
                old.append(i)
                if s in starts[key][mixed]:
                    new.append(i)
                else:
                    rejected['nine-class'] += 1
        if regenerate:
            rebuilt = []
            for (support,choices),by_graph in starts.items():
                for graph,nines in historical.graph_profiles(support,choices):
                    if graph not in ((), ((3,5),)):
                        continue
                    if nines in by_graph.get(bool(graph), ()):
                        rebuilt.append(json.loads(json.dumps((support,choices,graph,nines))))
            if rebuilt != [rows[i] for i in new]:
                raise ValueError(f'Regenerated reservations differ at rank {rank}')
        reports.append({'rank':rank,'prime_profiles':profile_counts,'original_rows':len(rows),
                        'class_only_rows':len(packets),'before_nines_indices':old,
                        'retained_indices':new,'discarded':rejected,
                        'retained_edgeless':sum(not rows[i][2] for i in new),
                        'retained_mixed':sum(bool(rows[i][2]) for i in new),
                        'archive':entry['path'],'archive_sha256':entry['sha256']})
    return {'schema':'odd-square-profile-census-v1','source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'baseline_source_manifest_sha256':PROVENANCE_SHA256,'screens':reports}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplement',type=Path,default=HERE.parents[1])
    parser.add_argument('--write-receipt',action='store_true')
    parser.add_argument('--regenerate-graphs',action='store_true',help='also generate graph reservations from the checked prime starts')
    args=parser.parse_args()
    data=run(args.supplement.resolve(), args.regenerate_graphs)
    receipt=HERE/'receipt.json'
    if args.write_receipt:
        receipt.write_text(json.dumps(data,indent=1,sort_keys=True)+'\n')
    elif data!=json.loads(receipt.read_text()):
        raise ValueError('Profile census receipt mismatch')
    for screen in data['screens']:
        print(screen['rank'],dict(screen['prime_profiles']),screen['original_rows'],len(screen['before_nines_indices']),
              len(screen['retained_indices']),screen['retained_edgeless'],screen['retained_mixed'])

if __name__=='__main__':
    main()
