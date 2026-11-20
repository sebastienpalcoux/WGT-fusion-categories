"""First open rank, through the disjoint component-population mass test."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RANK = 217


def run(root, reservations=False):
    census = root / "verification/profile_census/check.py"
    sys.path[:0] = [str(census.parent), str(root / "arithmetic")]
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("profile_census", census)
    original = module_from_spec(spec)
    spec.loader.exec_module(original)
    import odd_square
    provenance = json.loads(original.read_verified(root / original.PROVENANCE, original.PROVENANCE_SHA256))
    for filename, digest in provenance["baseline_source_sha256"].items():
        original.read_verified(root / filename, digest)
    odd_square.RANK = RANK
    primes = tuple(p for p in range(3,(RANK+1)//2+1,2) if all(p%k for k in range(3,isqrt(p)+1,2)))
    counts = Counter()
    component_counts = Counter()
    screened = Counter()
    for size in range(3,len(primes)+1):
        if 2*sum(primes[:size])-1>RANK:
            break
        for support in combinations(primes,size):
            if 2*sum(support)-1>RANK:
                continue
            opts = [tuple(odd_square.profiles(p,support)) for p in support if p>=5]
            for choices in product(*opts):
                counts["original_prime_starts"] += 1
                if not all(original.possible_class(p,c) for p,c in zip((p for p in support if p>=5),choices)):
                    continue
                counts["classified_starts"] += 1
                graphs = [g for g in (False,True) if original.fits(RANK,support,choices,g)]
                if not graphs:
                    continue
                counts["component_starts_open_nines"] += 1
                ns = range(0,77,2) if 3 in support else (0,)
                retained = False
                for mixed in graphs:
                    nines = [s for s in ns if original.fits(RANK,support,choices,mixed,s)]
                    if nines:
                        retained = True
                        component_counts["mixed" if mixed else "edgeless"] += len(nines)
                counts["component_starts_closed_nines"] += int(retained)
                if retained and reservations:
                    for graph, s in odd_square.graph_profiles(support,choices):
                        if graph not in ((),((3,5),)):
                            continue
                        mixed = bool(graph)
                        if mixed in graphs and original.fits(RANK,support,choices,mixed,s):
                            screened["mixed" if mixed else "edgeless"] += 1
    result = {"rank":RANK,"starts":counts,"component_reservations":component_counts,
            "sha256":{"verification/profile_census/check.py":sha256(census.read_bytes()).hexdigest(),
                      "check.py":sha256(Path(__file__).read_bytes()).hexdigest()},
            "baseline_provenance":original.PROVENANCE_SHA256}
    if reservations:
        result["graph_reservations"] = screened
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement",type=Path,default=HERE.parents[1])
    parser.add_argument("--write-receipt",action="store_true")
    parser.add_argument("--graphs",action="store_true",help="also apply the graph reservation tests to the 227 surviving starts")
    args=parser.parse_args()
    receipt=HERE/("graphs.json" if args.graphs else "receipt.json")
    result=run(args.supplement.resolve(),args.graphs)
    if args.write_receipt:
        receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    elif result!=json.loads(receipt.read_text()):
        raise SystemExit("Frontier count does not match receipt.")
    print(dict(result["starts"]),dict(result["component_reservations"]),dict(result.get("graph_reservations",{})))
