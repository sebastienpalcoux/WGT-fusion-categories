"""Check saved odd-square reservations against neutral-orbit constraints.

This is a short exact necessary-condition screen. It runs no profile
regeneration and no final Gauss searches. Retained reservations are not
realizable modular data. The rank-209 input is exploratory.

Run from the supplement root:
    python3 verification/neutral_orbits/check.py
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt, prod
from pathlib import Path

SCHEMA = "odd-neutral-orbit-screen-v1"
HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parent.parent

def twopart(n):
    return n & -n


def screen(row, rank):
    support, classes, edges, nines = row
    marked = {p: tuple(v) for p, v in zip((p for p in support if p >= 5), classes) if v[0]}
    edges = tuple(map(tuple, edges))
    if edges not in ((), ((3, 5),)):
        return {'reason': 'unhandled_graph', 'keep': True}
    isolated = [p for p in support if all(p not in e for e in edges)]
    base = [(1, prod(support)), (nines, 9)]
    base += [(v[1], p) for p, v in marked.items()]
    # Pure active sectors of absent primes: odd:isolated-dyadic-cost.
    active = {p: twopart(p-1)*(p-1) for p in isolated if p != 3 and p not in marked}
    base += [(n, p*p) for p, n in active.items()]
    if edges:
        base.append((80, 15))  # odd:edge35-reservation, below triple activity.

    def mass(records):
        slots = rank-sum(n for n, b in records)
        return slots, sum((Q(n, b*b) for n, b in records), Q(0))+Q(slots, 225)

    def fits(records):
        slots, upper = mass(records)
        return slots >= 0 and upper >= 1

    if not fits(base):
        return {'reason': 'active_cost', 'keep': False, 'upper': str(mass(base)[1])}
    fields = []
    pure = []
    sparse = []
    for p in isolated:
        if p == 3:
            continue
        d = twopart(p-1)
        if p in marked:
            val, total, act = marked[p]
            if val != 1:
                continue
            neutral = total-act
            if neutral <= 0:
                continue  # No inference from an unverified pure-neutral hypothesis.
            if act < d*(p-1):
                return {'reason': 'present_active_cost', 'keep': False, 'prime': p}
            neutral_total = max(2*d, ((neutral+2*d-1)//(2*d))*(2*d))
            if p % 4 == 3 and (act//(p-1)) % 4 == 0:
                neutral_total = max(neutral_total, 8)
            fields.append((neutral_total-neutral, 3*p))
            sparse.append(p)
        else:
            # Even with these many p-active labels, ordinary mass or slots fail.
            # Thus actual A_p < p(p-1), guaranteeing an inactive pure p orbit.
            dense = base+[(max(0, p*(p-1)-active[p]), p*p)]
            if not fits(dense):
                fields.append((2*d, p*p))
                pure.append((d, p*p))
                sparse.append(p)
    # By field localization, these moved neutral sets are mutually disjoint.
    if not fits(base+fields):
        return {'reason': 'neutral_fields', 'keep': False,
                'sparse_primes': sparse, 'upper': str(mass(base+fields)[1])}
    # The pure neutral orbits are disjoint from all mixed-component neutrals.
    components = [{3, 5}]+[{p} for p in support if p not in (3, 5)] if edges else [{p} for p in support]
    c = len(components)
    mixed_floor = min(p*q for i, C in enumerate(components) for E in components[i+1:] for p in C for q in E)
    mixed = (2*(c-1), max(15, mixed_floor))
    if not fits(base+pure+[mixed]):
        return {'reason': 'mixed_neutral_orbits', 'keep': False,
                'sparse_primes': sparse, 'upper': str(mass(base+pure+[mixed])[1])}
    return {'reason': 'retained', 'keep': True, 'sparse_primes': sparse}



def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def validate_row(row, rank):
    require(isinstance(row, list) and len(row) == 4, "invalid reservation record")
    support, classes, edges, nines = row
    require(isinstance(support, list) and len(support) >= 3,
            "at least three dimension primes are required")
    require(support == sorted(set(support)) and all(type(p) is int and p >= 3
            and p % 2 and all(p % d for d in range(3, isqrt(p) + 1, 2))
            for p in support), "invalid odd prime support")
    require(len(classes) == len([p for p in support if p >= 5]),
            "prime-class metadata has the wrong length")
    for p, value in zip((p for p in support if p >= 5), classes):
        require(isinstance(value, list) and len(value) == 3
                and all(type(n) is int and n >= 0 for n in value),
                "invalid prime-class triple")
        valuation, total, active = value
        require((value == [0, 0, 0]) or (valuation in (1, 2)
                and 0 < total < rank and 0 <= active <= total),
                "invalid class counts or valuation marker")
        if valuation == 1:
            require(active > 0 and active % (p - 1) == 0,
                    "invalid valuation-one activity")
    require(type(nines) is int and 0 <= nines < rank and nines % 2 == 0
            and (3 in support or nines == 0), "invalid denominator-nine count")
    require(isinstance(edges, list), "invalid twist-graph edge list")
    require(all(isinstance(edge, list) and len(edge) == 2
                and edge[0] < edge[1] and all(p in support for p in edge)
                for edge in edges), "invalid twist-graph edge")
    require(len(set(map(tuple, edges))) == len(edges), "duplicate twist-graph edge")


def evaluate(rank, root, manifest):
    metadata = manifest["inputs"][str(rank)]
    path = (root / metadata["path"]).resolve()
    require(path.is_relative_to(root), "input path leaves the supplement")
    require(file_hash(path) == metadata["sha256"], "saved input fingerprint changed")
    original = json.loads(path.read_text())
    if rank == 201:
        rows = original
    else:
        require(original["count"] == metadata["prime_profiles"],
                "exploratory profile count changed")
        rows = original["rows"]
    require(len(rows) == metadata["reservations"], "saved reservation count changed")
    counts = Counter()
    outcomes = []
    retained = []
    for index, row in enumerate(rows):
        validate_row(row, rank)
        verdict = screen(row, rank)
        counts[verdict["reason"]] += 1
        outcomes.append({"index": index, **verdict})
        if verdict["keep"]:
            retained.append(index)
    return {
        "schema": SCHEMA,
        "scope": ("Necessary-condition screen of the specified saved arithmetic "
                  "reservations. Retention does not establish realizability; "
                  "the rank-209 input is exploratory, and no rank-209 exclusion "
                  "or complete modular-data census is certified."),
        "rank": rank,
        "baseline_prime_profiles": metadata["prime_profiles"],
        "input_reservations": len(rows),
        "input_path": metadata["path"],
        "input_sha256": file_hash(path),
        "checker_sha256": file_hash(Path(__file__).resolve()),
        "provenance_sha256": file_hash(HERE / "provenance.json"),
        "baseline_receipt_sha256": manifest["baseline_receipt_sha256"],
        "counts": dict(counts),
        "retained": len(retained),
        "retained_indices": retained,
        "outcomes": outcomes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rank", choices=("201", "209", "all"), default="all")
    parser.add_argument("--supplement-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--write-receipts", action="store_true",
                        help="replace only this checker's two derived receipts")
    args = parser.parse_args()
    root = args.supplement_root.resolve()
    manifest = json.loads((HERE / "provenance.json").read_text())
    require(manifest["schema"] == "odd-neutral-screen-inputs-v1", "unknown input schema")
    require(file_hash(root / manifest["baseline_receipt"])
            == manifest["baseline_receipt_sha256"], "baseline receipt changed")
    for relative, expected in manifest["baseline_source_sha256"].items():
        require(file_hash(root / relative) == expected,
                "baseline source changed: " + relative)
    ranks = (201, 209) if args.rank == "all" else (int(args.rank),)
    summaries = []
    for rank in ranks:
        result = evaluate(rank, root, manifest)
        destination = HERE / ("rank%d-receipt.json" % rank)
        if args.write_receipts:
            destination.write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
        else:
            require(json.loads(destination.read_text()) == result,
                    "neutral-orbit receipt mismatch at rank %d" % rank)
        summaries.append({key: result[key] for key in
                          ("rank", "input_reservations", "counts", "retained")})
    print(json.dumps({"status": "passed", "screens": summaries,
                      "final_gauss_searches": 0}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
