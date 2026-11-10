"""Project the septic and fourteen-dyadic mass theorems onto the original rank-30 leaves."""
import argparse
from fractions import Fraction
import json
from math import isqrt, lcm
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]


def prefix_flags(prefix, q, rank):
    if q != 1 or prefix.count(4) != 14:
        return False, False
    septic = 7 in prefix
    # At rank 30 the codegree-sixteen saturation bound permits at most 14 fours.
    # The complement must have mass <1/16; the pivot and the vacuum are 4-divisible.
    exterior = sum((Fraction(1, b*b) for b in prefix if b % 4), Fraction())
    if 16 * exterior >= 1:
        return septic, True
    # Bound all possible contributions to the >1/16 required interior mass,
    # including the vacuum with denominator M at least 5400 and divisible by L.
    interior = sum((Fraction(1, b*b) for b in prefix if b != 4 and b % 4 == 0), Fraction())
    next_b = 4 * ((max(prefix[-1], 5) + 3)//4)
    common = lcm(*prefix)
    min_root = common * ((5400 + common - 1)//common)
    h = rank - 1 - len(prefix)
    upper = interior + Fraction(h, next_b*next_b) + Fraction(1, min_root*min_root)
    return septic, 16 * upper <= 1


def row_flags(degrees, q):
    if q != 1:
        return False, False
    D = sum(d*d for d in degrees)
    M = isqrt(D)
    if M*M != D or any(M % d for d in degrees):
        raise ValueError("invalid square-dimension row")
    bs = [M//d for d in degrees]
    if bs.count(4) != 14:
        return False, False
    exterior = sum(d*d for d,b in zip(degrees,bs) if b % 4)
    # The list includes the vacuum. Both mass formulations must agree.
    interior = sum(d*d for d,b in zip(degrees,bs) if b != 4 and b % 4 == 0)
    if 16*(exterior+interior) != 2*D:
        raise ValueError("incorrect 14-four complement")
    return 7 in bs, 16*exterior >= D


def totals():
    return dict(leaves=0, rows=0, native_nodes=0, gmp_nodes=0, gmp_handoffs=0)


def accumulate(counter, leaf):
    counter["leaves"] += 1
    counter["rows"] += len(leaf["rows"])
    for field in ("native_nodes", "gmp_nodes"):
        counter[field] += leaf[field]
    counter["gmp_handoffs"] += leaf["fallbacks"]


def audit(path):
    names = ("septic", "dyadic_mass", "both", "union")
    rows = dict.fromkeys(names, 0)
    prefixes = {name: totals() for name in names}
    physical = totals()
    fourteen = totals()
    all_rows_hit = totals()
    with ZipFile(path) as z:
        spec = json.loads(z.read("spec.json"))
        if spec["rank"] != 30:
            raise ValueError("only the original rank-30 return is in scope")
        for line in z.open("receipts.jsonl"):
            parent = json.loads(line)
            for leaf in parent.get("children", [parent]):
                accumulate(physical, leaf)
                if leaf["q"] == 1 and leaf["prefix"].count(4) >= 14:
                    accumulate(fourteen, leaf)
                s, d = prefix_flags(leaf["prefix"], leaf["q"], spec["rank"])
                for name, hit in zip(names,(s,d,s and d,s or d)):
                    if hit:
                        accumulate(prefixes[name], leaf)
                all_hit = bool(leaf["rows"])
                for row in leaf["rows"]:
                    a,b = row_flags(row, leaf["q"])
                    if (s and not a) or (d and not b):
                        raise ValueError("prefix consequence failed on emitted row")
                    for name, hit in zip(names,(a,b,a and b,a or b)):
                        rows[name] += hit
                    all_hit &= a or b
                if all_hit:
                    accumulate(all_rows_hit, leaf)
    return {"fourteen_fours_whole_leaves": fourteen, "physical_leaf_totals": physical, "row_exclusions": rows,
            "prefix_certified_leaves": prefixes, "nonempty_leaves_whose_rows_all_fail_union": all_rows_hit}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ROOT/"data/revision/rank30_completed.zip")
    parser.add_argument("--verify", type=Path, help="require exact agreement with this JSON report")
    args = parser.parse_args()
    result = audit(args.archive)
    if args.verify and result != json.loads(args.verify.read_bytes()):
        raise ValueError("rank-30 mass projection disagrees")
    print(json.dumps(result, indent=2, sort_keys=True))
