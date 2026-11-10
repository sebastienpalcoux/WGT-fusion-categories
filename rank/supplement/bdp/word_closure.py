"""Exact saturation of the exposed-adjoint intersection calculus on F_s.

This checks the abstract word obstruction at a target prime. The underlying
common-constituent and exposure implications are proved in the manuscript.
An open result is a closed state of these inference rules only.
"""

import argparse
import json
from functools import lru_cache
from itertools import combinations
from math import comb, isqrt
from pathlib import Path


def least_overlap(h):
    j = (isqrt(1 + 4*h) - 1) // 2
    return j + (j*(j+1) < h)


def reduction(s, k, ell):
    if s < 3 or s % 2 == 0 or any(s % d == 0 for d in range(3, isqrt(s) + 1, 2)):
        raise ValueError("target must be an odd prime")
    if not 2 <= k <= s or ell not in (2, 3):
        raise ValueError("require 2 <= k <= s and ell in {2,3}")
    H = (s-1)//2
    E = comb(k, 2)
    if ell * (1 + 2*E) >= (ell-1)*s:
        return {"status": "initial exposure undecided"}
    h0 = next((h for h in range(1, E+1)
               if (H-1)*least_overlap(h) <= h*(h-1)), None)
    if h0 is None:
        return {"status": "excluded by scalar inequality"}
    b = comb(k-1, 2)
    upper = max(E**2, s*b*b + (E-b)**2)
    lower = s*H*least_overlap(h0)
    status = "one coset" if upper < lower else "higher conductor undecided"
    return {"status": status, "minimum_differences": h0,
            "scaled_overlap_upper": upper, "scaled_overlap_lower": lower}


class Field:
    def __init__(self, s, ell, h0):
        self.s, self.ell, self.h0 = s, ell, h0
        self.H = (s-1)//2
        self.mask = (1 << self.H)-1
        for gen in range(1, s):
            xs = [pow(gen, j, s) for j in range(self.H)]
            if len({min(x, s-x) for x in xs}) == self.H:
                break
        self.xs = xs
        self.bits = [0]*s
        for i,x in enumerate(xs):
            self.bits[x] = self.bits[-x] = 1 << i

    def rotate(self, mask, j):
        return ((mask << j) & self.mask) | (mask >> (self.H-j))

    @lru_cache(maxsize=65536)
    def canonical(self, mask):
        return min(self.rotate(mask, j) for j in range(self.H))

    def support_differences(self, a):
        out = 0
        for x,y in combinations(a, 2):
            out |= self.bits[(x-y) % self.s]
        return out

    @lru_cache(maxsize=65536)
    def difference(self, d):
        xs = [x for i,x in enumerate(self.xs) if d & (1 << i)]
        out = d
        for x in xs:
            out |= self.bits[2*x % self.s]
        for x,y in combinations(xs, 2):
            out |= self.bits[(x-y) % self.s] | self.bits[(x+y) % self.s]
        return out

    def exposed(self, mask):
        return self.ell * (1 + 2*mask.bit_count()) < (self.ell-1)*self.s

    def certificate(self, seed):
        if seed.bit_count() < self.h0:
            return {"steps": [], "end": "seed too small"}
        known = [seed]
        seen = {seed}
        steps = []
        # Appending a new mask gives a finite fair queue. When its turn
        # comes, it is paired with every earlier mask and with itself.
        for i,left in enumerate(known):
            for j,right in enumerate(known[:i+1]):
                for u in range(self.H):
                    common = left & self.rotate(right, u)
                    nxt = self.difference(common)
                    if nxt.bit_count() < self.h0:
                        return {"steps": steps, "end": [i, j, u]}
                    if self.exposed(nxt):
                        nxt = self.canonical(nxt)
                        if nxt not in seen:
                            seen.add(nxt)
                            steps.append([i, j, u, nxt])
                            known.append(nxt)
        return {"steps": steps, "end": "closed", "masks": known}


def affine_class(a, s):
    """Lexicographically least representative containing both 0 and 1."""
    return min(tuple(sorted((v-x)*pow(y-x, -1, s) % s for v in a))
               for x in a for y in a if x != y)


def seeds(field, k):
    counts = {}
    for tail in combinations(range(2, field.s), k-2):
        a = (0, 1, *tail)
        seed = field.canonical(field.support_differences(a))
        counts[seed] = counts.get(seed, 0) + 1
    return counts


def solve(s, k, ell):
    red = reduction(s, k, ell)
    result = {"prime": s, "positions": k, "ell": ell, "reduction": red}
    if red["status"] != "one coset":
        return result
    field = Field(s, ell, red["minimum_differences"])
    counts = seeds(field, k)
    records = {str(seed): field.certificate(seed) for seed in sorted(counts)}
    unresolved = {seed for seed in counts if records[str(seed)]["end"] == "closed"}
    result.update(records=records, normalized_supports=sum(counts.values()),
                  unresolved_supports=sum(counts[x] for x in unresolved))
    return result


def checked_record(field, seed, record):
    """Independent direct set arithmetic for replay, including closure checks."""
    s, H = field.s, field.H

    @lru_cache(maxsize=32768)
    def residue_set(mask):
        xs = [x for i,x in enumerate(field.xs) if mask & (1 << i)]
        return frozenset([0, *xs, *(s-x for x in xs)])

    @lru_cache(maxsize=32768)
    def consequence(left, right, u):
        rotated = {field.xs[u]*x % s for x in residue_set(right)}
        a = residue_set(left) & rotated
        diffs = {(x-y) % s for x in a for y in a}
        return sum({field.bits[x] for x in diffs})

    def parents(indices, known):
        if len(indices) != 3 or any(type(i) is not int for i in indices):
            raise ValueError("invalid inference")
        i,j,u = indices
        if not (0 <= i < len(known) and 0 <= j < len(known) and 0 <= u < H):
            raise ValueError("invalid premise or multiplier")
        return consequence(known[i], known[j], u)

    if seed.bit_count() < field.h0:
        return record == {"steps": [], "end": "seed too small"}
    if not field.exposed(seed):
        raise ValueError("initial difference bound is unexposed")
    known = [seed]
    seen = {seed}
    for step in record["steps"]:
        if len(step) != 4:
            return False
        nxt = parents(step[:3], known)
        mask = field.canonical(nxt)
        if (type(step[3]) is not int or step[3] != mask or mask in seen
                or nxt.bit_count() < field.h0 or not field.exposed(nxt)):
            return False
        seen.add(mask)
        known.append(mask)
    end = record["end"]
    if isinstance(end, list):
        return parents(end, known).bit_count() < field.h0
    if end != "closed" or record["masks"] != known:
        return False
    for i,left in enumerate(known):
        for right in known[:i+1]:
            for u in range(H):
                nxt = consequence(left, right, u)
                if nxt.bit_count() < field.h0:
                    return False
                if field.exposed(nxt) and field.canonical(nxt) not in seen:
                    return False
    return True


def verify(result):
    s, k, ell = (result[x] for x in ("prime", "positions", "ell"))
    red = reduction(s, k, ell)
    if red != result["reduction"]:
        return False
    if red["status"] != "one coset":
        return set(result) == {"prime", "positions", "ell", "reduction"}
    field = Field(s, ell, red["minimum_differences"])
    counts = seeds(field, k)
    records = result["records"]
    if set(records) != {str(x) for x in counts}:
        return False
    if not all(checked_record(field, x, records[str(x)]) for x in counts):
        return False
    unresolved = sum(n for x,n in counts.items() if records[str(x)]["end"] == "closed")
    return (result["normalized_supports"] == comb(s-2, k-2)
            and result["unresolved_supports"] == unresolved)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check")
    check.add_argument("s", type=int)
    check.add_argument("k", type=int)
    check.add_argument("--ell", type=int, choices=(2,3), default=3)
    check.add_argument("--write", type=Path, help="save the full replayable JSON certificate")
    check.add_argument("--classes", action="store_true", help="enumerate open affine classes")
    replay = sub.add_parser("verify")
    replay.add_argument("file", type=Path)
    args = parser.parse_args()
    if args.command == "verify":
        try:
            valid = verify(json.loads(args.file.read_text()))
        except (ValueError, TypeError, KeyError, IndexError):
            valid = False
        print("Certificate verified." if valid else "Invalid certificate.")
        return 0 if valid else 1
    result = solve(args.s, args.k, args.ell)
    if args.write:
        args.write.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    summary = {k: v for k,v in result.items() if k != "records"}
    if "records" in result:
        summary["initial_mask_classes"] = len(result["records"])
        if args.classes:
            field = Field(args.s, args.ell, result["reduction"]["minimum_differences"])
            classes = set()
            for tail in combinations(range(2,args.s), args.k-2):
                a = (0,1,*tail)
                d = field.canonical(field.support_differences(a))
                if result["records"][str(d)]["end"] == "closed":
                    classes.add(affine_class(a, args.s))
            summary["open_affine_classes"] = sorted(classes)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
