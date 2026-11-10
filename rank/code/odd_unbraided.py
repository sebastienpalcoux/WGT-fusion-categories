"""Exact paired formal-codegree certificate for the arbitrary odd-dimensional case.

Place beside catalogue.py. Run python3 code/odd_unbraided.py --verify data/odd_unbraided.json.
A standalone development run can point --catalogue at the directory containing catalogue.py.
"""
import argparse
import json
import sys
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm, prod
from pathlib import Path


def enumerate_pairs(rank, factors):
    """Emit sorted f_1,...,f_s,D with 2 sum(1/f_i)+1/D=1 and f_i | D."""
    def visit(num, den, floor, multiple, prefix, left):
        if left == 1:
            divisors = [1]
            for prime, exponent in factors(den):
                divisors += [u * prime ** j for u in divisors for j in range(1, exponent + 1)]
            for k in sorted(divisors):
                if (2 * k + 1) % num:
                    continue
                dimension = den * ((2 * k + 1) // num)
                last = dimension // k
                if dimension % 2 and dimension % multiple == 0 and last >= floor:
                    yield (*prefix, last, dimension)
            return
        lower = max(floor, (2 * den) // num + 1) | 1
        upper = (2 * left + 1) * den // num
        for current in range(lower, upper + 1, 2):
            a, b = num * current - 2 * den, den * current
            common = gcd(a, b)
            yield from visit(a // common, b // common, current, lcm(multiple, current),
                             (*prefix, current), left - 1)
    yield from visit(1, 1, 3, 1, (), (rank - 1) // 2)


def two_squares(n):
    if n <= 1:
        return False
    return any(isqrt(n - x * x) ** 2 == n - x * x for x in range(1, isqrt(n - 1) + 1))


def repeated_square(n):
    return any(isqrt(n - 2 * x * x) ** 2 == n - 2 * x * x
               for x in range(1, isqrt((n - 1) // 2) + 1))


def classify(row, rank, covered):
    *fs, dimension = row
    if covered(dimension):
        return "center", None
    if rank != 11:
        return "open", None
    if fs[0] == 3:
        if dimension % 9 != 3:
            return "stabilizer-3", None
        u = [(dimension // f - 1) % 3 for f in fs]
        if sum(u) != 8:
            return "cubic-induction-norm", None
        vacuum_gap = dimension // 3 - sum(x * (dimension // f) for x, f in zip(u, fs))
        if vacuum_gap == 1 and dimension not in fs and covered(dimension // 3):
            return "central-C3", None
    if len(set(fs)) == len(fs) and isqrt(prod(row)) ** 2 != prod(row):
        return "real-determinant", None
    if fs[:3] == [5, 5, 11] and fs.count(11) == 1:
        target, remainder = divmod(dimension - 55, 6050)
        if remainder or not two_squares(target):
            return "dimension-55", row
    if fs[:2] == [3, 7] and fs.count(7) == 1:
        target, remainder = divmod(dimension - 21, 882)
        if remainder:
            return "dimension-21", row
        if (fs == [3, 7, 55, 179, 24323] and dimension == 5028658635
                and vacuum_gap == 1 and dimension not in fs
                and not repeated_square(target) and dimension % 147 == 21
                and covered(dimension // 21)):
            return "two-central-descents", row
    return "open", row


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=Path(__file__).resolve().parent)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--verify", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.catalogue.resolve()))
    from catalogue import certificate
    from families import factors

    @lru_cache(maxsize=None)
    def covered(dimension):
        return dimension < 5400 or bool(certificate(dimension * dimension))

    report = {}
    for rank in (3, 5, 7, 9, 11):
        counts = Counter()
        last_cases = []
        for row in enumerate_pairs(rank, factors):
            rule, residual = classify(row, rank, covered)
            counts[rule] += 1
            if residual is not None:
                last_cases.append({"row": list(row), "rule": rule})
        if counts.get("open"):
            raise ValueError("incomplete at rank " + str(rank))
        report[str(rank)] = {"total": sum(counts.values()), "rules": dict(sorted(counts.items())),
                             "last_cases": sorted(last_cases, key=lambda item: item["row"])}
    if args.verify and report != json.loads(args.verify.read_text()):
        raise ValueError("receipt mismatch")
    contents = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.write:
        args.write.write_text(contents)
    print(contents, end="")


if __name__ == "__main__":
    run()
