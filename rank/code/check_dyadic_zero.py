#!/usr/bin/env python3
"""Count the strengthened dyadic lcm obstruction in the original rank-30 receipts."""
import argparse
import json
from collections import Counter
from math import isqrt, lcm
from pathlib import Path
from zipfile import ZipFile


def default_archive():
    roots = [Path(__file__).resolve().parent, *Path(__file__).resolve().parents, Path.cwd()]
    names = ("rank30_completed.zip", "data/ordinary/rank30_completed.zip",
             "data/revision/rank30_completed.zip", "upload/rank30_completed.zip")
    for root in roots:
        for name in names:
            path = root / name
            if path.is_file():
                return path
    raise FileNotFoundError("Pass the path of rank30_completed.zip")


def count(archive):
    stats = Counter()
    with ZipFile(archive) as bundle, bundle.open("receipts.jsonl") as receipts:
        for line in receipts:
            shard = json.loads(line)
            if shard["q"] != 1 or shard["prefix"][0] != 4:
                continue
            for degrees in shard["rows"]:
                dimension = sum(d * d for d in degrees)
                root = isqrt(dimension)
                if root * root != dimension or degrees.count(1) != 1:
                    raise ValueError("Non-square or imperfect archived input")
                denominators = [root // d for d in degrees if d != 1]
                if any(root % d for d in degrees):
                    raise ValueError("Invalid archived degree")
                if root % 8:
                    continue
                stats["eligible"] += 1
                if lcm(*(b for b in denominators if b % 2 == 0)) != root:
                    continue
                stats["old_lcm_pass"] += 1
                if lcm(*(b for b in denominators if b % 4 == 0)) != root:
                    stats["additional_lcm_exclusions"] += 1
    return dict(stats)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path)
    args = parser.parse_args()
    results = count(args.archive or default_archive())
    print(json.dumps(results, indent=2, sort_keys=True))
    assert results == {"eligible": 92943, "old_lcm_pass": 80876,
                       "additional_lcm_exclusions": 29873}
