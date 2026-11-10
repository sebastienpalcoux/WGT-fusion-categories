#!/usr/bin/env python3
"""Check the rank-21/22/23 candidate lists, their exclusions and optional complete native regenerations."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from time import monotonic

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
from orbit_screen import first_obstruction, small_codegree_obstruction
from screen import run as screen


def require(test, message):
    if not test:
        raise ValueError(message)


def verify(rank, executable=None, timeout=None):
    start = monotonic()
    base = ROOT / "data/ordinary"
    raw = (base / f"rank{rank}.jsonl").read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    require(all(isinstance(row, list) and len(row) == rank and
                all(type(d) is int and d > 0 for d in row) for row in rows), "malformed degree vector")
    require(rows == sorted(rows) and len({tuple(row) for row in rows}) == len(rows), "duplicate or unordered list")
    primary = screen(rows)
    excluded = Counter()
    remaining = []
    for row in primary["survivors"]:
        if rank >= 23:
            f = small_codegree_obstruction(row)
            if f is not None:
                excluded[str(f)] += 1
                continue
        failure = first_obstruction(row, principal=(rank >= 22))
        if failure is None:
            remaining.append(row)
        else:
            excluded[str(failure[1])] += 1
    calculated = {"rank":rank, "candidates":len(rows), "primary":primary["counts"],
                  "final_codegrees":dict(excluded), "survivors":remaining}
    expected = json.loads((base / f"rank{rank}_summary.json").read_text())
    require(all(expected[key] == value for key, value in calculated.items()),
            f"rank-{rank} report disagrees with the candidate lists or the screening")
    require(not remaining, f"unexpected survivor at rank {rank}")
    if executable:
        command = [str(executable), str(rank), "8"]
        run = subprocess.run(command, check=True, capture_output=True, timeout=timeout)
        lines = run.stderr.decode().splitlines()
        marker = re.fullmatch(r"Complete rank (\d+), nodes (\d+), types (\d+), seconds [0-9.eE+\-]+",
                             lines[-1] if lines else "")
        rerun = sorted(json.loads(line) for line in run.stdout.splitlines())
        require(marker and int(marker[1]) == rank and int(marker[2]) > 0 and int(marker[3]) == len(rows)
                and rerun == rows, f"rank-{rank} regeneration incomplete or different")
    mode = "profiles regenerated and screened" if executable else "stored profiles screened"
    print(f"Rank {rank}: {len(rows)} {mode}, {len(remaining)} unresolved; {monotonic()-start:.2f}s.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, nargs="+", choices=(21, 22, 23), default=[21, 22, 23])
    parser.add_argument("--replay", action="store_true", help="also regenerate the complete raw candidate list")
    parser.add_argument("--timeout", type=float, help="optional seconds allowed for each native replay")
    args = parser.parse_args()
    if args.timeout is not None and (not args.replay or args.timeout <= 0):
        parser.error("--timeout requires --replay and a positive number of seconds")
    with tempfile.TemporaryDirectory() as tmp:
        executable = None
        if args.replay:
            executable = Path(tmp) / "ordinary_native"
            subprocess.run(["c++", "-O3", "-std=c++17", str(ROOT / "code/ordinary_native.cpp"), "-o",
                            str(executable)], check=True)
        for rank in dict.fromkeys(args.ranks):
            verify(rank, executable, args.timeout)


if __name__ == "__main__":
    main()
