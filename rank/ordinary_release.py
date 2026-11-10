#!/usr/bin/env python3
"""Check or regenerate the completed ordinary simple-perfect modular search."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import next_rank as engine

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SMALL = DATA / "ordinary_release_6_27.zip"


def large_archive(rank):
    if rank == 30:
        return DATA / "revision/rank30_completed.zip"
    return DATA / "revision" / f"rank{rank}_optimized.zip"


def validate_rows(rows, rank):
    ordered = sorted(tuple(row) for row in rows)
    if len(ordered) != len(set(ordered)):
        raise ValueError(f"duplicate profile at rank {rank}")
    digest = sha256()
    for row in ordered:
        raw = engine.encoded(row)
        engine.degree_vector(raw, rank)
        dimension = sum(d*d for d in row)
        if any(dimension == 9*d*d for d in row):
            raise ValueError("ordinary codegree-nine root admitted")
        digest.update(raw + b"\n")
    return ordered, digest.hexdigest()


def small_generation(rank, executables, jobs):
    # The manifest kernel includes the no-nine theorem at the first denominator.
    keys = engine.uniform_manifest(executables["exact"], rank, 1)
    if any(q*b[0]**2 == 9 for q, b in keys):
        raise ValueError("codegree-nine root in ordinary manifest")
    def run(key):
        return engine.compute(executables, rank, key, use_theory=False)
    with ThreadPoolExecutor(jobs) as pool:
        receipts = list(pool.map(run, keys))
    return [row for receipt in receipts for row in receipt["rows"]]


def check_small(rank, executables, jobs, replay, gmp):
    with ZipFile(SMALL) as z:
        rows = [json.loads(line) for line in z.read(f"{rank}.jsonl").splitlines()]
    ordered, digest = validate_rows(rows, rank)
    if replay:
        regenerated, _ = validate_rows(small_generation(rank, executables, jobs), rank)
        if regenerated != ordered:
            raise ValueError(f"regenerated partition differs at rank {rank}")
    if gmp:
        proc = subprocess.run([executables["exact"], "--all", str(rank)],
                              capture_output=True, check=True)
        lines = proc.stdout.splitlines()
        marker = re.fullmatch(rf"COMPLETE search rank={rank} nodes=\d+ types={len(lines)} shards=0\n",
                              proc.stderr.decode("ascii"))
        exact, _ = validate_rows([json.loads(line) for line in lines], rank)
        if not marker or exact != ordered:
            raise ValueError(f"independent GMP-only traversal differs at rank {rank}")
    result = engine.screen(ordered, rank)
    return dict(result, degree_lists_sha256=digest)


def check_large(rank, executables, hashes, jobs, replay):
    archive = large_archive(rank)
    # Rank 30 has a second, independently enumerated nonsquare partition.
    if rank == 30:
        from verify_rank30 import summarize
        audit = summarize(archive, DATA / "revision")
        if audit != json.loads((DATA / "revision/rank30_audit.json").read_bytes()):
            raise ValueError("rank-30 checkpoint or separate nonsquare partition differs")
        if replay:
            subprocess.run([sys.executable, ROOT / "verify_nonsquare.py", "--ranks",
                            "30", "--replay"], check=True)
    engine.verify(archive, executables, hashes, replay, jobs)
    with ZipFile(archive) as z:
        rows = [row for line in z.open("receipts.jsonl") for row in json.loads(line)["rows"]]
        result = json.loads(z.read("screen.json"))
    _, digest = validate_rows(rows, rank)
    return dict(result, degree_lists_sha256=digest)


def main():
    if not __debug__:
        raise RuntimeError("Run the release verification with Python assertions enabled.")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ranks", type=int, nargs="+", choices=range(6,31),
                   default=list(range(6,31)), metavar="R")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--replay", action="store_true",
                   help="regenerate the actual candidate lists and compare every row")
    p.add_argument("--gmp", action="store_true",
                   help="also compare a direct all-GMP traversal, for ranks 6–27")
    args = p.parse_args()
    if args.gmp and max(args.ranks) > 27:
        p.error("--gmp applies to ranks 6–27 only")
    if args.jobs < 1:
        p.error("at least one worker is required")
    ledger = json.loads((DATA / "ordinary_release.json").read_bytes())
    with TemporaryDirectory() as build:
        executables, hashes = engine.executables(Path(build))
        if hashes != ledger["source_sha256"]:
            raise ValueError("release kernel mismatch")
        for rank in dict.fromkeys(args.ranks):
            if rank <= 27:
                got = check_small(rank, executables, args.jobs, args.replay, args.gmp)
            else:
                got = check_large(rank, executables, hashes, args.jobs, args.replay)
            if got != ledger["ranks"][str(rank)] or got["survivors"]:
                raise ValueError(f"incomplete exclusions or ledger mismatch at rank {rank}")
            mode = "regenerated" if args.replay or args.gmp else "verified from archive"
            print(f"Rank {rank}: {got['profiles']} lists {mode}, zero survivors.", flush=True)


if __name__ == "__main__":
    main()
