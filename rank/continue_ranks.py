#!/usr/bin/env python3
"""Run consecutive ordinary modular ranks, resume completed shards, and package their certificates."""
import argparse
import fcntl
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
from rank_audit import fast_report, read_rows

BASE_FILES = ("manifest.tsv", "search.json", "complete.json", "candidates.jsonl", "family_result.json")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def package(rank, work):
    """Include the evidence for each shard, including a record for each empty output."""
    complete = json.loads((work / "complete.json").read_text())
    if complete["rank"] != rank or complete["source_sha256"] != sha256(ROOT / "continuation/search.cpp"):
        raise ValueError("completion receipt disagrees with rank or search source")
    for key, name in (("candidate_sha256", "candidates.jsonl"), ("manifest_sha256", "manifest.tsv")):
        if sha256(work / name) != complete[key]:
            raise ValueError(f"{name} disagrees with completion receipt")
    result = work / f"rank{rank}_result.zip"
    temporary = result.with_suffix(".zip.partial")
    with ZipFile(temporary, "w", compression=ZIP_DEFLATED) as archive:
        for name in BASE_FILES:
            archive.write(work / name, name)
        # Each row is [manifest index, node count, count, original output SHA-256].
        # Whole-cluster archives need only a single file for the empty shard receipts.
        with archive.open("shard_receipts.jsonl", "w") as receipts:
            for index in range(complete["shards"]):
                stem = work / "parts" / f"{index:08d}"
                receipt = json.loads(stem.with_suffix(".receipt").read_text())
                fields = dict(field.split("=") for field in receipt["completion"].split()[2:])
                row = [index, int(fields["nodes"]), int(fields["types"]), receipt["sha256"]]
                receipts.write((json.dumps(row, separators=(",", ":")) + "\n").encode())
        # Preserve original emission order so all returned outputs can be hashed independently.
        with archive.open("shard_candidates.jsonl", "w") as outputs:
            for index in range(complete["shards"]):
                path = work / "parts" / f"{index:08d}.jsonl"
                with path.open("rb") as file:
                    for row in file:
                        outputs.write(str(index).encode() + b"\t" + row)
    temporary.replace(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, nargs="+", default=[28], metavar="R",
                        help="exact ranks to run in order (default: 28)")
    parser.add_argument("--verify-only", action="store_true",
                        help="check existing shards and repackage; fail if any shard is missing or invalid")
    parser.add_argument("--jobs", type=int, default=12, help="maximum simultaneous searches (default: 12)")
    parser.add_argument("--depth", type=int, help="split after this many nonunit denominators (default: min(12, rank-3))")
    parser.add_argument("--work", type=Path, default=Path("rank-work"), help="working directory (default: rank-work)")
    args = parser.parse_args()
    if args.jobs < 1 or len(set(args.ranks)) != len(args.ranks) or any(r < 6 or r > 48 for r in args.ranks):
        parser.error("jobs must be positive and the ranks must be distinct integers from 6 through 48")
    if args.depth is not None and not 1 <= args.depth <= min(args.ranks) - 3:
        parser.error("depth must lie between 1 and the smallest requested rank minus 3")
    base = args.work.resolve()
    base.mkdir(parents=True, exist_ok=True)
    with (base / ".series.lock").open("a") as series_lock:
        try:
            fcntl.flock(series_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.exit(2, f"Another coordinator is using {base}.\n")
        archives = []
        for rank in args.ranks:
            work = base / f"rank{rank}"
            depth = args.depth if args.depth is not None else min(12, rank - 3)
            command = [sys.executable, str(ROOT / "continuation/run_shards.py"), str(rank),
                       "--jobs", str(args.jobs), "--depth", str(depth), "--work", str(work)]
            if args.verify_only:
                command.append("--verify-only")
            print(f"Starting or resuming rank {rank}, split depth {depth}.", flush=True)
            subprocess.run(command, check=True)
            # The coordinator has already checked every receipt and output before merging.
            with (work / ".lock").open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                rows = read_rows(work / "candidates.jsonl", rank)
                final = fast_report(rows, rank)
                (work / "family_result.json").write_text(json.dumps(final, indent=2) + "\n")
                print(f"After the family and block tests, {len(final['survivors'])} profiles remain.", flush=True)
                archives.append(package(rank, work))
            print(f"Rank {rank} finished. Returnable result: {archives[-1]}", flush=True)
        if len(archives) == 1:
            return
        combined = base / ("ranks_" + "_".join(map(str, args.ranks)) + "_results.zip")
        temporary = combined.with_suffix(".zip.partial")
        with ZipFile(temporary, "w") as archive:
            for result in archives:
                archive.write(result, result.name)
        temporary.replace(combined)
        print(f"Return this combined file: {combined}", flush=True)


if __name__ == "__main__":
    main()
