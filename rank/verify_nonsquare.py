#!/usr/bin/env python3
"""Verify the complete nonsquare-dimension searches at ranks 29, 30 and 31."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "revision" / "nonsquare"
SOURCES = ROOT / "current" / "nonsquare"
sys.path.insert(0, str(ROOT / "code"))
from current_screen import screen
from rank_audit import validate

SOURCE_HASHES = {
    "complementary_packets.cpp": "1f6994f2a8894024357762023edc66ba2a3c7d4b7762f21d85a5a87294b93d2b",
    "parity.cpp": "39dfefc5da62aed2c038ede22a5735222620007b1be66e0713ef55905c156cac",
}
MANIFEST_HASHES = {
    29: "3362732617a9022b0062da0ed0331973ce217c7dc1422da89319b20948749347",
    30: "dda7d394ecd32323f655e9862ff914026ea868b0c663076a84b21afdd1f41504",
    31: "a147dc677d16cf19e11289f06c55db8d1e9317d060d9d844a73a7ccdfc34b5b9",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(blob):
    return sha256(blob).hexdigest()


def canonical(rows):
    return "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in sorted(rows)).encode("ascii")


def root_map(rank):
    """Map a possible largest-object codegree to its unique nonsquare root (q,b)."""
    roots = {}
    for q in range(2, rank):
        divisors = range(2, isqrt(q) + 1)
        if any(q % (i*i) == 0 for i in divisors):
            continue
        prime = all(q % i != 0 for i in divisors)
        for b in range(1, isqrt((rank-1)//q) + 1):
            f = q*b*b
            if f >= 8 and not (b == 1 and prime):
                roots[f] = (q, b)
    return roots


def replay_job(rank, receipt, expected, executables, build_dir):
    source = receipt["source"]
    if source not in executables:
        executable = Path(build_dir) / Path(source).stem
        subprocess.run(["c++", "-O3", "-std=c++17", str(SOURCES/source),
                        "-o", str(executable)], check=True)
        executables[source] = executable
    job = subprocess.run([str(executables[source]), str(rank), "8",
                          str(receipt["q"]), str(receipt["b"])],
                         check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    match = re.fullmatch(rb"Complete rank (\d+), nodes (\d+), types (\d+), seconds (\S+)",
                         job.stderr.splitlines()[-1])
    context = f"rank {rank}, q={receipt['q']}, b={receipt['b']}"
    require(match is not None, f"missing completion marker at {context}")
    require(tuple(map(int, match.groups()[:3])) == (rank, receipt["nodes"], receipt["types"]),
            f"regenerated node or type count differs at {context}")
    regenerated = sorted(json.loads(line) for line in job.stdout.splitlines())
    require(regenerated == expected, f"regenerated degree list differs at {context}")


def verify(rank, replay, executables, build_dir):
    blob = (DATA/f"rank{rank}.json").read_bytes()
    require(digest(blob) == MANIFEST_HASHES[rank], f"rank-{rank} manifest changed")
    manifest = json.loads(blob)
    require(manifest["rank"] == rank and manifest["scope"] == "q > 1" and manifest["complete"] is True,
            f"incorrect scope at rank {rank}")
    root_by_f = root_map(rank)
    jobs = manifest["roots"]
    require([(r["q"],r["b"]) for r in jobs] == sorted(root_by_f.values()),
            f"missing, duplicate or unexpected root at rank {rank}")
    require(all(r["complete"] is True for r in jobs), f"incomplete root at rank {rank}")
    sources = {r["source"] for r in jobs}
    require(manifest["sources"] == {name: SOURCE_HASHES[name] for name in sources},
            f"unexpected generator at rank {rank}")
    require(manifest["data"] == f"rank{rank}.jsonl", f"wrong data filename at rank {rank}")
    raw = (DATA/manifest["data"]).read_bytes()
    require(digest(raw) == manifest["sha256"], f"rank-{rank} data digest differs")
    rows = [json.loads(line) for line in raw.splitlines()]
    require(len(rows) == manifest["types"] and len(set(map(tuple,rows))) == len(rows)
            and canonical(rows) == raw, f"rank-{rank} data not unique and canonical")
    grouped = defaultdict(list)
    for row in rows:
        dimension = validate(row, rank)
        top = row[-1]
        f = dimension//(top*top)
        require(f in root_by_f, f"type has no admissible root at rank {rank}")
        q, b = root_by_f[f]
        M = b*top
        require(dimension == q*M*M and all(M % y == 0 for y in row),
                f"inconsistent common denominator at rank {rank}")
        grouped[q,b].append(row)
    for job in jobs:
        root_rows = grouped[job["q"],job["b"]]
        require(len(root_rows) == job["types"] and digest(canonical(root_rows)) == job["sha256"],
                f"rank-{rank} root digest or count differs for {(job['q'],job['b'])}")
        if replay:
            replay_job(rank, job, root_rows, executables, build_dir)
    report = screen(rows, rank)
    require(report["profiles"] == len(rows) and sum(report["first_exclusions"].values()) == len(rows)
            and not report["survivors"] and manifest["survivors"] == 0,
            f"current screen leaves an unresolved rank-{rank} type")
    print(f"Rank {rank}: {len(jobs)} nonsquare roots, {len(rows)} profiles, zero survivors.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", nargs="+", choices=sorted(MANIFEST_HASHES), type=int,
                        default=sorted(MANIFEST_HASHES))
    parser.add_argument("--replay", action="store_true", help="recompile and rerun every selected root")
    args = parser.parse_args()
    for name, expected in SOURCE_HASHES.items():
        require(digest((SOURCES/name).read_bytes()) == expected, f"generator digest differs for {name}")
    with tempfile.TemporaryDirectory() as build_dir:
        executables = {}
        for rank in args.ranks:
            verify(rank, args.replay, executables, build_dir)


if __name__ == "__main__":
    main()
