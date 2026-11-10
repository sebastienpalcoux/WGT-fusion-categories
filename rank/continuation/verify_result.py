#!/usr/bin/env python3
"""Verify a returned rank search, regenerate its manifest and independently repeat the Python screening."""
import argparse
import hashlib
import io
import json
from collections import Counter
from math import isqrt
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from screen import run as screen
from rank_audit import fast_report
from orbit_screen import screen_survivors


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def degree_vector(line, rank):
    d = json.loads(line)
    require(isinstance(d, list) and len(d) == rank and all(type(x) is int and x > 0 for x in d),
            "malformed dimension vector")
    require(d == sorted(d) and d.count(1) == 1, "invalid type order or multiplicity of the unit")
    total = sum(x*x for x in d)
    require(all(total % (x*x) == 0 for x in d), "square divisibility failed")
    return tuple(d)


def verify_one(archive, compiled_source, replay, roots=None):
    members = archive.namelist()
    require(len(members) == len(set(members)), "duplicated certificate member")
    search = json.loads(archive.read("search.json"))
    source_digest = search["source_sha256"]
    exe, historical = compiled_source(source_digest)
    require(search.get('roots') == roots, "certificate root scope disagrees")
    complete = json.loads(archive.read("complete.json"))
    rank, depth = search["rank"], search["depth"]
    require(type(rank) is int and 6 <= rank <= 48 and type(depth) is int and 1 <= depth <= rank - 3,
            "invalid rank or split depth")
    require(source_digest == search["source_sha256"] and all(complete[k] == v for k, v in search.items()),
            "source or configuration mismatch")
    generated = subprocess.run([exe, "--manifest", str(rank), str(depth)], capture_output=True, check=True)
    marker = re.fullmatch(rb"COMPLETE manifest rank=(\d+) nodes=\d+ types=0 shards=(\d+)\n", generated.stderr)
    require(marker and int(marker[1]) == rank and int(marker[2]) == len(generated.stdout.splitlines()), "bad manifest completion")
    expected = generated.stdout
    if roots is not None:
        selected = set(map(tuple, roots))
        def selected_root(raw):
            q,prefix = raw.split(b"\t")
            return (int(q),int(prefix.split(b",")[0])) in selected
        expected = b"".join(raw for raw in expected.splitlines(keepends=True) if selected_root(raw))
    manifest = archive.read("manifest.tsv")
    require(expected == manifest and sha(manifest) == complete["manifest_sha256"], "manifest mismatch")
    coordinates = []
    for line in manifest.splitlines():
        q, prefix = line.split(b"\t")
        coordinates.append((int(q), tuple(map(int, prefix.split(b",")))))
    require(len(coordinates) == complete["shards"] and len(set(coordinates)) == len(coordinates),
            "duplicate or missing manifest coordinates")
    raw = archive.read("candidates.jsonl")
    require(sha(raw) == complete["candidate_sha256"], "candidate digest mismatch")
    ordered = [degree_vector(line, rank) for line in raw.splitlines()]
    types = set(ordered)
    require(len(types) == len(ordered) == complete["candidate_count"] and ordered == sorted(ordered),
            "duplicate, unsorted or missing candidate")
    if "family_result.json" in members:
        final = fast_report([list(t) for t in ordered], rank)
        require(final == json.loads(archive.read("family_result.json")), "family and block audit disagrees")
        final_survivors = final["survivors"]
    else:
        trace = json.loads(archive.read("trace_result.json"))
        calculated = screen([list(t) for t in ordered])
        require(calculated == trace, "trace audit disagrees with reported screening")

        if not historical:
            require("orbit_result.json" in members, "modern result is missing its orbit screening")
        final_survivors = trace["survivors"]
        if "orbit_result.json" in members:
            orbit = json.loads(archive.read("orbit_result.json"))
            require(screen_survivors(final_survivors, principal=True) == orbit, "orbit audit disagrees with result")
            final_survivors = orbit["survivors"]
    if "shard_receipts.jsonl" in members:
        hashes = {}
        counts = Counter()
        seen = set()
        with archive.open("shard_candidates.jsonl") as lines:
            for line in lines:
                label, payload = line.split(b"\t", 1)
                idx = int(label)
                require(0 <= idx < len(coordinates), "invalid shard index")
                q, prefix = coordinates[idx]
                t = degree_vector(payload, rank)
                require(t in types and t not in seen, "unexpected shard type")
                seen.add(t)
                D = sum(x*x for x in t)
                M = isqrt(D//q)
                require(q*M*M == D and all(M%x == 0 for x in t) and
                        tuple(M//x for x in reversed(t[-depth:])) == prefix, "type assigned to wrong shard")
                if idx not in hashes:
                    hashes[idx] = hashlib.sha256()
                hashes[idx].update(payload)
                counts[idx] += 1
        require(seen == types, "missing types in shard streams")
        empty_digest = sha(b"")
        total_nodes = 0
        n_receipts = 0
        with archive.open("shard_receipts.jsonl") as lines:
            for line in lines:
                idx, nodes, n, digest = json.loads(line)
                require(idx == n_receipts and n == counts[idx] and type(nodes) is int and nodes >= 0,
                        "invalid shard receipt")
                require(digest == (hashes[idx].hexdigest() if idx in hashes else empty_digest), "shard digest mismatch")
                total_nodes += nodes
                n_receipts += 1
        require(n_receipts == complete["shards"] and total_nodes == complete["nodes_summed_over_shards"],
                "shard receipt totals disagree")
        receipt_note = "per-shard receipts verified"
    else:
        require(historical and "shard_candidates.jsonl" not in members, "missing shard receipts")
        receipt_note = "legacy archive: aggregate receipt only"
    if replay and roots is not None:
        rerun = []
        for q,b in roots:
            result = subprocess.run([exe, "--shard", str(rank), str(q), str(b)], capture_output=True, check=True)
            mark = re.fullmatch(rb"COMPLETE search rank=(\d+) nodes=\d+ types=(\d+) shards=0\n", result.stderr)
            part = [degree_vector(line, rank) for line in result.stdout.splitlines()]
            require(mark and int(mark[1]) == rank and int(mark[2]) == len(part), "incomplete root replay")
            rerun.extend(part)
        require(len(rerun) == len(types) and set(rerun) == types, "GMP root replay disagrees")
        receipt_note += ", both GMP roots independently regenerated"
    elif replay:
        proc = subprocess.run([exe, "--all", str(rank)], capture_output=True, check=True)
        marker = re.fullmatch(rb"COMPLETE search rank=(\d+) nodes=\d+ types=(\d+) shards=0\n", proc.stderr)
        rerun = [degree_vector(line, rank) for line in proc.stdout.splitlines()]
        require(marker and int(marker[1]) == rank and int(marker[2]) == len(rerun) == len(ordered)
                and set(rerun) == types, "monolithic replay disagrees")
        receipt_note += ", all profiles independently regenerated without shards"
    print(f"Rank {rank}: {len(ordered)} candidates, {len(final_survivors)} survive; {receipt_note}.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path, help="one rank*_result.zip or a combined ranks_*_results.zip")
    parser.add_argument("--replay", action="store_true", help="also regenerate every candidate in one process per rank")
    args = parser.parse_args()
    names = ("search.cpp", "rank24_search.cpp", "rank20_search.cpp")
    sources = {}
    for name in names:
        source = Path(__file__).resolve().with_name(name)
        if source.is_file():
            sources[sha(source.read_bytes())] = source
    with tempfile.TemporaryDirectory() as tmp, ZipFile(args.result) as archive:
        def compiled_source(digest):
            require(digest in sources, "no bundled search source matches the certificate SHA-256")
            exe = Path(tmp) / digest
            if not exe.exists():
                subprocess.run(["c++", "-O3", "-std=c++17", str(sources[digest]),
                                "-lgmpxx", "-lgmp", "-o", str(exe)], check=True)
            return exe, sources[digest].name == "rank20_search.cpp"

        if "search.json" in archive.namelist():
            verify_one(archive, compiled_source, args.replay)
        else:
            entries = archive.namelist()
            require(entries and len(entries) == len(set(entries))
                    and all(re.fullmatch(r"rank\d+_result\.zip", name) for name in entries),
                    "unexpected combined archive layout")
            for name in entries:
                with ZipFile(io.BytesIO(archive.read(name))) as nested:
                    require(name == f'rank{json.loads(nested.read("search.json"))["rank"]}_result.zip',
                            "nested filename does not match its rank")
                    verify_one(nested, compiled_source, args.replay)


if __name__ == "__main__":
    main()
