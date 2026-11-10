#!/usr/bin/env python3
"""Verify the completed rank-30 return, its checkpoint, and its separate nonsquare search."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def coordinate(item):
    return item["q"], tuple(item["prefix"])


def summarize(archive, data):
    """Cross-check original and completed returns, then count physical leaf jobs once."""
    with ZipFile(archive) as z:
        spec = json.loads(z.read("spec.json"))
        require((spec["rank"], spec["depth"], spec["format"]) == (30, 11, 1),
                "unexpected rank-30 configuration")
        parents = [json.loads(line) for line in z.open("receipts.jsonl")]
        current = json.loads(z.read("screen.json"))
        payload_digests = {name: sha(z.read(name)) for name in sorted(z.namelist())}

    completed = {coordinate(item): item for item in parents}
    require(len(completed) == len(parents), "duplicate parent in completed return")
    checkpoint = data / "rank30_checkpoint.zip"
    with ZipFile(checkpoint) as z:
        require(json.loads(z.read("spec.json")) == spec, "checkpoint generator disagreement")
        original = [json.loads(line) for line in z.open("receipts.jsonl")]
        stored = json.loads(z.read("screen.json"))
        pending = {(q, tuple(prefix)) for q, prefix in json.loads(z.read("pending.json"))}
    old_keys = {coordinate(item) for item in original}
    require(len(old_keys) == len(original) and pending == set(completed) - old_keys,
            "checkpoint pending partition disagrees")
    require(all(completed.get(coordinate(item)) == item for item in original),
            "an original parent receipt changed")
    require(all("children" in completed[key] for key in pending),
            "unexpected representation of a formerly pending parent")

    leaves = [leaf for parent in parents for leaf in parent.get("children", [parent])]
    square_classes, largest_codegrees, rules = Counter(), Counter(), Counter()
    nonsquare = set()
    distinct_rows = set()
    for parent in parents:
        q, prefix = coordinate(parent)
        for row in parent["rows"]:
            key = tuple(row)
            require(key not in distinct_rows, "duplicate rank-30 degree sequence")
            distinct_rows.add(key)
            square_classes[q] += 1
            largest_codegrees[q * prefix[0] ** 2] += 1
            if q > 1:
                nonsquare.add(key)
    separate = [tuple(json.loads(line)) for line in (data / "nonsquare/rank30.jsonl").read_bytes().splitlines()]
    require(len(separate) == len(nonsquare) and set(separate) == nonsquare,
            "separate nonsquare enumeration produced a different set")
    for leaf in leaves:
        if "rule" in leaf:
            rules[leaf["rule"]] += 1

    recorded = current["first_exclusions"]
    previous = stored["first_exclusions"]
    require(not stored["survivors"] and not current["survivors"], "unresolved degree sequence")
    delta = {name: recorded.get(name, 0) - previous.get(name, 0) for name in recorded.keys() | previous.keys()}
    added = sum(len(completed[key]["rows"]) for key in pending)
    require(len(distinct_rows) == current["profiles"] == sum(recorded.values()), "completed screen total")
    require(sum(len(item["rows"]) for item in original) == stored["profiles"] == sum(previous.values()),
            "checkpoint screen total")
    require(all(value >= 0 for value in delta.values()) and sum(delta.values()) == added,
            "incremental first exclusions disagree")
    result = {
        "rank": 30, "completed_parents": len(parents), "checkpoint_parents": len(original),
        "formerly_pending_parents": len(pending), "new_child_jobs": sum(len(completed[key]["children"]) for key in pending),
        "leaves": len(leaves), "searched_leaves": len(leaves) - sum(rules.values()),
        "theorem_leaves": dict(sorted(rules.items())), "nonempty_leaves": sum(bool(item["rows"]) for item in leaves),
        "native_nodes": sum(item["native_nodes"] for item in leaves),
        "gmp_nodes": sum(item["gmp_nodes"] for item in leaves), "gmp_handoffs": sum(item["fallbacks"] for item in leaves),
        "profiles": len(distinct_rows), "checkpoint_profiles": stored["profiles"], "additional_profiles": added,
        "by_square_class": {str(k): v for k, v in sorted(square_classes.items())},
        "by_largest_degree_codegree": {str(k): v for k, v in sorted(largest_codegrees.items())},
        "nonsquare_independent_matches": len(nonsquare),
        "first_exclusions": recorded, "additional_first_exclusions": {k: v for k, v in sorted(delta.items()) if v},
        "survivors": 0, "member_sha256": payload_digests,
        "checkpoint_sha256": sha(checkpoint.read_bytes()), "submitted_zip_sha256": sha(archive.read_bytes()),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent,
                        help="distribution directory containing next_rank.py")
    parser.add_argument("--archive", type=Path, help="completed ZIP (default data/revision/rank30_completed.zip)")
    parser.add_argument("--report", type=Path, help="expected JSON (default data/revision/rank30_audit.json)")
    parser.add_argument("--replay", action="store_true",
                        help="also rerun the 1,938 searched leaves, all theorem leaves and the independent nonsquare roots")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("jobs must be positive")
    root = args.root.resolve()
    data = root / "data/revision"
    archive = args.archive or data / "rank30_completed.zip"
    expected_path = args.report or data / "rank30_audit.json"
    sys.path.insert(0, str(root))
    import next_rank
    with TemporaryDirectory() as tmp:
        executables, hashes = next_rank.executables(Path(tmp))
        next_rank.verify(archive, executables, hashes, args.replay, args.jobs)
    if args.replay:
        subprocess.run([sys.executable, str(root / "verify_nonsquare.py"),
                        "--ranks", "30", "--replay"], check=True)
    stats = summarize(archive, data)
    require(stats == json.loads(expected_path.read_bytes()), "rank-30 audit report disagrees")
    print(f"Checkpoint {stats['checkpoint_profiles']:,} + new {stats['additional_profiles']:,} "
          f"= {stats['profiles']:,} distinct degree sequences; "
          f"{stats['nonsquare_independent_matches']} nonsquare matches; zero survivors.")


if __name__ == "__main__":
    main()
