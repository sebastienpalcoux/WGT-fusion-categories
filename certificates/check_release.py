#!/usr/bin/env python3
"""Check the frozen file inventory and rerun any selected arithmetic comparisons."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    if not __debug__:
        raise RuntimeError("Run the release verification with Python assertions enabled.")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases", nargs="+", metavar="ID",
                   help="e.g. D00 D24, or all for every dimension/boundary comparison")
    args = p.parse_args()
    manifest = json.loads((HERE/"manifest.json").read_bytes())
    for name, expected in manifest["file_sha256"].items():
        raw = (ROOT/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f"damaged release file: {name}")
    ordinary = next(item for item in manifest["rank_calculations"] if item["id"] == "R01")
    supplement = ROOT / ordinary["cwd"]
    certificate = supplement / ordinary["certificate"]
    with gzip.open(certificate, "rt") as stream:
        record = json.load(stream)
    if record["manifest"]["ranks"] != list(range(9, 32)):
        raise ValueError("incomplete ordinary rank range")
    observed = {reason for receipt in record["receipts"] for reason in receipt["excluded"]}
    if observed != set(manifest["ordinary_recorded_terminal_reasons"]):
        raise ValueError("incomplete theorem map for ordinary row exclusions")
    if record["summary"] != ordinary["summary"]:
        raise ValueError("ordinary rank totals disagree with the release manifest")
    if record["manifest"]["sha256"] != ordinary["source_sha256"]:
        raise ValueError("ordinary source inventory disagrees with the release manifest")
    if record["manifest"].get("families") is not True:
        raise ValueError("full-article certificate must enable the dimension families")
    subprocess.run([sys.executable, *ordinary["command"][1:]], cwd=supplement, check=True)
    independent = ordinary["independent_mode"]
    with gzip.open(supplement / independent["certificate"], "rt") as stream:
        standalone = json.load(stream)
    if (standalone["manifest"].get("families") is not False
            or standalone["manifest"]["sha256"] != independent["source_sha256"]
            or standalone["summary"] != independent["summary"]):
        raise ValueError("self-contained proof's mode, sources or totals disagree")
    subprocess.run([sys.executable, *independent["command"][1:]], cwd=supplement, check=True)
    for stage in manifest["ordinary_filter_map"]:
        for name in stage["implementation"]:
            if name not in manifest["file_sha256"]:
                raise ValueError(f"missing source hash for filter implementation: {name}")
    print(f"Verified {len(manifest['file_sha256'])} program, input and receipt digests.", flush=True)
    if args.cases is None:
        return
    entries = {item["id"]: item for item in manifest["calculations"]}
    ids = list(entries) if args.cases == ["all"] else list(dict.fromkeys(args.cases))
    if any(k not in entries for k in ids):
        p.error("choose all or an ID from D00 through D26")
    with TemporaryDirectory() as temp:
        for key in ids:
            item = entries[key]
            cmd = list(item["command"])
            if cmd[0] == "./strip_extra":
                exe = Path(temp)/"strip_extra"
                subprocess.run(["c++", "-std=c++20", "-O3", ROOT/"code/strip_extra.cpp",
                                "-o", exe], check=True)
                cmd[0] = exe
            completed = subprocess.run(cmd, cwd=ROOT/"code", capture_output=True, text=True, check=True)
            actual = completed.stdout.strip()
            if cmd[1:2] == ["four_primes/new_odd.py"]:
                actual = re.sub(r"; [0-9.]+s$", "", actual)
            if actual != item["expected_stdout"].strip():
                raise ValueError(f"unexpected output for {key}: {actual}")
            print(f"{key}: expected output reproduced.", flush=True)


if __name__ == "__main__":
    main()
