"""Reproduce the seven-test comparison on the complete archived denominator census.

Fresh enumerations call the very same terminal classifier. This checker
validates the archive in terms of its records and hashes, without loading
or executing any of its embedded producer programs.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
import argparse
import gzip
import json

import audit
import terminal

HERE = Path(__file__).resolve().parent
RECEIPTS = HERE.parent / "verification" / "rank_tables"
STAGES = ("degree_slices", "principal_overlap", "fusion_phases", "prime_square_trace",
          "block_entry", "simultaneous_block_trace", "defect_one")


def replay(archive):
    integrity = audit.audit(archive, recheck_prime=False)
    if not integrity["empty_complete_census"]:
        raise ValueError("reference tasks or outcomes are incomplete")
    rows, ledger = [], []
    with ZipFile(archive) as saved:
        for rank in range(1, 45):
            prefix = f"ranks/{rank:02d}/"
            tasks = json.loads(saved.read(prefix+"manifest.json"))["tasks"]
            counts, first = Counter(), Counter()
            for raw in gzip.decompress(saved.read(prefix+"raw-census.jsonl.gz")).splitlines():
                r, index, backend, _, degrees = json.loads(raw)
                if rank != r:
                    raise ValueError("incorrect rank in source record")
                task = tasks[index]
                reason = terminal.reason_for(degrees, task, rank)
                column = "prime_square_trace" if reason == "four_prime_isolation" else reason
                if column not in STAGES:
                    raise ValueError(f"the stated seven stages do not exclude this record: {r}, {index}, {reason}")
                counts["square" if task["q"] == 1 else "nonsquare"] += 1
                first[column] += 1
                ledger.append([rank, index, sha256(raw).hexdigest(), reason])
            remaining = counts["square"] + counts["nonsquare"]
            record = {"rank": rank, "square": counts["square"], "nonsquare": counts["nonsquare"]}
            for stage in STAGES:
                remaining -= first[stage]
                record[stage] = remaining
            rows.append(record)
    encoded = b"".join((json.dumps(line,separators=(",", ":"))+"\n").encode() for line in ledger)
    sources = {str(path.relative_to(HERE)): sha256(path.read_bytes()).hexdigest()
               for path in sorted(HERE.rglob("*")) if path.suffix in audit.SOURCE_SUFFIXES and "__pycache__" not in path.parts}
    result = {"archive_sha256": sha256(archive.read_bytes()).hexdigest(),
              "current_sources": sources, "reference_totals": integrity["totals"],
              "ledger_sha256": sha256(encoded).hexdigest(), "table": rows,
              "totals": {s: sum(row[s] for row in rows) for s in ("square","nonsquare",*STAGES)}}
    return result, encoded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate the release records")
    args = parser.parse_args()
    result, records = replay(HERE/"reference-census.zip")
    receipt, ledger = RECEIPTS/"table1.json", RECEIPTS/"dispositions.jsonl.gz"
    if args.write:
        receipt.write_text(json.dumps(result, indent=2)+"\n")
        ledger.write_bytes(gzip.compress(records,mtime=0))
    elif json.loads(receipt.read_text()) != result or gzip.decompress(ledger.read_bytes()) != records:
        raise ValueError("Seven-test counts or row-level dispositions disagree with the saved results")
    print(json.dumps(result["totals"]))


if __name__ == "__main__":
    main()
