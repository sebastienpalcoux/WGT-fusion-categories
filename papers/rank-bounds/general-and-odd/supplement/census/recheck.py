"""Replay the general-rank census and verify its individual and aggregate certificates."""
import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import audit
from terminal_arithmetic import STAGES, obstruction

HERE = Path(__file__).resolve().parent


def replay(archive):
    checked = audit.audit(archive, recheck_prime=False)
    if not checked["empty_complete_census"]:
        raise ValueError("The archive does not contain a complete census.")
    table = []
    ledger = bytearray()
    with ZipFile(archive) as saved:
        for rank in range(1, 45):
            prefix = f"ranks/{rank:02d}/"
            tasks = json.loads(saved.read(prefix+"manifest.json"))["tasks"]
            emitted, reasons = Counter(), Counter()
            for raw in gzip.decompress(saved.read(prefix+"raw-census.jsonl.gz")).splitlines():
                r, index, _, _, degrees = json.loads(raw)
                # audit_rank has verified each record's rank, task, degree identities,
                # denominator prefix, top multiplicity and binding to its completion receipt.
                stage = obstruction(degrees, tasks[index])
                if stage not in STAGES:
                    raise ValueError(f"Unexcluded list at rank {r} in task {index}: {degrees}")
                emitted["square" if tasks[index]["q"] == 1 else "nonsquare"] += 1
                reasons[stage] += 1
                record = [rank, index, sha256(raw).hexdigest(), stage]
                ledger.extend(json.dumps(record, separators=(",", ":")).encode()+b"\n")
            line = {"rank": rank, "square": emitted["square"], "nonsquare": emitted["nonsquare"]}
            survivors = sum(emitted.values())
            for stage in STAGES:
                survivors -= reasons[stage]
                line[stage] = survivors
            table.append(line)
    totals = {key: sum(line[key] for line in table) for key in ("square","nonsquare",*STAGES)}
    counts = {"archive_sha256": sha256(archive.read_bytes()).hexdigest(),
              "ledger_sha256": sha256(ledger).hexdigest(), "stages": list(STAGES),
              "totals": totals, "table": table}
    return counts, bytes(ledger)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate both certificates")
    parser.add_argument("--archive", type=Path, default=HERE/"reference-census.zip")
    args = parser.parse_args()
    counts, rows = replay(args.archive)
    receipt = HERE/"rank-counts.json"
    ledger = HERE/"rank-dispositions.jsonl.gz"
    if args.write:
        receipt.write_text(json.dumps(counts, indent=2)+"\n")
        ledger.write_bytes(gzip.compress(rows, mtime=0))
    elif json.loads(receipt.read_text()) != counts or gzip.decompress(ledger.read_bytes()) != rows:
        raise ValueError("The per-rank or per-record certificate differs from the replay.")
    print(json.dumps(counts["totals"]))


if __name__ == "__main__":
    main()
