#!/usr/bin/env python3
"""Index exact destinations, printed numbers and pages from the compiled article."""
import argparse
import csv
import json
from pathlib import Path
import re


def groups(text):
    """Read braced TeX groups without expanding commands or mistaking escaped braces."""
    result = []
    depth = 0
    start = 0
    for match in re.finditer(r"\\[{}\\]|[{}]", text):
        token = match[0]
        if token == "{":
            if depth == 0:
                start = match.end()
            depth += 1
        elif token == "}":
            depth -= 1
            if depth < 0:
                raise ValueError("unmatched closing brace")
            if depth == 0:
                result.append(text[start:match.start()])
    if depth:
        raise ValueError("unmatched opening brace")
    return result


def index(aux):
    rows = {}
    for line in aux.read_text().splitlines():
        if not line.startswith(r"\newlabel{"):
            continue
        key, body = groups(line)
        fields = groups(body)
        if len(fields) < 4:
            raise ValueError(f"incomplete hyperref label {key}")
        number, page, title, dest = fields[:4]
        if key in rows:
            raise ValueError(f"duplicate label {key}")
        rows[key] = (dest, number if not dest.startswith("section*.") else "", page, title)
    if not rows:
        raise ValueError("no hyperref labels found; run latexmk first")
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aux", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--manifest", type=Path, default=Path(__file__).with_name("manifest.json"))
    args = p.parse_args()
    rows = index(args.aux)
    manifest = json.loads(args.manifest.read_bytes())
    used = set()
    for group in ("calculations", "rank_calculations"):
        for item in manifest[group]:
            if not item.get("edition"):
                used.update(item["source_labels"])
                used.update(item.get("derivation_labels", []))
    for item in manifest.get("ordinary_filter_map", []):
        used.update(item["labels"])
    for labels in manifest.get("ordinary_recorded_terminal_reasons", {}).values():
        used.update(labels)
    missing = sorted(used - rows.keys())
    if missing:
        raise ValueError("missing theorem or filter labels: " + ", ".join(missing))
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_label", "pdf_destination", "printed_number", "page", "title_TeX"])
        for label, fields in rows.items():
            writer.writerow([label, *fields])
    print(f"{len(rows)} labels indexed; {len(used)} distinct computational anchors resolved.")


if __name__ == "__main__":
    main()
