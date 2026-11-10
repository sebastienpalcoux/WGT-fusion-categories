#!/usr/bin/env python3
"""Verify ordinary ranks 24–27, the partial rank 28, and the rank-specific support bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
sys.path.insert(0, str(ROOT / "code"))
from rank_audit import audit, fast_report, read_rows
from principal_cover import row_cover, dyadic_trace_possible


def check(condition, message):
    if not condition:
        raise ValueError(message)


def source_name(rank):
    return "ordinary_rank24_25.cpp" if rank <= 25 else f"ordinary_rank{rank}.cpp"


def check_prefix_boxes():
    for file, program, args in (
        ("three-prime-roots.json", "rank_family.py", list(range(24, 29))),
        ("three-prime-prefix.json", "prefix_boxes.py", ["--cut",25,*range(24,29)]),
        ("q10-through36.json", "prefix_boxes.py", ["--cut",50,"--q",10,*range(2,37)]),
    ):
        expected = json.loads((DATA / "strategy" / file).read_text())
        run = subprocess.run([sys.executable, str(ROOT / "code" / program), *map(str,args)],
                             check=True, capture_output=True)
        check(json.loads(run.stdout) == expected, f"{file}: support-box replay differs")


def verify(rank, replay, tmp):
    data = DATA / "ordinary"
    rows = read_rows(data / f"rank{rank}.jsonl", rank)
    analysis = audit(rows, rank)
    check(analysis == json.loads((DATA/"strategy"/f"rank{rank}.json").read_text()),
          f"rank {rank}: audit differs")
    check(not analysis["survivors"] and "open" not in analysis["family_strategy_first_exclusion"],
          f"rank {rank}: at least one pipeline leaves an open type")
    pattern = rf"Complete rank {rank}, nodes (\d+), types (\d+), seconds [\d.eE+-]+"
    saved = re.fullmatch(pattern, (data / f"rank{rank}.log").read_text().splitlines()[-1])
    check(saved and int(saved[2]) == len(rows), "invalid stored completion marker")
    if rank in (26, 27):
        sq16 = [v for v in rows if sum(d*d for d in v) == 16*v[-1]**2]
        check(sum(row_cover(v) == "lcm" for v in sq16) == {26:9,27:85}[rank],
              "principal lcm count differs")
        covered = [v for v in sq16 if row_cover(v) == 'passes']
        norms = [v for v in covered if dyadic_trace_possible(v,use_modular=False)]
        final = [v for v in norms if dyadic_trace_possible(v)]
        check((len(covered)-len(norms),len(norms)-len(final),len(final)) ==
              {26:(10,2,0),27:(305,41,1)}[rank], "dyadic row tests differ")
    if replay:
        name = source_name(rank)
        exe = tmp / Path(name).stem
        if not exe.exists():
            subprocess.run(["c++","-O3","-std=c++17",str(ROOT/"code"/name),"-o",str(exe)],check=True)
        run = subprocess.run([str(exe),str(rank),"8"], capture_output=True,check=True)
        marker = re.fullmatch(pattern,run.stderr.decode().splitlines()[-1])
        generated = sorted(json.loads(row) for row in run.stdout.splitlines())
        check(marker and marker.groups()==saved.groups() and generated==rows, "complete enumeration differs")
    print(f"Rank {rank}: {len(rows)} profiles, zero survivors, {saved[1]} recursion states.",flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks",type=int,nargs="+",choices=range(24,28),default=list(range(24,28)))
    parser.add_argument("--replay",action="store_true",help="regenerate all selected profiles from the rank alone")
    args = parser.parse_args()
    for line in (DATA/"frontier_sha256.txt").read_text().splitlines():
        digest, path = line.split("  ",1)
        check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest, f"digest differs: {path}")
    with tempfile.TemporaryDirectory() as d:
        for r in args.ranks:
            verify(r,args.replay,Path(d))
    partial=DATA/"frontier"
    rows=read_rows(partial/"rank28_partial.jsonl",28)
    receipt=json.loads((partial/"rank28_partial.json").read_text())
    check(receipt["complete"] is False and sum(t["complete"] for t in receipt["roots"])==13,
          "rank-28 frontier was expected to be partial")
    incomplete = {(t["q"],t["b"]) for t in receipt["roots"] if not t["complete"]}
    check(incomplete=={(1,3),(1,4)} and len(rows)==414, "unexpected rank-28 frontier")
    check(fast_report(rows,28) == json.loads((partial/"rank28_screen.json").read_text()),"partial screening differs")
    check_prefix_boxes()
    subprocess.run([sys.executable,str(ROOT/"code/q10_crosscheck.py")],check=True)
    print("Verified the partial rank-28 screen and the rank-bounded three-prime families.",flush=True)


if __name__ == "__main__":
    main()
