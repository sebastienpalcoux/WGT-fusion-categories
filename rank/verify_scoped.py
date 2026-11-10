#!/usr/bin/env python3
"""Independent entry points for the reduced-table, paired-odd and formal-codegree searches."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from verify import CODE, DATA, compare, compile_tool, reduced_verdict, require

def reduced(replay, work):
    base = json.loads((DATA/"reduced14.json").read_bytes())
    short = base["row_feasible"]
    require(len(short) == 8 and all(1 + sum(x*x for x in ds) == q*m*m and
            2*q*m*m < 29160000 for q,m,ds in short), "reduced bounds through 14")
    if replay:
        path = work/"small.json"
        subprocess.run([sys.executable, CODE/"reduced_types.py", "14", path], check=True)
        require(json.loads(path.read_bytes()) == base, "reduced recursion through 14")
        exe = compile_tool("reduced_four", work)
    for rank, expected in [(15, {"BDP":6}), (16, {"BDP":133,"row gcd":104,"four":18,"three":1})]:
        name = f"four{rank}.jsonl"
        rows = [json.loads(line) for line in (DATA/name).read_bytes().splitlines()]
        require(all(len(d) == rank and d == sorted(d) and d.count(1)==1 for d in rows)
                and len({tuple(x) for x in rows}) == len(rows), "reduced profile format")
        require(Counter(map(reduced_verdict, rows)) == expected, "reduced terminal screen")
        if replay:
            p = subprocess.run([exe, str(rank)], capture_output=True, check=True)
            require(p.stderr.decode().splitlines()[-1].startswith(f"Complete rank {rank},"),
                    "missing reduced traversal completion marker")
            compare(p.stdout, name)
    print("Reduced half-ranks through 14: 8 profiles; 15: 6; 16: 256. All excluded.")

def paired(work):
    for rank, name, args in [(56,"odd56.json",[]),(74,"nonsquare74.json",["--omit-square-class-one"])]:
        out = work/name
        subprocess.run([sys.executable, CODE/"odd_frontier.py",str(rank),out,*args],check=True)
        data=json.loads(out.read_bytes())
        require(data == json.loads((DATA/"revision"/name).read_bytes()) and not data["row_feasible"],
                f"paired odd search through {rank}")
        print(f"{name}: {data['visited']} prefixes, zero degree lists.")

def formal(replay, work):
    name="formal_residual.txt"
    if replay:
        exe=compile_tool("formal_sieve",work)
        compare(subprocess.check_output([exe]),name)
    subprocess.run([sys.executable, CODE/"formal_check.py", DATA/name], check=True)

def main():
    if not __debug__:
        raise RuntimeError("Run the release verification with Python assertions enabled.")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("suite",choices=["reduced","paired-odd","formal"])
    p.add_argument("--replay",action="store_true",help="regenerate reduced or formal input; paired-odd always regenerates")
    args=p.parse_args()
    with TemporaryDirectory() as tmp:
        work=Path(tmp)
        if args.suite=="reduced":reduced(args.replay,work)
        elif args.suite=="formal":formal(args.replay,work)
        else:paired(work)

if __name__=="__main__":
    main()
