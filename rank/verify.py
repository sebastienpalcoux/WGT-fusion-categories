"""Shared I/O and compiler helpers for the scoped rank calculations."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CODE, DATA = ROOT/'code', ROOT/'data'
sys.path.insert(0, str(CODE))
from check_reduced import check as reduced_verdict


def require(ok,message):
    if not ok:
        raise ValueError(message)


def run(cmd):
    return subprocess.run(list(map(str,cmd)),check=True,capture_output=True,text=True)


def compile_tool(name,work):
    exe=work/name
    run(['c++','-O3','-std=c++17',CODE/(name+'.cpp'),'-o',exe])
    return exe


def compare(source,name):
    require(source==(DATA/name).read_bytes(),f'Regenerated output differs: {name}')


if __name__ == "__main__":
    raise SystemExit(
        "Use python3 papers/rank-bound-supplement/verify.py for the unified search; "
        "this module supplies the scoped-suite helpers."
    )
