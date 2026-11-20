"""Fresh, source-bound square-dimension enumeration at ranks 1 through 44.

Public interface used by the local runner:

* ``tasks(rank, builddir)`` regenerates the complete native prefix partition.
* ``command(builddir, task, rank)`` returns one native search command.
* ``exact_rows(item, task, rank)`` exhausts one arbitrary-precision handoff.
* ``row_failure(row, task, rank)`` has no additional terminal obstruction;
  the runner's terminal module owns the shared and square terminal tests.

Compile ``NATIVE_SOURCE`` as C++17, with the package's ``engine`` directory
on the include path, to ``builddir / BINARY_NAME``. The native protocol is
newline-delimited JSON: complete degree lists or {"exact": prefix} records,
followed by the mandatory DONE receipt on stderr. A native interruption or
an unfinished Python continuation never supplies completed task coverage.

No historical degree list, completion receipt, timing archive, or prior
coverage file is opened here. Structural preprocessing is a separate fresh
preflight. The optimized rank-44 continuation uses its high-four input;
lower ranks use the conservative unbounded continuation in this package.
"""

import json
from math import isqrt
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for directory in (ROOT / "engine", HERE):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

MAX_RANK = 44
BINARY_NAME = "square"
NATIVE_SOURCE = HERE / "square44.cpp"
PLANNER_DEPTH = 24
PLANNER_GROUP = 32
EXPECTED_TASK_COUNT = 1452
OBSTRUCTIONS = frozenset()


def _rank(rank):
    if type(rank) is not int or not 1 <= rank <= MAX_RANK:
        raise ValueError("the square adapter supports ranks 1 through 44")


def start_pairs(rank):
    """Arithmetic square starts before the necessary prefix exclusions."""
    _rank(rank)
    if rank < 18:
        return ()  # B>=4 and B^2<=rank-2.
    return tuple((b, 1) for b in range(4, isqrt(rank - 2) + 1))


def _prefix(prefix, rank):
    if (not isinstance(prefix, list) or not 1 <= len(prefix) <= rank - 3
            or any(type(b) is not int or b <= 0 for b in prefix)
            or prefix != sorted(prefix)
            or (prefix[0], 1) not in start_pairs(rank)):
        raise ValueError("invalid square task prefix")
    return prefix


def _task(task, rank):
    _rank(rank)
    if not isinstance(task, dict) or task.get("q") != 1:
        raise ValueError("expected a square task")
    prefix = _prefix(task.get("prefix"), rank)
    if type(task.get("B")) is not int or task["B"] != prefix[0]:
        raise ValueError("square task root differs from its prefix")
    return prefix


def tasks(rank, builddir):
    """Regenerate the complete disjoint partition from the native source.

    The native planner performs a complete depth-24 traversal and groups
    adjacent child subtrees in batches of at most 32 frontier leaves. Its
    completion trailer is mandatory. Partial output is never a partition.
    """
    _rank(rank)
    if not start_pairs(rank):
        return []
    executable = Path(builddir).resolve() / BINARY_NAME
    result = subprocess.run(
        [str(executable), "--rank", str(rank), "--tasks-cost",
         str(PLANNER_DEPTH), str(PLANNER_GROUP)],
        capture_output=True, text=True, check=True,
    )
    jobs = []
    for line in result.stdout.splitlines():
        item = json.loads(line)
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError("invalid square planner record")
        prefix, cost = item
        _prefix(prefix, rank)
        if type(cost) is not int or cost <= 0:
            raise ValueError("invalid square planner cost")
        jobs.append((tuple(prefix), cost))
    if result.stderr.strip() != f"TASKS {len(jobs)}":
        raise ValueError("missing or inconsistent square planner completion")
    if rank == 44 and len(jobs) != EXPECTED_TASK_COUNT:
        raise ValueError("square partition differs from the rank-44 source")
    previous = None
    for prefix, _ in sorted(jobs):
        if previous is not None and prefix[:len(previous)] == previous:
            raise ValueError("duplicate or overlapping square task prefixes")
        previous = prefix
    return [dict(q=1, B=prefix[0], prefix=list(prefix), cost=cost)
            for prefix, cost in sorted(jobs, key=lambda item: (-item[1], item[0]))]


def command(builddir, task, rank):
    """Native command for precisely one prefix subtree."""
    prefix = _task(task, rank)
    return [str(Path(builddir).resolve() / BINARY_NAME), "--rank", str(rank),
            "--prefix", ",".join(map(str, prefix))]


def exact_rows(item, task, rank):
    """Exhaust one width handoff using Python's unbounded integers."""
    _task(task, rank)
    if not isinstance(item, dict) or "exact" not in item:
        raise ValueError("invalid square exact-handoff record")
    if rank == 44:
        from exact_square44 import complete
    else:
        from exact_square_legacy import complete
    yield from complete(item, task, rank)


def row_failure(row, task, rank):
    """No adapter-specific rejection; terminal.py owns all terminal checks."""
    _task(task, rank)
    return None
