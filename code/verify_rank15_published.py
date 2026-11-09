"""Verify the tally used solely in the additional rank-15 corollary."""
import ast
import hashlib
from collections import Counter
from pathlib import Path

data = (Path(__file__).resolve().parents[1] / "data" / "published-rank15-perfect.txt").read_bytes()
assert hashlib.sha256(data).hexdigest() == "04a25855339d84fed23fb0e67cdef00b7783e2f1816d83fe679e3e143ac75dbc"
types = [ast.literal_eval(line) for line in data.decode().splitlines()]
assert len(types) == len({tuple(d) for d in types}) == 2341
counts = Counter()
exceptions = []
for d in types:
    assert len(d) == 15 and d == sorted(d) and d[:2] == [1, d[1]] and d[1] > 1
    D = sum(d_i**2 for d_i in d)
    assert all(D % (d_i*d_i) == 0 for d_i in d)
    f = D // max(d)**2
    counts[f] += 1
    if f > 4:
        exceptions.append((D, d))
assert counts == {2: 459, 3: 24, 4: 1857, 5: 1}
assert exceptions == [(720, [1, 2, 2, 2, 3, 3, 3, 4, 4, 6, 6, 12, 12, 12, 12])]
print(counts, exceptions)
