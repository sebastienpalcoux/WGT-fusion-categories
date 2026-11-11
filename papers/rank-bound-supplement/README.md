# Reproducing the rank bound

The verifier enumerates degree lists for a hypothetical non-pointed simple integral modular category. The paper excludes ranks 1–9 analytically. Its computational table includes rank 9 as a cross-check. For each rank from 9 through 31, the verifier generates the square and nonsquare global-dimension trees, checks each completed list, and counts it under its first obstruction. A surviving list would be written out in full.

## Run from this directory

Use Python 3.9 or later on Linux (including WSL) or macOS, and a C++17 compiler with unsigned 128-bit integers, such as GCC or Clang. The Python modules use only the standard library. Set `CXX` if your compiler is named something other than `c++`.

To reproduce the entire computation, choose a fresh directory name for `--work` and adjust `--jobs` to the number of processor cores you want to use.

```sh
python3 verify.py --jobs 4 --work rank-results
```

The driver compiles the generators, prints progress, and writes one checkpoint per completed task. To resume after an interruption, repeat the same command. The final table gives the number of degree lists entering the terminal tests and the number surviving. A zero in either dimension column means its enumeration generated no complete profile at that rank. The complete output is `rank-results/certificate.json.gz`. The published run covered 3,046 tasks and tested 453,041 profiles (452,592 square and 449 nonsquare). Every profile was excluded.

The same driver accepts any positive rank. For example, the following command explores rank 32 and stores its checkpoints separately.

```sh
python3 verify.py --rank 32 --jobs 4 --work rank-32
```

Higher ranks are exploratory. The archive accompanying the theorem covers ranks 9–31. Larger searches can require substantially more time and memory, and any surviving profiles are candidates for further analysis. To run a bounded number of tasks in one sitting, add `--limit N`. A final certificate is produced when every task has a receipt.

## Check the published record

```sh
python3 verify.py --check certificate.json.gz
```

The SHA-256 digest of the published compressed archive is `384d93fa55992826bf3c9566f09087fd21aa6ed2161997246ef1043e5ba69c96`. This quick check recompiles the generators, regenerates the task partition, checks all nine source hashes, and reconstructs the table from the archived receipts. To verify the enumerations themselves, use the full command above. A regenerated archive can also be compared with the published one in Python, independently of the gzip encoding.

```sh
python3 -c "import gzip, json; a,b=(json.load(gzip.open(p)) for p in ('certificate.json.gz','rank-results/certificate.json.gz')); assert a == b"
```

## Files and record format

- [verify.py](verify.py) constructs the manifest, dispatches tasks, checks incoming lists and saves the receipts. [primes.py](primes.py) supplies rank-dependent prime lists to the Python generators.
- [square.cpp](square.cpp) and [nonsquare.cpp](nonsquare.cpp) enumerate in the two square classes of the global dimension. Both receive the rank as an argument. Square tasks fix an initial denominator prefix. Nonsquare tasks fix the squarefree part, the first denominator and, for the two most expensive roots, the maximal-degree multiplicity.
- [square_exact.py](square_exact.py) and [nonsquare_exact.py](nonsquare_exact.py) use Python integers for branches needing more precision. They also handle whole tasks when the rank exceeds a native mask capacity (64 for square dimension and 136 for nonsquare dimension). A native generator emits one JSON record per line, either a sorted degree list or an object with an `exact` denominator prefix. Nonsquare continuations also record the prime support. The final diagnostic reports the native node, profile and handoff counts. An inconsistent stream triggers exact replay of the whole affected task.
- [screen.py](screen.py) applies the degree, crossed-fiber, activity, balance and dyadic tests, calling [trace.py](trace.py) for orbit traces and [multiplication.py](multiplication.py) for multiplication constraints.
- [certificate.json.gz](certificate.json.gz) is a single gzip-compressed JSON object. Its `manifest` records the ranks, SHA-256 source hashes and ordered tasks. Each task has the form `[rank, case, key]`. The corresponding entry in `receipts` records its native nodes, exact handoffs, profile count, counts by first obstruction, and all survivors. A null node count means that the whole task was computed with Python integers. The `summary` gives the totals by rank. An empty exclusion count means that no generated degree list reached a terminal obstruction in that task.
