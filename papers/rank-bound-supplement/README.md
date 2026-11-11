# Verification supplement for the rank bound

This directory accompanies *A rank bound for simple integral modular categories*. The single entry point verifies ranks through 31. Ranks 2–8 are excluded in the paper without enumeration. At ranks 9–30 it regenerates all degree profiles and checks the fixed reference counts. At rank 31 it uses the faster, resumable prefix partition and exact integer continuations. Both routes use the terminal necessary conditions in `screen.py`, `trace.py` and `multiplication.py`.

Run from this directory with Python 3.9 or later and a C++17 compiler (GCC or Clang). For the entire computation, including rank 31, use

```sh
python3 verify.py --jobs 6
```

To verify one rank, use `python3 verify.py --rank 27 --jobs 4` or `python3 verify.py --rank 31 --jobs 6`. For a short rank-31 smoke test in a fresh directory use `python3 verify.py --rank 31 --limit 15 --jobs 4`. A partial rank-31 run produces checkpoints and reports how many jobs remain. Rerun the command without `--limit` to complete it. Use `--work DIRECTORY` to select its checkpoint directory. The default is `rank31-results-fast/` here. Do not run two copies against the same directory at once. `--case square` or `--case nonsquare` selects a rank-31 dimension class.

At ranks 9–30, `square.cpp` and `nonsquare.cpp` implement the original rank dependent search. At rank 31, `square31.cpp` and `nonsquare31.cpp` implement the stronger prefix pruning described in `square-pruning.md` and `nonsquare-pruning.md`. When either reaches its native integer width, `square_exact.py` or `nonsquare_exact.py` completes the branch using Python integers. These specialized internal kernels share one interface and screening logic. Their different prefix cuts preserve the speed of the completed rank-31 run.

The rank-31 engine (`run_rank31.py`) writes a manifest with source hashes, one `job-*.jsonl` per completed task, `summary.json`, and `survivors.jsonl`. It checks that a task stream has a complete receipt before recording it. The public verifier additionally checks the full summary against the rank-31 row of `counts.tsv` and requires an empty survivor file. This takes place only after every job has completed. The 1,167 jobs comprise 1,137 square and 30 nonsquare jobs, producing 434,509 and 428 terminal profiles respectively, with no survivors. The completed job receipts, manifest and summary are included in `rank31-results-fast-output.tar.gz`. To inspect that completed run and verify its receipts, without recomputing finished jobs, use

```sh
tar -xzf rank31-results-fast-output.tar.gz
python3 verify.py --rank 31 --jobs 4
```

`check_six()` in `verify.py` gives an additional consistency check on the codegree-six subcase, which also has a short hand proof in the manuscript. The mathematical reasons for each generator cut and terminal exclusion are given in the manuscript and its linked pruning notes. The 31 run uses no input database or fixed list of degree profiles.
