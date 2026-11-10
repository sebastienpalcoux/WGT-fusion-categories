"""Run a complete exact-rank search on multiple cores; preserve verified shards."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix+".new")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True)+"\n")
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("rank", type=int)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--depth", type=int, help="split after this many nonunit denominators")
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--prepare-only", action="store_true", help="produce the list of expected shards")
    p.add_argument("--verify-only", action="store_true", help="verify every shard and combine")
    p.add_argument("--only-q", type=int, help="compute one multiplier; global completion still requires every shard")
    p.add_argument("--root", nargs=2, type=int, action="append", metavar=("Q", "B"),
                   help="restrict the certificate to this (q, first denominator) root; repeatable")
    args = p.parse_args()
    if args.jobs < 1 or not 6 <= args.rank <= 48:
        p.error("jobs must be positive and rank must lie in 6..48")
    depth = args.depth if args.depth is not None else max(1,args.rank-7)
    if not 1 <= depth <= args.rank-3:
        p.error("depth must lie in 1..rank-3")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    # One coordinator per directory, including on restarts.
    lock = (work/".lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    source = Path(__file__).with_name("search.cpp").resolve()
    source_sha = digest(source)
    manifest = work/"manifest.tsv"
    config = {"rank":args.rank,"depth":depth,"source_sha256":source_sha}
    roots = sorted(set(map(tuple, args.root or [])))
    if roots:
        if any(q < 1 or b < 1 for q, b in roots):
            p.error("root coordinates must be positive")
        config["roots"] = list(map(list, roots))
    config_path = work/"search.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError("rank, split depth or source changed: choose a new work directory")
    atomic_json(config_path,config)
    exe = work/"search"
    build_file = work/"build.json"
    build = json.loads(build_file.read_text()) if build_file.exists() else {}
    if not exe.exists() or build.get("source_sha256") != source_sha or digest(exe) != build.get("binary_sha256"):
        tmp = work/"search.new"
        subprocess.run(["c++","-O3","-std=c++17",str(source),"-lgmpxx","-lgmp","-o",str(tmp)],check=True)
        tmp.replace(exe)
        atomic_json(build_file, {"source_sha256":source_sha,"binary_sha256":digest(exe)})
    new_manifest = work/"manifest.new"
    manifest_mark = re.compile(rf"COMPLETE manifest rank={args.rank} nodes=(\d+) types=0 shards=(\d+)\n")
    with new_manifest.open("w") as out:
        made = subprocess.run([exe,"--manifest",str(args.rank),str(depth)],
                              stdout=out,stderr=subprocess.PIPE,text=True,check=True)
    match = manifest_mark.fullmatch(made.stderr)
    rows = []
    seen = set()
    for line in new_manifest.read_text().splitlines():
        q, prefix = line.split("\t")
        coords = tuple(map(int,prefix.split(",")))
        if len(coords) != depth or (q,prefix) in seen or list(coords) != sorted(coords):
            raise ValueError("malformed or duplicate shard in regenerated manifest")
        seen.add((q,prefix))
        rows.append((q,prefix))
    if not match or int(match[2]) != len(rows):
        raise ValueError("manifest generation did not finish")
    if roots:
        rows = [(q,prefix) for q,prefix in rows if (int(q),int(prefix.split(",")[0])) in roots]
        new_manifest.write_text("".join(f"{q}\t{prefix}\n" for q,prefix in rows))
    if manifest.exists() and digest(manifest) != digest(new_manifest):
        raise ValueError("saved manifest disagrees with regenerated one")
    new_manifest.replace(manifest)
    print(f"Manifest contains {len(rows)} expected shards.",flush=True)
    if args.prepare_only:
        return
    shard_dir = work/"parts"
    shard_dir.mkdir(exist_ok=True)
    mark = re.compile(rf"COMPLETE search rank={args.rank} nodes=(\d+) types=(\d+) shards=0\n")

    def record(i):
        return shard_dir/f"{i:08d}.jsonl", shard_dir/f"{i:08d}.receipt"

    def checked(i):
        target, receipt = record(i)
        if not target.exists() or not receipt.exists():
            return False
        try:
            r = json.loads(receipt.read_text())
            q,prefix = rows[i]
            if (r["q"],r["prefix"],r["source_sha256"]) != (q,prefix,source_sha):
                return False
            if r["sha256"] != digest(target):
                return False
            m = mark.fullmatch(r["completion"])
            if not m:
                return False
            with target.open("rb") as f:
                return sum(1 for _ in f) == int(m[2])
        except (ValueError, TypeError, KeyError, OSError):
            return False

    def run_one(i):
        q,prefix = rows[i]
        target,receipt = record(i)
        tmp = target.with_suffix(".partial")
        with tmp.open("w") as out:
            proc = subprocess.run([exe,"--shard",str(args.rank),q,prefix],
                                  stdout=out,stderr=subprocess.PIPE,text=True,check=True)
        m = mark.fullmatch(proc.stderr)
        if not m:
            raise ValueError(f"missing completion marker for shard {i}")
        with tmp.open("rb") as f:
            if sum(1 for _ in f) != int(m[2]):
                raise ValueError(f"candidate count mismatch in shard {i}")
        tmp.replace(target)
        atomic_json(receipt, {"q":q,"prefix":prefix,"source_sha256":source_sha,
                              "completion":proc.stderr,"sha256":digest(target)})

    if not args.verify_only:
        selected = [i for i,(q,_) in enumerate(rows)
                    if args.only_q is None or int(q)==args.only_q]
        todo = [i for i in selected if not checked(i)]
        if todo:
            (work / "complete.json").unlink(missing_ok=True)
        done = len(selected) - len(todo)
        print(f"Rank {args.rank}: {done}/{len(selected)} shards already verified.", flush=True)
        last_report = time.monotonic()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending = [pool.submit(run_one, i) for i in todo]
            try:
                for finished in as_completed(pending):
                    finished.result()
                    done += 1
                    if time.monotonic() - last_report >= 15:
                        print(f"Rank {args.rank}: {done}/{len(selected)} shards verified.", flush=True)
                        last_report = time.monotonic()
            except BaseException:
                for task in pending:
                    task.cancel()
                raise
    missing = [i for i in range(len(rows)) if not checked(i)]
    if missing:
        (work/"complete.json").unlink(missing_ok=True)
        print(f"INCOMPLETE: {len(missing)} missing or invalid shards; first indices {missing[:20]}")
        if args.verify_only:
            sys.exit(1)
        return
    raw = set()
    nodes = 0
    for i,(q,prefix) in enumerate(rows):
        target,receipt = record(i)
        r = json.loads(receipt.read_text())
        nodes += int(mark.fullmatch(r["completion"])[1])
        with target.open() as f:
            for line in f:
                t = tuple(json.loads(line))
                if len(t)!=args.rank or t[0]!=1 or tuple(sorted(t))!=t or t in raw:
                    raise ValueError(f"malformed or repeated degree vector in shard {i}")
                raw.add(t)
    combined = work/"candidates.jsonl"
    temporary = combined.with_suffix(".jsonl.new")
    with temporary.open("w") as out:
        for t in sorted(raw):
            out.write(json.dumps(t,separators=(",",":"))+"\n")
    temporary.replace(combined)
    atomic_json(work/"complete.json",dict(config,shards=len(rows),candidate_count=len(raw),
                candidate_sha256=digest(combined),manifest_sha256=digest(manifest),
                nodes_summed_over_shards=nodes))
    scope = f"selected roots {roots} at rank" if roots else "rank"
    print(f"COMPLETE {scope} {args.rank}: {len(rows)} verified shards; {len(raw)} degree vectors.")


if __name__ == "__main__":
    main()
