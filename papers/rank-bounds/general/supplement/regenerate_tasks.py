#!/usr/bin/env python3
"""Regenerate the forty-four general census partitions and compare every task."""
import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile


def verify(archive_path):
    with ZipFile(archive_path) as archive, TemporaryDirectory() as directory:
        temporary = Path(directory)
        source_dir = temporary / "source"
        expected = json.loads(archive.read("preparation.json"))["source_sha256"]
        for name, digest in expected.items():
            relative = PurePosixPath(name)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Invalid census source path")
            data = archive.read("source/" + name)
            if sha256(data).hexdigest() != digest:
                raise ValueError(f"Source mismatch at {name}")
            source = source_dir / relative
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_bytes(data)
        sys.path.insert(0, str(source_dir))
        import run

        build = temporary / "build"
        run.build(build)
        total = 0
        for rank in range(1, 45):
            run.RANK = rank
            tasks = run.plan(build)
            record = json.loads(archive.read(f"ranks/{rank:02d}/manifest.json"))
            if tasks != record["tasks"] or run.digest(tasks) != record["partition_sha256"]:
                raise ValueError(f"Task partition mismatch at rank {rank}")
            total += len(tasks)
            print(f"Rank {rank}: {len(tasks)} tasks", flush=True)
        print(f"All 44 partitions regenerated, {total} tasks.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    verify(parser.parse_args().archive)
