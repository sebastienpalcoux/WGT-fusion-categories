"""Rebuild the odd-square census at rank ceiling 209 and compare every decision."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

from profiles import RANK, run, witnesses
from isolated import separate
from mixed import mixed_rows

HERE = Path(__file__).resolve().parent


def classify(row):
    if row[2]:
        return mixed_rows([row])[0]
    return separate(row, RANK, 1500000, witnesses)


def replay(jobs):
    starts, rows = run()
    if jobs == 1:
        outcomes = [classify(row) for row in rows]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            outcomes = list(pool.map(classify, rows, chunksize=4))
    unresolved = [(i, reason) for i, (reason, _) in enumerate(outcomes)
                  if reason != 'excluded']
    summary = {
        'rank_ceiling': RANK,
        'prime_profiles': starts,
        'reservations': len(rows),
        'isolated': sum(not row[2] for row in rows),
        'isolated_nodes': sum(nodes for row, (_, nodes) in zip(rows, outcomes) if not row[2]),
        'mixed_nodes': sum(nodes for row, (_, nodes) in zip(rows, outcomes) if row[2]),
        'unresolved': unresolved,
    }
    print(json.dumps(summary), flush=True)
    if unresolved:
        raise ValueError('There are unresolved reservations')
    return {'census.json': rows, 'certificate.json': dict(summary, results=outcomes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='write the computed census and certificate')
    parser.add_argument('--jobs', type=int, default=1, help='number of workers')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    for name, body in replay(args.jobs).items():
        data = json.dumps(body, separators=(',', ':')) + '\n'
        path = HERE / name
        if args.write:
            path.write_text(data)
        elif path.read_text() != data:
            raise ValueError(f'Certificate mismatch at {name}')


if __name__ == '__main__':
    main()
