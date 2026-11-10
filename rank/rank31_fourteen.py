#!/usr/bin/env python3
"""Verify the exact rank-31 root with fourteen leading denominator fours."""
import argparse
from collections import Counter
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import next_rank as nr
from generator_activity import mass_obstruction

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'/'rank31'
KEY = (1, (4,)*14)


def stored(sources=None):
    path = DATA/'fourteen_31.jsonl'
    ledger = json.loads((DATA/'fourteen_31.json').read_bytes())
    rows = nr.emit_rows(path.read_bytes().splitlines(),31,KEY)
    nr.require(ledger['rank']==31 and ledger['square_class']==1 and ledger['prefix']==list(KEY[1]),'wrong root')
    nr.require(ledger['raw_sha256']==nr.digest(path.read_bytes()),'root list checksum disagrees')
    nr.require(rows==sorted(set(rows)) and len(rows)==ledger['profiles']==251,'invalid root list')
    reasons = Counter()
    for row in rows:
        # The input has rank 31, square dimension and these fourteen leading
        # denominators, as independently checked by emit_rows above.
        hit = mass_obstruction(row)
        nr.require(hit is not None,'unexcluded row in the fourteen-fours root')
        reasons[str(hit[0])] += 1
    nr.require(dict(reasons)==ledger['activity_exclusions']=={'3':133,'5':118},'activity classifications disagree')
    if sources is not None:
        nr.require(ledger['sources']==sources,'mismatching root generator')
    return ledger,rows


def replay(exe, sources):
    ledger,expected = stored(sources)
    nodes,raw = nr.exact_run(exe['exact'],'shard',31,['1',nr.prefix_text(KEY[1])])
    rows = sorted(nr.emit_rows(raw,31,KEY))
    nr.require(nodes==ledger['gmp_nodes'] and rows==expected,'GMP replay disagrees')
    return ledger


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',action='store_true',help='regenerate the entire GMP root')
    args=ap.parse_args()
    if args.replay:
        import frontier31
        with TemporaryDirectory() as tmp:
            exe, hashes=frontier31.build_second(Path(tmp));ledger=replay(exe,hashes)
    else:
        ledger,_=stored()
    print(f"Rank 31, fourteen fours: {ledger['profiles']} types, 133 excluded at three and 118 at five.")


if __name__=='__main__':
    main()
