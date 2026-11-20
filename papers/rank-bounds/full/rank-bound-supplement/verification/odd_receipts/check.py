#!/usr/bin/env python3
"""Check stored core odd execution records and regenerate the square census.

No final odd continuation is executed. The certificate binds source/result
bytes and verifies record accounting, not the mathematical exclusions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def normalized(value):
    return json.loads(json.dumps(value))


def main():
    if not __debug__:
        raise SystemExit('Run without -O or -OO.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    sys.path.insert(0,str(root/'odd_square'))
    from profiles import RANK, run
    starts, rows = run()
    census = json.loads((root/'odd_square/census.json').read_text())
    receipt = json.loads((root/'odd_square/certificate.json').read_text())
    assert RANK==209 and starts==16512 and normalized(rows)==census
    assert len(rows)==len(receipt['results'])==receipt['reservations']==1461
    assert len({json.dumps(row,sort_keys=True) for row in census})==len(census)
    assert receipt['rank_ceiling']==RANK and receipt['prime_profiles']==starts
    assert receipt['unresolved']==[]
    assert all(reason=='excluded' and type(nodes) is int and nodes>=0
               for reason,nodes in receipt['results'])
    isolated = sum(not row[2] for row in census)
    isolated_nodes = sum(n for row,(_,n) in zip(census,receipt['results']) if not row[2])
    mixed_nodes = sum(n for row,(_,n) in zip(census,receipt['results']) if row[2])
    assert isolated==receipt['isolated']==1151
    assert isolated_nodes==receipt['isolated_nodes']==6721348
    assert mixed_nodes==receipt['mixed_nodes']==2753944
    nonsquare = json.loads((root/'odd_nonsquare/receipt.json').read_text())
    assert nonsquare['rank_cap']==211 and nonsquare['solutions']==[]
    assert len(nonsquare['counts'])==14
    assert sum(v[0] for v in nonsquare['counts'].values())==685040
    assert all(all(type(v) is int and v>=0 for v in counts)
               and counts==sorted(counts,reverse=True)
               for counts in nonsquare['counts'].values())
    assert len(nonsquare['height_survivors'])==4
    assert nonsquare['terminal_denominators']=={'36':[[15,21,33]],'38':[],'40':[],'42':[]}
    assert sum(nonsquare['tail_nodes_by_remaining_pairs'])==72
    paths = [path for folder in ('odd_square','odd_nonsquare')
             for path in sorted((root/folder).iterdir())
             if path.suffix in ('.py','.json')]
    record = dict(format='odd-record-accounting-v1',
        scope='Source/result binding and regenerated square profile census; final odd recursions not rerun.',
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_and_result_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in paths},
        square=dict(rank_ceiling=209,prime_profiles=starts,reservations=len(rows),
                    isolated=isolated,mixed=len(rows)-isolated,
                    stored_isolated_nodes=isolated_nodes,stored_mixed_nodes=mixed_nodes,
                    stored_unresolved=0),
        nonsquare=dict(rank_ceiling=211,supports=14,stored_starting_profiles=685040,
                       stored_height_survivors=4,stored_tail_nodes=72,stored_solutions=0))
    text = json.dumps(record,sort_keys=True,indent=2)+'\n'
    destination = here/'certificate.json'
    if args.write:
        destination.write_text(text)
    else:
        assert destination.read_text()==text, 'Odd source/result binding or accounting differs'
    print('Regenerated 16,512 prime profiles and all 1,461 square reservations;')
    print('validated stored odd outcome accounting and all core source/result hashes.')
    print('The final odd continuation recursions were not repeated.')


if __name__=='__main__':
    main()
