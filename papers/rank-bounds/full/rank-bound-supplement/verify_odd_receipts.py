#!/usr/bin/env python3
"""Check every recorded odd-square outcome and regenerate its reservation census.

This validates the source-bound execution record. It does not replace the
proofs of the pruning conditions or re-execute the final arithmetic searches.
"""
from hashlib import sha256
from pathlib import Path
import json
import sys


def verify():
    root = Path(__file__).resolve().parent
    arithmetic = root / 'arithmetic'
    nonsquare = json.loads((arithmetic / 'odd-nonsquare-receipt.json').read_text())
    assert nonsquare['survivors'] == 0 and nonsquare['global_nonsquare_rank_cap'] == 208
    assert sha256((root / nonsquare['source_file']).read_bytes()).hexdigest() == nonsquare['source_sha256']
    assert sha256((root / 'verification/odd_theory/nonsquare-replay.log').read_bytes()).hexdigest() == nonsquare['stdout_sha256']
    receipt = json.loads((arithmetic / 'odd-square-receipt.json').read_text())
    assert receipt['status'] == 'complete'
    assert receipt['unresolved_reservations'] == []
    for name, digest in receipt['source_sha256'].items():
        assert sha256((root / name).read_bytes()).hexdigest() == digest, name
    case_file = arithmetic / 'odd-square-cases.json'
    outcome_file = arithmetic / 'odd-square-outcomes.jsonl'
    assert sha256(case_file.read_bytes()).hexdigest() == receipt['case_list_sha256']
    assert sha256(outcome_file.read_bytes()).hexdigest() == receipt['outcomes_sha256']
    cases = json.loads(case_file.read_text())
    outcomes = [json.loads(line) for line in outcome_file.read_text().splitlines()]
    count = receipt['reservations']
    assert len(cases) == len(outcomes) == count == receipt['excluded_reservations']
    assert len({json.dumps(case, sort_keys=True) for case in cases}) == count
    assert sorted(row['index'] for row in outcomes) == list(range(count))
    for row in outcomes:
        assert row['status'] == 'excluded', row
        assert row['kind'] == int(bool(cases[row['index']][2])), row
    assert sum(not case[2] for case in cases) == receipt['isolated_reservations']
    assert sum(bool(case[2]) for case in cases) == receipt['mixed_reservations']
    sys.path.insert(0, str(arithmetic))
    import odd_square
    profiles, regenerated = odd_square.run()
    assert profiles == receipt['prime_profiles']
    assert odd_square.RANK == receipt['rank_ceiling']
    assert json.loads(json.dumps(regenerated)) == cases
    print(f'Verified all {count} indexed exclusions, regenerated {profiles:,} '
          'prime profiles, and checked every production source hash.')


if __name__ == '__main__':
    verify()
