"""Verify the current uniform odd-square reservation screen.

This reads the original saved inputs, preserves their production sources and
receipts, and performs no profile regeneration or Gauss search. The rank209
input remains exploratory. Run from the supplement root:

    python3 verification/odd_uniform/check.py

After an intentional checker change, use --write-receipt to update only this
checker's derived receipt.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
from component import component_envelope
from inertia import BASELINE, HERE, SUPPLEMENT, extra_screen, file_hash, old


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical_result(value):
    """Use JSON key/container types in both write and read-only replay."""
    return json.loads(json.dumps(value))


def evaluate():
    provenance = BASELINE.parent/'provenance.json'
    manifest = json.loads(provenance.read_text())
    require(manifest['schema'] == 'odd-neutral-screen-inputs-v1', 'unknown baseline schema')
    require(file_hash(SUPPLEMENT/manifest['baseline_receipt']) == manifest['baseline_receipt_sha256'],
            'original odd-square receipt changed')
    for name, digest in manifest['baseline_source_sha256'].items():
        require(file_hash(SUPPLEMENT/name) == digest, 'original source changed: '+name)
    screens = []
    for rank in (201, 209):
        metadata = manifest['inputs'][str(rank)]
        path = SUPPLEMENT/metadata['path']
        require(file_hash(path) == metadata['sha256'], 'saved input changed')
        data = json.loads(path.read_text())
        if rank == 209:
            require(data['count'] == metadata['prime_profiles'], 'profile count changed')
        rows = data if rank == 201 else data['rows']
        require(len(rows) == metadata['reservations'], 'reservation count changed')
        counts, outcomes, retained = Counter(), [], []
        for index, row in enumerate(rows):
            old.validate_row(row, rank)
            verdict = extra_screen(row, rank)
            if verdict['keep']:
                verdict = component_envelope(row, rank)
            reason = 'retained' if verdict['keep'] else verdict['reason']
            counts[reason] += 1
            outcomes.append({'index': index, **verdict})
            if verdict['keep']:
                retained.append(index)
        graph_counts = Counter('edge35' if rows[i][2] else 'edgeless' for i in retained)
        screens.append({'rank': rank, 'input_path': metadata['path'],
                        'input_sha256': file_hash(path),
                        'baseline_prime_profiles': metadata['prime_profiles'],
                        'input_reservations': len(rows), 'counts': dict(counts),
                        'retained': len(retained), 'retained_graphs': dict(graph_counts),
                        'retained_indices': retained, 'outcomes': outcomes})
    return canonical_result({
        'schema': 'odd-uniform-reservation-screen-v2',
        'scope': ('Exact necessary-condition screens of the original saved reservations. '
                  'No profile regeneration, Gauss search, rank209 exclusion, or '
                  'realizability assertion is made.'),
        'current_module_sha256': {name: file_hash(HERE/name) for name in (
            'check.py', 'inertia.py', 'component.py', 'residue.py',
            'global_weyl.py', 'coupled_weyl.py', 'test_conditions.py')},
        'baseline_checker_sha256': file_hash(BASELINE),
        'baseline_provenance_sha256': file_hash(provenance),
        'baseline_receipt_sha256': manifest['baseline_receipt_sha256'],
        'baseline_source_sha256': manifest['baseline_source_sha256'],
        'final_gauss_searches': 0, 'screens': screens,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    result = evaluate()
    path = HERE/'current-inertia-and-component-receipt.json'
    if args.write_receipt:
        path.write_text(json.dumps(result, indent=1, sort_keys=True)+'\n')
    else:
        require(json.loads(path.read_text()) == result, 'uniform odd receipt mismatch')
    print(json.dumps({'status': 'passed', 'screens': [
        {key: s[key] for key in ('rank', 'input_reservations', 'counts', 'retained', 'retained_graphs')}
        for s in result['screens']], 'final_gauss_searches': 0}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
