"""Check the septic cutoff constants and count its overlap with the rank-30 upload."""
import argparse
import json
import math
import zipfile
from collections import Counter
from fractions import Fraction


def cutoff(rank, q, denominators):
    return rank <= 34 and q == 1 and denominators.count(4) >= 14 and 7 in denominators


def verify_constants():
    m = 5400
    p = 7
    c = Fraction(4*(p-1)**2, p**4)
    for j, inactive in [(2, 2), (4, 1), (5, 1)]:
        assert inactive*c > 1 - Fraction(14, 16) - Fraction(j, 49)
    assert 3*c > 1 - Fraction(14,16) - Fraction(3,49)
    assert Fraction(3,4) + Fraction(2,9) < 1
    assert 3*5 + 3 > 34 - 14 - 3 - 1
    assert 14*Fraction(1,16) + 3*Fraction(1,49) + 12*Fraction(1,196) == 1 - Fraction(1,392)
    assert Fraction(11, 196) + Fraction(2, 441) + Fraction(1, m*m) < Fraction(3, 49)
    delta = 1 - Fraction(14, 16) - Fraction(6, 49)
    assert delta == Fraction(1, 392)
    error = delta + Fraction(1, m)
    assert max(32*error, 98*error) < 1
    assert [(a, b) for a in range(-14, 15) for b in [-1, 6]
            if abs(Fraction(a, 16) + Fraction(b, 49)) <= error] == [(-2, 6)]
    assert Fraction(9, 196) > delta
    assert (26**2+47)//48 == 15


def project(filename):
    counts = Counter()
    shard_counts = Counter()
    profile_counts = Counter()
    with zipfile.ZipFile(filename) as archive:
        with archive.open('spec.json') as file:
            rank = json.load(file)['rank']
        if rank != 30:
            raise ValueError("Expected the rank-30 completed archive.")
        with archive.open('receipts.jsonl') as receipts:
            for line in receipts:
                shard = json.loads(line)
                q = shard['q']
                prefix = shard['prefix']
                hit = cutoff(rank, q, prefix)
                counts['shards'] += 1
                counts['native_nodes'] += shard['native_nodes']
                if hit:
                    counts['whole_shards'] += 1
                    counts['whole_shard_native_nodes'] += shard['native_nodes']
                    counts['whole_shard_gmp_nodes'] += shard['gmp_nodes']
                    shard_counts[prefix.count(7)] += 1
                for degrees in shard['rows']:
                    counts['profiles'] += 1
                    if q != 1:
                        assert not hit
                        continue
                    dimension = sum(d*d for d in degrees)
                    M = math.isqrt(dimension)
                    assert M*M == dimension and all(M % d == 0 for d in degrees)
                    denominators = [M//d for d in degrees]
                    applies = cutoff(rank, q, denominators)
                    if hit:
                        assert applies
                        counts['whole_shard_profiles'] += 1
                    if applies:
                        counts['cutoff_profiles'] += 1
                        profile_counts[denominators.count(7)] += 1
    return {'counts': dict(sorted(counts.items())),
            'full_shards_by_prefix_sevens': dict(sorted(shard_counts.items())),
            'profiles_by_exact_sevens': dict(sorted(profile_counts.items()))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rank30_zip', nargs='?')
    args = parser.parse_args()
    verify_constants()
    if args.rank30_zip:
        print(json.dumps(project(args.rank30_zip), indent=2))
    else:
        print("Septic cutoff constants verified.")
