"""Two finite checks for the interaction of prime-denominator packets."""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys


def determinant(matrix):
    if not matrix:
        return 1
    return sum((-1)**j * entry * determinant([row[:j] + row[j + 1:] for row in matrix[1:]])
               for j, entry in enumerate(matrix[0]))


def joint_gram():
    signs = list(product((-1, 1), repeat=2))
    packets = [(80, (1, -1, -2, 0)), (112, (1, 0, 0, -1))]
    tested = 0
    survivors = []
    for n1 in range(8):
        for n2 in range(8 - n1):
            for n3 in range(8 - n1 - n2):
                counts = (n1, n2, n3, 7 - n1 - n2 - n3)
                seven = [(35 * n, (1, a, 0, u)) for n, (a, u) in zip(counts, signs)]
                for a, y, u in product((-1, 1), (-2, 0, 2), (-1, 1)):
                    tested += 1
                    vectors = packets + seven + [(35, (1, a, y, u))]
                    gram = [[560 * (i == j) - sum(n * v[i] * v[j] for n, v in vectors)
                             for j in range(4)] for i in range(4)]
                    minors = (determinant([[gram[i][j] for j in ix] for i in ix])
                              for k in range(1, 5) for ix in combinations(range(4), k))
                    if all(d >= 0 for d in minors):
                        survivors.append(((a, y, u), counts))
                        kernel = (6, -1, -7, 1)
                        assert all(sum(x * z for x, z in zip(row, kernel)) == 0 for row in gram)
    assert tested == 1440 and survivors == [((-1, 2, -1), (0, 1, 2, 4))]


def mixed_five(source):
    sys.path.insert(0, str(source))
    from large_four import profiles
    for fours, count in [(10, 5), (11, 1)]:
        lists = list(profiles(44, fours, False, (3, 5, 7), (0, 1, 2)))
        assert len(lists) == count
        for odd, even, other, mass in lists:
            assert odd == (1, 2, 4) and 5 in even and 6 in even
            if fours == 11:
                assert even == (3, 5, 6) and not other
                assert 1 + 11 + 16 + 7 + 3 + 6 == 44
                assert sum(bool(mask & 4) for mask in odd + even) == 3 < 6
    floors = (216, 40, 56, 756, 140, 7560)
    bound = Fraction(10, 16) + Fraction(16, 100) + Fraction(3, 81) + Fraction(1, 7) + Fraction(2, 64)
    bound += sum(Fraction(1, b*b) for b in floors)
    assert bound == Fraction(9498569, 9525600) < 1


if __name__ == '__main__':
    joint_gram()
    mixed_five(Path(sys.argv[1]).resolve())
    print('Five and seven Gram matrices and both mixed-five signatures checked.')
