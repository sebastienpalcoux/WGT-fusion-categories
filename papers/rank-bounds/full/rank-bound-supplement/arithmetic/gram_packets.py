"""Check the column Gram matrices for simultaneous five- and seven-packets."""
from itertools import combinations, product


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


if __name__ == '__main__':
    joint_gram()
    print('Checked all 1,440 Gram matrices and the unique kernel witness.')
