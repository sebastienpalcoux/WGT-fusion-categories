"""Regression checks for the corner cases in the four/five-prime catalogue."""
from itertools import permutations
from many_prime_families import S3, S4, many_prime_family, dyadic_character_box, odd_residue_periods


def prime(n):
    return n >= 2 and all(n % d for d in range(2, n) if d*d <= n)


def test_quadratic_family():
    """Compare the closed condition with the divisor character test directly."""
    for t in range(3, 121):
        if not prime(t) or t % 8 != 1:
            continue
        residues = {i*i % t for i in range(1, t)}
        allowed = [j for j in range(1, t) if (t-1) % j == 0 or (t+1) % j == 0]
        q = next(q for q in range(3, 100, 2) if prime(q) and q in residues)
        # Pick w = 3 mod 8 so only the singleton t can trigger the test.
        w = next(w for w in range(3, 300, 8) if prime(w) and w != t
                 and w % t not in residues)
        assert dyadic_character_box(q, w, t) == all(j in residues for j in allowed)
    assert dyadic_character_box(19, 43, 73)
    assert 26**2 % 73 == 19 and pow(43, 36, 73) == 72
    assert not dyadic_character_box(43, 19, 73)


def test_catalogue():
    assert len(S3) == 9 and len(S4) == 39 and S3 < S4 and max(S4) == 283
    for exps in set(permutations((7, 5, 3, 2))):
        assert many_prime_family(zip((2, 3, 5, 7), exps))
    for exps in set(permutations((4, 3, 3, 3, 2))):
        assert many_prime_family(zip((2, 3, 5, 7, 11), exps))
    assert many_prime_family([(2, 100), (19, 100), (43, 3), (73, 3)])
    assert many_prime_family([(2, 100), (3, 100), (43, 3), (73, 2)])
    assert many_prime_family([(2, 100), (3, 100), (43, 3), (73, 3)]) is None
    assert many_prime_family([(2, 4), (3, 3), (5, 3), (7, 3), (11, 3)]) is None
    assert many_prime_family([(2, 3), (3, 3), (5, 4), (7, 3), (11, 3)])
    assert many_prime_family([(2, 3), (3, 2), (5, 3), (7, 5), (11, 3)]) is None
    expected = 7**5 * (2*3*5*11)**2
    visited = []
    def prior(n):
        visited.append(n)
        return n == expected
    assert many_prime_family([(2, 3), (3, 2), (5, 3), (7, 5), (11, 3)], prior)
    assert visited == [expected]


if __name__ == '__main__':
    assert odd_residue_periods(23, 13, 449, 463) == (168, (42, 46, 126, 130))
    assert odd_residue_periods(23, 67, 449, 463) == (21, ())
    assert odd_residue_periods(227, 3, 113, 127) == (336, (44, 100))
    test_quadratic_family()
    test_catalogue()
    print('four/five-prime catalogue checks passed')
