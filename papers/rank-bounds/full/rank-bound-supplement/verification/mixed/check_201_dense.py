"""Check the historic dense boundary with the shared mixed phase intervals."""
from fractions import Fraction as F
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
root = HERE.parents[1]
sys.path[:0] = [str(HERE), str(root / "arithmetic")]
from dense_screen import dense_possible
from odd_square import witnesses


def check():
    support = (3,5,7,11)
    marked = {11: (1,22,20)}
    tested = []
    for nines in (34,36,38,40):
        for fives in range(40, 81-nines, 4):
            # 201-(1+22+32+12+42+12) = 80.
            if (nines,fives) == (36,44):
                upper = (F(36,81) + F(22,121) + F(76,225) + F(12,45**2)
                         + F(42,49**2) + F(2,99**2) + F(2,55**2)
                         + F(2,3773**2) + F(6,343**2) + F(1,105**2))
                assert upper == F(4070637023,4118135175) < 1
            else:
                assert not dense_possible(209, support, marked, nines,
                                          (32,12,fives,0,12,0), 7, witnesses)
            tested.append((nines,fives))
    assert tested == [(34,40),(34,44),(36,40),(36,44),(38,40),(40,40)]
    return tested


if __name__ == "__main__":
    print("Excluded dense boundary populations:", check())
