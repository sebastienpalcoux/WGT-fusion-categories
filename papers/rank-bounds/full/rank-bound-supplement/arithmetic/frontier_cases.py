"""The high-four omission signatures at ranks 41–43, before the final exact root test."""
from argparse import ArgumentParser
import json

from large_four import MIXED, PRIME, preliminary, profiles
from phase_deficits import pure_check


def cases(rank):
    if rank not in (41,42,43):
        raise ValueError('the prime-support estimates cover 41 through 43')
    for rows, five, support, modes in preliminary(rank):
        # A prime >=23, the two smallest other odd primes, nine fours, the
        # vacuum and three odd hits consume at least 23+3+5+9+1+3 seats.
        if next(profiles(rank,rows,five,support,modes),None) is None:
            continue
        roots = None
        if MIXED not in modes:
            neutral = (6 if five and modes.count(PRIME)>1 else 4 if five else
                       3 if PRIME in modes else 2)
            roots = pure_check(rank,rows,support,modes,neutral)
            if not roots:
                continue
        yield dict(rows=rows, five=five, support=support, modes=modes, roots=roots)


def verify():
    """Regenerate all omission tables used by the native search."""
    from pathlib import Path
    import re
    engine = Path(__file__).resolve().parents[1] / "census/engine_src/square44.cpp"
    source = engine.read_text()
    primes = (2, 3, 5, 7, 11, 13, 17, 19)
    for rank, expected in ((41, 16), (42, 21), (43, 30)):
        name = "four_cases" if rank == 41 else f"four_cases{rank}"
        table = source.split(name + "[]={", 1)[1].split("};", 1)[0]
        native = {tuple(map(int, row)) for row in
                  re.findall(r"\{(\d+),(\d+),(\d+),(\d+)\}", table)}
        encoded, roots = set(), set()
        mixed = 0
        for case in cases(rank):
            values = case["roots"]
            if values is None:
                mixed += 1
                values = (0,)
                if case["modes"][0] != MIXED or MIXED in case["modes"][1:]:
                    raise ArithmeticError(("mixed support", case))
            else:
                roots.update(values)
            support = 1 + sum(1 << primes.index(p) for p in case["support"])
            prime = sum(1 << primes.index(p) for p, mode in
                        zip(case["support"], case["modes"]) if mode == PRIME)
            encoded.update((case["rows"], support, prime, value) for value in values)
        expected_roots = {15482880} if rank == 41 else {136080, 15482880}
        if mixed != expected or roots != expected_roots or encoded != native:
            raise ArithmeticError((rank, mixed, roots, encoded ^ native))
        print(f"Rank {rank}: {mixed} mixed signatures and {len(roots)} roots; native table agrees.")


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('rank', nargs='?', type=int, choices=(41, 42, 43))
    args = parser.parse_args()
    if args.rank is None:
        verify()
    else:
        print(json.dumps(list(cases(args.rank)), indent=2))
