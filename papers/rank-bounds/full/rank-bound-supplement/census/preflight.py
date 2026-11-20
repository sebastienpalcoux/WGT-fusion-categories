"""Recompute the structural prerequisites through rank44; no old receipts."""
from fractions import Fraction
import json
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'arithmetic'))


def require(test, message):
    if not test:
        raise RuntimeError(message)


def _signature(item):
    return (item['rows'], item['five'], tuple(item['support']), tuple(item['modes']))


def _encoded_case(item, root=0):
    import large_four
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23)
    support = 1 + sum(1 << primes.index(p) for p in item['support'])
    primary = sum(1 << primes.index(p) for p, mode
                  in zip(item['support'], item['modes']) if mode == large_four.PRIME)
    return item['rows'], support, primary, root


def _native_cases(source, name):
    marker = name + '[]={'
    require(source.count(marker) == 1, 'Missing or duplicate native structural table: '+name)
    body = source.split(marker, 1)[1].split('};', 1)[0]
    rows = [tuple(map(int, match)) for match in
            re.findall(r'\{(\d+),(\d+),(\d+),(\d+)\}', body)]
    require(rows and len(rows) == len(set(rows)), 'Empty or duplicate native structural table: '+name)
    return set(rows)


def check_lower_rank_envelope(mixed, removed, pure, primary, progress=None):
    """Fresh mixed tables at40--43, covered by the freshly solved44 envelope.

The omission bounds are monotone in their rank budget: append ordinary
zero-mask slots to an earlier omission profile. The pure deficit search
allows at most its slot budget, and the root closer allows zero padding.
Thus the fresh rank44 primary calculation covers every smaller rank.
Only the sharper lower-rank mixed restrictions need separate regeneration.
"""
    import large_four
    source = (ROOT/'engine_src/square44.cpp').read_text()
    envelope = {_signature(item) for item in mixed+removed+pure}
    removed_signatures = {_signature(item) for item in removed}
    primary_signatures = {_signature(item) for item in pure}
    fresh_roots = {_encoded_case(entry['profile'], root)
                   for entry in primary for root in entry['roots']}
    native44 = _native_cases(source, 'four_cases44')
    require({case for case in native44 if case[3] == 0}
            == {_encoded_case(item) for item in mixed},
            'Fresh rank44 mixed signatures disagree with the native table.')
    require({case for case in native44 if case[3] != 0} == fresh_roots,
            'Fresh rank44 roots disagree with the native structural table.')
    records = []
    for rank in (40, 41, 42, 43):
        lower_mixed, lower_pure, lower_removed = [], [], []
        for rows, five, support, modes in large_four.preliminary(rank):
            if next(large_four.profiles(rank, rows, five, support, modes), None) is None:
                continue
            item = dict(rows=rows, five=five, support=list(support), modes=list(modes))
            signature = _signature(item)
            require(signature in envelope, 'Lower-rank structural case escaped the rank44 envelope.')
            if large_four.MIXED not in modes:
                require(signature in primary_signatures,
                        'Lower-rank primary case is not covered by the fresh root calculation.')
                lower_pure.append(item)
            elif signature in removed_signatures:
                lower_removed.append(item)
            else:
                require(modes[0] == large_four.MIXED
                        and large_four.MIXED not in modes[1:],
                        'Unexpected lower-rank mixed-prime case.')
                lower_mixed.append(item)
        native_roots = []
        if rank == 40:
            require(len(lower_mixed) == 10,
                    'Unexpected number of mixed signatures through rank40.')
            for item in lower_mixed:
                rows, five, support, modes = _signature(item)
                seven = modes[support.index(7)] == large_four.PRIME if 7 in support else False
                require(rows <= 12 and support == (3, 5, 7)
                        and (rows > 11 or five or seven)
                        and (not (rows == 12 or rows >= 11 and five) or not seven),
                        'A rank40 signature violates the native through40 reduction.')
        else:
            name = 'four_cases' if rank == 41 else 'four_cases'+str(rank)
            native = _native_cases(source, name)
            require({case for case in native if case[3] == 0}
                    == {_encoded_case(item) for item in lower_mixed},
                    f'Fresh rank{rank} mixed signatures disagree with the native table.')
            positive = {case for case in native if case[3] != 0}
            require(positive <= fresh_roots,
                    'A discarded lower-rank primary root lacks a fresh rank44 closure.')
            native_roots = sorted(case[3] for case in positive)
        records.append(dict(rank=rank, mixed_signatures=lower_mixed,
                            primary_profiles_in_rank44_envelope=lower_pure,
                            mixed_cases_closed_by_rank44_bounds=lower_removed,
                            native_positive_roots_covered=native_roots,
                            native_agreement=True))
        if progress is not None:
            progress(f'Structural checks: rank{rank} has {len(lower_mixed)} mixed signatures; native agreement verified.', flush=True)
    return dict(status='complete', fresh=True, ranks=list(range(1, 45)),
                mixed_checks=records,
                primary_coverage='At-most44 slot budget; final divisor sums permit zero padding.',
                lower_omission_coverage='Rank40 inequalities cover smaller ranks by adding zero-mask ordinary slots.')


def check_two_five(signatures):
    """Recompute the six omission profiles and their exact closure bounds.

    The six omission profiles and their closure inequalities are
    regenerated directly. No prior calculation or completion JSON is read.
    The structural implications of these inequalities are the manuscript's
    rank-44 mixed two-five lemma; the arithmetic is regenerated here.
    """
    import large_four

    expected = {
        (10, False, (3, 5, 7), (0, 1, 2)): 5,
        (11, False, (3, 5, 7), (0, 1, 2)): 1,
    }
    actual = [(s['rows'], s['five'], tuple(s['support']), tuple(s['modes']))
              for s in signatures]
    require(len(actual) == len(expected) and set(actual) == set(expected),
            'Incomplete or unexpected mixed two-five signature set.')
    F = Fraction
    maxima, profiles = {}, []
    for rows, five, support, modes in sorted(actual):
        choices = list(large_four.profiles(44, rows, five, support, modes))
        require(len(choices) == expected[rows, five, support, modes],
                'Unexpected mixed two-five omission profile count.')
        maxima[rows] = max(item[3] for item in choices)
        for odd, even, other, mass in choices:
            require(odd == (1, 2, 4) and 5 in even and 6 in even,
                    'Mixed two-five witness shape changed.')
            if rows == 11:
                require(even == (3, 5, 6) and not other,
                        'Eleven-four mixed two-five shape changed.')
            profiles.append(dict(rows=rows, odd=list(odd), even=list(even),
                                 other=list(other), mass=str(mass)))
        excess = maxima[rows] - 1
        require(excess < 2 * (F(1, 36) - F(1, 100)),
                'Additional five-active pair is not excluded.')
        require(excess < 4 * (F(1, 100) - F(1, 400)),
                'Higher dyadic height in the five-fibers is not excluded.')
        require(excess < F(3, 81) - F(3, 729),
                'Higher pure ternary level is not excluded.')
    require(maxima == {10: F(77237, 75600), 11: F(233293, 226800)},
            'Mixed two-five upper masses changed.')
    require(maxima[10] - 1 < 2 * (F(1, 36) - F(1, 81)),
            'Two additional ternary labels are not excluded.')
    require(maxima[10] - 1 < F(1, 36) - F(1, 729),
            'An additional higher-level ternary label is not excluded.')

    def base(rows, nines):
        return F(rows, 16) + F(16, 100) + F(nines, 81) + F(1, 7)

    require(base(11, 3) > 1, 'Eleven-four two-five mass does not overshoot.')
    old_neutral = sum((F(1, b * b) for b in (72, 40, 56, 252, 140)),
                      F(1, 2520**2))
    require(old_neutral == F(13, 10800), 'Neutral witness mass changed.')
    four_nine_residual = 1 - base(10, 4)
    require(F(1, 36) > four_nine_residual > F(1, 64) + old_neutral,
            'Four-nine remainder is not excluded.')
    residual = 1 - base(10, 3)
    require(residual == F(1327, 37800), 'Three-nine remainder changed.')
    gap = residual - F(1, 36) - F(1, 144)
    require(gap == F(29, 75600), 'Required neutral gap changed.')
    large_five = sum((F(1, b * b) for b in (216, 200, 56, 756, 700)),
                     F(1, 37800**2))
    small_five = sum((F(1, b * b) for b in (864, 160, 224, 756, 140)),
                     F(1, 30240**2))
    require(large_five == F(43949, 119070000) < gap,
            'Neutral witnesses at five-height at least two are not excluded.')
    require(small_five == F(2873, 25401600) < gap,
            'Neutral witnesses at five-height one are not excluded.')
    require(F(2, 64) + max(large_five, small_five) < residual,
            'Remainder without denominator six is not excluded.')
    require(F(2, 36) > residual and F(1, 36) + F(1, 121) > residual,
            'Remainder with two small denominators is not excluded.')
    require([a for a in range(3, 9) if (135 * pow(2, a, 7)) % 7 in (1, 6)] == [5, 8],
            'Seven-packet two-adic congruence changed.')
    return dict(status='complete', fresh=True, omission_profiles=profiles,
                maximum_masses={str(m): str(v) for m, v in maxima.items()},
                three_nine_residual=str(residual),
                four_nine_residual=str(four_nine_residual),
                required_gap=str(gap), neutral_upper_c_ge2=str(large_five),
                neutral_upper_c1=str(small_five), surviving_modes=0)


def run(progress=print):
    """Fresh structural checks for every rank from1 through44."""
    if sys.flags.optimize:
        raise RuntimeError('Run without Python -O: structural proof assertions must remain enabled.')
    import large_four
    from int_deficits import Stats, pure_check
    import gauss_divisors

    started = time.monotonic()
    large_four.PRIMES = (3, 5, 7, 11, 13, 17, 19, 23)
    mixed, removed, pure = [], [], []
    for rows, five, support, modes in large_four.preliminary(44):
        if next(large_four.profiles(44, rows, five, support, modes), None) is None:
            continue
        item = dict(rows=rows, five=five, support=list(support), modes=list(modes))
        if large_four.MIXED in modes:
            if 5 in support and modes[support.index(5)] == large_four.MIXED:
                require((rows, five, support, modes) in (
                    (10, False, (3, 5, 7), (0, 1, 2)),
                    (11, False, (3, 5, 7), (0, 1, 2))),
                    'Unexpected mixed two-five profile; no exclusion is certified.')
                removed.append(item)
            else:
                require(modes[0] == large_four.MIXED and large_four.MIXED not in modes[1:],
                        'Unexpected mixed-prime signature.')
                mixed.append(item)
        else:
            old = (6 if five and modes.count(large_four.PRIME) > 1 else 4 if five
                   else 3 if large_four.PRIME in modes else 2)
            pure.append(dict(item, neutral=max(old, 3)))

    table = json.loads((ROOT/'arithmetic/high-four44-mixed.json').read_text())
    require(table['rank'] == 44 and mixed == table['signatures'],
            'Regenerated mixed signatures disagree with the native/exact input table.')
    require(len(mixed) == 39 and len(removed) == 2 and len(pure) == 46,
            'Unexpected structural profile count.')
    two_five = check_two_five(removed)
    progress('Structural checks: regenerated 39 mixed signatures and 46 primary profiles.', flush=True)
    primary, closures = [], []
    total_nodes = 0
    for index, profile in enumerate(pure):
        stats = Stats()  # No deadline or numerical search cap.
        roots = list(pure_check(44, profile['rows'], tuple(profile['support']),
                                tuple(profile['modes']), profile['neutral'], stats))
        total_nodes += stats.nodes
        primary.append(dict(index=index, profile=profile, roots=roots, **stats.info()))
        for root in roots:
            require(profile['rows'] == 9 and profile['support'] == [3, 5, 7]
                    and profile['modes'] == [0, 2, 2],
                    'A primary profile requires a new proof; no exclusion is certified.')
            require(Fraction(9, 16) + Fraction(1, 5) + Fraction(12, 49) > 1,
                    'Invalid prime-packet multiplicity inequality.')
            counts = gauss_divisors.verify_root(root, rank=44, fours=9)
            closures.append(dict(profile=index, root=root, counts=list(counts), status='excluded'))
        if (index+1) % 10 == 0 or index+1 == len(pure):
            progress(f'Structural checks: {index+1}/{len(pure)} primary profiles complete.', flush=True)
    lower = check_lower_rank_envelope(mixed, removed, pure, primary, progress)
    progress('Structural checks: fresh mixed reductions through40 and at41,42,43 agree with native rules.', flush=True)
    return dict(format=2, rank=44, ranks=list(range(1, 45)), status='complete', fresh=True,
                mixed_signatures=mixed, mixed_two_five_excluded=removed,
                mixed_two_five_check=two_five,
                primary_profiles=primary, root_closures=closures,
                lower_rank_structural_checks=lower,
                primary_nodes=total_nodes, all_primary_survivors=0,
                seconds=time.monotonic()-started,
                method='Regenerated all mixed reductions through44 and recomputed the entire44 primary envelope and divisor closures; no prior receipts read.')
