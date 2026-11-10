"""Audit the newly generated exact-rank lists and compare the marginal yields of tests."""
import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from families import factors, dimension_family
from defect_mass import first_excess, joint_obstruction
from support_certificates import missing_pair
from catalogue import certificate
from principal_cover import dyadic_trace_possible, row_cover
from odd_blocks import odd_block_obstruction
from screen import screen
from orbit_screen import first_obstruction, small_codegree_obstruction


def validate(row, rank):
    if (not isinstance(row, list) or len(row) != rank or row != sorted(row)
            or any(type(x) is not int or x < 1 for x in row) or row.count(1) != 1):
        raise ValueError("expected positive sorted degrees of the requested rank and a single unit")
    D = sum(x*x for x in row)
    if any(D % (x*x) for x in row):
        raise ValueError("failed formal-codegree divisibility")
    return D


def character_verdict(row):
    first = screen(row)
    if first != 'open':
        return first
    f = small_codegree_obstruction(row)
    if f is not None:
        return 'codegree ' + str(f)
    orbit = first_obstruction(row)
    if orbit is not None:
        return 'orbit at ' + str(orbit[1])
    return 'open'


def fast_verdict(row, D=None):
    """First exclusion in the continuation; expensive row tests come last."""
    if D is None:
        D = sum(x*x for x in row)
    family = certificate(D)
    if family:
        return 'divisor family: ' + family
    cap = first_excess(row)
    if cap is not None:
        return f'defect-one f={cap}'
    block = odd_block_obstruction(row)
    if block is not None:
        return block
    profile = factors(D)
    prime = joint_obstruction(row, profile)
    if prime is not None:
        return f'joint mass at {prime}'
    pair = missing_pair(row, D, tuple(p for p, e in profile))
    if pair is not None:
        return 'principal pair'
    if D == 16*row[-1]**2:
        if row_cover(row) != 'passes' or not dyadic_trace_possible(row):
            return 'dyadic row at 16'
    return character_verdict(row)


def fast_report(rows, rank):
    counts = Counter()
    survivors = []
    for row in rows:
        reason = fast_verdict(row, validate(row, rank))
        counts[reason] += 1
        if reason == 'open':
            survivors.append(row)
    return {'rank':rank, 'profiles':sum(counts.values()), 'counts':dict(sorted(counts.items())),
            'survivors':survivors}


def audit(rows, rank):
    proof, strategies, joint, top = Counter(), Counter(), Counter(), Counter()
    qfb = Counter()
    open_rows = []
    for row in rows:
        D = validate(row, rank)
        fac = factors(D)
        q = 1
        for p, e in fac:
            if e % 2:
                q *= p
        basic = dimension_family(D)
        extended = certificate(D)
        cap = first_excess(row)
        mass = joint_obstruction(row, fac)
        pair = missing_pair(row, D, tuple(p for p, _ in fac))
        block = odd_block_obstruction(row)
        joint['odd block or slice'] += block is not None
        joint['catalogue including legacy'] += bool(extended)
        joint['legacy dimension families'] += bool(basic)
        joint['defect-one multiplicity'] += cap is not None
        joint['joint mass'] += mass is not None
        joint['principal pair'] += pair is not None
        after_basic = not basic
        f = D//row[-1]**2
        m = row.count(row[-1])
        if after_basic:
            form = ' '.join(f'{p}^{e}' for p, e in fac)
            top[form] += 1
            qfb[f'q={q}, f={f}, m={m}'] += 1
        if basic:
            verdict = 'legacy dimension families'
        elif cap is not None:
            verdict = f'defect-one f={cap}'
        else:
            verdict = character_verdict(row)
        if verdict == 'open' and block:
            verdict = block
        if verdict == 'open' and pair is not None:
            verdict = 'principal pair'
        proof[verdict] += 1
        if verdict == 'open':
            open_rows.append(row)
        if basic:
            reason = 'legacy dimension families'
        elif extended:
            reason = 'additional divisor family'
        elif cap is not None:
            reason = 'defect-one multiplicity'
        elif block is not None:
            reason = block
        elif mass is not None:
            reason = 'joint mass'
        elif pair is not None:
            reason = 'principal pair'
        else:
            reason = character_verdict(row)
        strategies[reason] += 1
    ordered = lambda c: dict(sorted(c.items(), key=lambda kv: (-kv[1],kv[0])))
    return {'rank':rank, 'profiles':len(rows), 'proof_first_exclusion':ordered(proof),
            'survivors':open_rows, 'family_strategy_first_exclusion':ordered(strategies),
            'overlapping_individual_tests':ordered(joint), 'dimensions_after_legacy':ordered(top),
            'roots_after_legacy':ordered(qfb)}


def read_rows(path, rank):
    rows = [json.loads(raw) for raw in path.read_bytes().splitlines()]
    if rows != sorted(rows) or len(rows) != len({tuple(row) for row in rows}):
        raise ValueError("duplicate or unordered candidate rows")
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('rank', type=int)
    ap.add_argument('input', type=Path)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = audit(read_rows(args.input,args.rank), args.rank)
    encoded = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({k:v for k,v in result.items() if k not in ('roots_after_legacy','dimensions_after_legacy')}))
