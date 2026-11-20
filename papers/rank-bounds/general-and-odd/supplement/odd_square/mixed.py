"""Classify the reservations having a nonempty twist graph."""
from highprofile import possibilities as high_order_possible
from sparse_profile import sparse_possible
import dense_budget
import edge_budget
import joint_gauss
from profiles import witnesses
from quadratic import two_low_possible

BUDGET = 2000000
FIRST_PASS = 15000

def mixed_rows(profiles):
    reasons = {}
    steps = {}
    for index, row in enumerate(profiles):
        support, choices, edges, nines = row
        if not edges:
            continue
        if tuple(map(tuple, edges)) != ((3, 5),):
            raise ValueError("Unexpected component", index, edges)
        marked = {p: tuple(c) for p, c in zip(support[1:], choices) if c[0]}
        isolated = tuple(p for p in support if p not in (3, 5) and p not in marked)
        shapes = set()
        edge_budget.compatible(209, support, marked, edges, nines, witnesses, shapes)
        groups = {}
        dense_unresolved = False
        for shape in sorted(shapes):
            if len(marked) == 2 and len(support) == 4:
                if not two_low_possible(209, nines, marked, shape):
                    continue
            if shape[0] in (96, 128) and not high_order_possible(support, marked, nines, shape, witnesses):
                continue
            if any(dense_budget.possible(209, support, marked, nines, shape, q, witnesses) for q in isolated):
                dense_unresolved = True
                continue
            if not sparse_possible(row, shape, witnesses):
                continue
            groups.setdefault((shape[0], shape[3], shape[1] in (12, 18)), []).append(shape)
        if dense_unresolved:
            reasons[index] = "unresolved dense isolated prime"
            continue
        used = 0
        outcome = "excluded"
        for family in groups.values():
            def search(candidates, limit, seven=None):
                nonlocal used
                if not candidates:
                    return "excluded"
                joint_gauss.remaining_visits = limit
                envelope = tuple(map(min, zip(*candidates)))
                try:
                    survives = joint_gauss.possible(209, tuple(support), nines, envelope,
                                                    witnesses, marked, profile_shapes=candidates,
                                                    seven=seven)
                except RuntimeError as error:
                    if error.args != ("budget",):
                        raise
                    return "node limit"
                except (ValueError, ArithmeticError) as error:
                    if error.args[0] not in {"joint phase order is unresolved",
                                             "unresolved dense isolated prime",
                                             "unresolved pure neutral quota",
                                             "unresolved joint support",
                                             "unexpected prime-denominator class"}:
                        raise
                    return error.args[0]
                finally:
                    used += limit-max(0,joint_gauss.remaining_visits)
                return "retained by arithmetic" if survives else "excluded"

            can_split = isolated == (7,) and all(s[0] in (32,48,64) for s in family)
            limit = min(FIRST_PASS, BUDGET-used) if can_split else BUDGET-used
            outcome = search(family, limit)
            if outcome == "excluded":
                continue
            if outcome != "node limit" or not can_split:
                break
            for total in (12, 24, 36):
                for square in range(0, total+1, 12):
                    for cube in range(0,total-square+1,12):
                        spec = (total, square, cube)
                        candidates = [s for s in family if sparse_possible(row,s,witnesses,seven=spec)]
                        outcome = search(candidates,BUDGET-used,spec)
                        if outcome != "excluded":
                            break
                    if outcome != "excluded":
                        break
                if outcome != "excluded":
                    break
            if outcome != "excluded":
                break
        reasons[index] = outcome
        steps[index] = used
    return {i: (v,steps.get(i,0)) for i,v in reasons.items()}
