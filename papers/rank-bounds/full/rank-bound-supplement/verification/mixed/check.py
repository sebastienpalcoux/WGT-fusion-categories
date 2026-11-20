"""Replay the bounded mixed-component screen. Exact rank-209 phase and mass checks."""
import argparse
import hashlib
import json
import os
from collections import Counter
from functools import lru_cache
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BUDGET = 550000
FIRST_PASS = 15000
UNIFORM = "verification/odd_uniform/current-inertia-and-component-receipt.json"
SPARSE = "verification/sparse_charge/receipt.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculate(root, use_fiber_theorem=True):
    original = root / "arithmetic"
    sys.path[:0] = [str(HERE), str(original)]
    os.environ["ODD_NATIVE"] = "0"
    import dense_budget
    from dense_screen import dense_possible
    from highprofile import possibilities as high_order_possible
    from sparse_profile import sparse_possible
    import edge_budget
    from full_blocks import count
    import joint_gauss
    import nine_gram
    from odd_square import witnesses
    from quadratic import two_low_possible

    @lru_cache(None)
    def free_count(p, nines):
        # Deletion supplies two nonunit dual pairs at every prime.
        return count(p, nines) if p <= 13 else 4

    old_dense = dense_budget.possible

    def dense(*args):
        return old_dense(*args) and dense_possible(*args)

    dense_budget.free_count = nine_gram.count = joint_gauss.free_count = free_count
    dense_budget.possible = dense
    uniform = json.loads((root / UNIFORM).read_text())
    screen, = (s for s in uniform["screens"] if s["rank"] == 209)
    profiles = json.loads((root / screen["input_path"]).read_text())["rows"]
    sparse = json.loads((root / SPARSE).read_text())
    already_excluded = {item["index"] for item in sparse["rejected"]}
    reasons = {}
    steps = {}
    shapes_checked = 0
    one_pair_shapes = 0
    for index in screen["retained_indices"]:
        row = profiles[index]
        support, choices, edges, nines = row
        if not edges:
            continue
        if edges != [[3, 5]]:
            raise ValueError("Unexpected component", index, edges)
        if index in already_excluded:
            reasons[index] = "sparse top charge"
            continue
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
            # All 64 joint labels have order 15. A zero fourth entry fixes
            # their denominator at 15. The eight solely-five labels then
            # supply one neutral dual pair at each fifth root. The local
            # moments force denominator 15 there. The one-pair theorem
            # requires four dimension primes, contradiction on this input.
            if (use_fiber_theorem and len(support) == 3 and
                    shape[0] == 64 and shape[2] == 8 and shape[3] == 0):
                one_pair_shapes += 1
                continue
            if any(dense(209, support, marked, nines, shape, q, witnesses) for q in isolated):
                dense_unresolved = True
                continue
            if not sparse_possible(row, shape, witnesses):
                continue
            groups.setdefault((shape[0], shape[3], shape[1] in (12, 18)), []).append(shape)
        shapes_checked += sum(map(len, groups.values()))
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
    rejected = {"excluded", "sparse top charge"}
    # Existing baseline arithmetic files are bound by the uniform receipt.
    sources = [HERE / name for name in ("check.py", "edges.py", "joint_gauss.py", "dense64.py",
                                         "dense_phase.py", "dense_screen.py", "highprofile.py",
                                         "sparse_profile.py","profile_window.py","pure_window.py")]
    evidence = {path.name: digest(path) for path in sources}
    for path in [root / UNIFORM, root / SPARSE, root / screen["input_path"],
                 original / "edges.py", original / "joint_gauss.py"]:
        evidence[str(path.relative_to(root))] = digest(path)
    return {"rank": 209, "node_budget_per_profile": BUDGET,
            "input_and_source_sha256": evidence, "mixed_profiles": len(reasons), "recursive_nodes": sum(steps.values()),
            "max_nodes_on_a_profile": max(steps.values()),
            "prior_top_charge_indices": sorted(i for i,v in reasons.items() if v=="sparse top charge"),
            "shapes_after_preliminary_screens": shapes_checked,
            "one_neutral_pair_exclusions": one_pair_shapes,
            "outcomes": dict(sorted(Counter(reasons.values()).items())),
            "retained_indices": [i for i, reason in reasons.items() if reason not in rejected]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement", type=Path, default=HERE.parents[1])
    parser.add_argument("--write-receipt", action="store_true")
    parser.add_argument("--independent-replay", action="store_true",
                        help="also omit the one-neutral-pair theorem and compare to its separate receipt")
    options = parser.parse_args()
    result = calculate(options.supplement.resolve(), not options.independent_replay)
    text = json.dumps(result, indent=2) + "\n"
    target = HERE / ("independent.json" if options.independent_replay else "receipt.json")
    if options.write_receipt:
        target.write_text(text)
    elif not target.exists() or target.read_text() != text:
        raise SystemExit("Mixed screen receipt mismatch")
    print(f"{result['mixed_profiles']} mixed profiles, "
          f"{len(result['retained_indices'])} retained. {result['outcomes']}")


if __name__ == "__main__":
    main()
