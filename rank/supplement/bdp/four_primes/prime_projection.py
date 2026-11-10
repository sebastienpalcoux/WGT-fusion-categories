"""Project the periodic dyadic test to the odd prime with unbounded exponent."""
import argparse
import unbounded_four as periodic


def projection(height):
    totals = [0,0,0]
    primes = set()
    for case in periodic.configurations(height):
        totals[0] += 1
        budgets = periodic.common_budgets(case)
        if budgets is None:
            continue
        totals[1] += 1
        if periodic.has_common_gram(case,budgets):
            totals[2] += 1
            primes.add(case.q)
    return totals,sorted(primes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("heights",type=int,nargs="+")
    args = parser.parse_args()
    for height in args.heights:
        if height<2:
            parser.error("the dyadic root height must be at least two")
        totals, primes = projection(height)
        print(f"h={height}: {totals}; heavy-prime projection {primes}",flush=True)
