#!/usr/bin/env python3
"""Reproduce the seven trace rows for total height at two higher primes."""

from check_four import axis_floor, axis_options, axis_trace_costs, common_terms, vacuum_cases


def excluded_rows(h, k, case):
    q,u,v,ju,jv,sign = case
    ends = ((u,v,ju), (v,u,jv))
    budgets = [j*(t+1-j) for t,w,j in ends]
    odd_axis = [axis_floor(q,k,2,t,w,j,sign,1,False) for t,w,j in ends]
    two_axis = [axis_floor(2,h,q,t,w,j,sign,1,False) for t,w,j in ends]
    caps = [budget-a-b for budget,a,b in zip(budgets,odd_axis,two_axis)]
    for i,l,zu,zv,m in common_terms(2,h,k,case,caps):
        if zu*zv >= 0:
            continue
        for (t,w,j),budget,odd_cost,z in zip(ends,budgets,odd_axis,(zu,zv)):
            available = budget-odd_cost-m*z*z
            options = axis_options(2,h,q,t,w,j,sign,1)
            residues = tuple(sorted(axis_trace_costs(options,q,available,1)))
            needed = sign*j % q
            if needed not in residues:
                yield (sign*q,u,v,ju,jv,t,available,residues,needed), (i,l)
                break
        else:
            raise ArithmeticError(f"unexcluded common degree at {(h,k,case,i,l)}")


def main():
    heights = ((2,2),(2,3),(3,2),(2,4),(4,2))
    totals = []
    rows = {}
    for h,k in heights:
        cases = list(vacuum_cases(2,h,k))
        totals.append(len(cases))
        for case in cases:
            for key,powers in excluded_rows(h,k,case):
                rows.setdefault((h,k,*key),[]).append(powers)
    print("Vacuum counts:", *totals)
    print("h k signed-q u v j_u j_v powers t K residues required")
    for (h,k,sq,u,v,ju,jv,t,cap,res,need),powers in sorted(
        rows.items(), key=lambda item: (item[0][0:2], item[0][3], item[0][2])
    ):
        print(h,k,sq,u,v,ju,jv,powers,t,cap,res,need)
    if totals != [12,14,92,18,365] or len(rows) != 7:
        raise ArithmeticError("counts differ from the printed theorem")


if __name__ == "__main__":
    main()
