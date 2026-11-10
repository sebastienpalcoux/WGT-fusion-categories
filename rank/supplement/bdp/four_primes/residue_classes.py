"""Cyclic quotient examples for the two higher exponents in Section 14.1."""


def local_indices(t):
    """Return the odd indices H of t-1 allowed by the one-endpoint bound."""
    return [h for h in range(5, t, 2)
            if (t - 1) % h == 0 and 3 * ((t - 1) // h) <= t + 1]


def two_endpoints(p, q, u, v):
    def product(x, y):
        return (x[0] * y[0] % u, x[1] * y[1] % v)

    allowed = set()
    for hu in local_indices(u):
        for hv in local_indices(v):
            ju, jv = (u - 1) // hu, (v - 1) // hv
            if ju + 3 <= hv or jv + 3 <= hu:
                for sign in (1, -1):
                    allowed.add((sign * hu * pow(v, -1, u) % u,
                                 sign * hv * pow(u, -1, v) % v))

    subgroup = set()
    x = (1, 1)
    while x not in subgroup:
        subgroup.add(x)
        x = product(x, (p, p))

    # E and the coset q^k <p> meet exactly when e/q^k lies in <p> for some e.
    inv_q = (pow(q, -1, u), pow(q, -1, v))
    inv_power = (1, 1)
    residues = []
    k = 0
    while True:
        if any(product(e, inv_power) in subgroup for e in allowed):
            residues.append(k)
        k += 1
        inv_power = product(inv_power, inv_q)
        if inv_power in subgroup:
            return k, residues


if __name__ == "__main__":
    for q in (13, 67):
        order, residues = two_endpoints(23, q, 449, 463)
        print(f"Q={q}: order of Q modulo <P> {order}, residues {residues}")
        if q == 13:
            assert (order, residues) == (168, [42, 46, 126, 130])
        else:
            assert not residues
