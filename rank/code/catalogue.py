"""Sufficient divisor-cone certificates for the dimension of the category.

Dependent/local tests requiring twist, rank or block information are stated
separately in the manuscript. A negative return always retains the family.
"""
from functools import lru_cache
from families import factors, dimension_family, family
from many_prime_families import many_prime_family
from three import covered_three
from support_certificates import extension_certificate


def direct(items):
    legacy = family(items)
    if legacy:
        return legacy
    if len(items) == 3:
        return covered_three(*zip(*items))
    return many_prime_family(items, covered=certificate)


@lru_cache(maxsize=65536)
def certificate(D):
    old = dimension_family(D)
    if old:
        return old
    profile = tuple((p, e) for p, e in factors(D) if e >= 2)
    proof = extension_certificate(profile, direct)
    return proof["rule"] if proof else None
