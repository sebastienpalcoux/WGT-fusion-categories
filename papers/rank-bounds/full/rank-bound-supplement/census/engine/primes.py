"""Prime supports for the rank-parameterized denominator searches."""
from functools import lru_cache
from math import isqrt


@lru_cache(maxsize=8)
def primes_up_to(n):
    """Return the primes at most n in increasing order."""
    if n < 2:
        return ()
    sieve = bytearray(b'\x01') * (n+1)
    for p in range(2, isqrt(n)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((n-p*p)//p+1)
    return tuple(p for p in range(2, n+1) if sieve[p])
