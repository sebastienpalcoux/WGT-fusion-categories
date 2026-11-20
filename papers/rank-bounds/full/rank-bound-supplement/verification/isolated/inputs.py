"""Regenerate the rank-209 input with the prime interval evaluated at rank 209."""
import json
from math import isqrt
from pathlib import Path
import sys

def main():
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root/'arithmetic'))
    import odd_square
    odd_square.RANK = 209
    odd_square.PRIMES = tuple(p for p in range(3, (209+1)//2+1, 2)
                              if all(p % q for q in range(3, isqrt(p)+1, 2)))
    count, rows = odd_square.run()
    fresh = json.dumps({'count':count,'rows':rows},indent=1).encode()
    if fresh != (root/'verification/neutral_orbits/rank209-input.json').read_bytes():
        raise ValueError('The regenerated rank-209 input differs from the release')
    print(f'{count} prime profiles, {len(rows)} reservations, identical to the saved input.')

if __name__ == '__main__':
    main()
