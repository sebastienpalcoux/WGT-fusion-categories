"""Recalculate both stages of the rank-45 prefix cover."""
from pathlib import Path
import json

from verify_sector import verify
from sector_prefix import run

def main():
    here = Path(__file__).resolve().parent
    if verify() != json.loads((here/'sector-check.json').read_text()):
        raise ValueError('Independent arithmetic checks disagree')
    result = run(2)
    if result != json.loads((here/'sector-prefix45-depth2.json').read_text()):
        raise ValueError('The two-class prefix cover differs from its receipt')
    print(f"{result['starts']} starts, {result['plans']} sector families; {result['prefix_tests']} prefixes checked.")

if __name__ == '__main__':
    main()
