"""Independently execute the nonsquare census and validate all recorded counts."""
from pathlib import Path
import json
import subprocess
import sys

def main():
    here = Path(__file__).resolve().parent
    current = subprocess.run([sys.executable, '-E', '-B', str(here/'census.py')],
                             check=True, text=True, capture_output=True)
    result = json.loads(current.stdout)
    if result != json.loads((here/'receipt.json').read_text()):
        raise ValueError('Nonsquare census differs from the published receipt')
    counts = list(result['counts'].values())
    print(f"{len(counts)} supports, {sum(x[0] for x in counts)} initial profiles, "
          f"{len(result['height_survivors'])} final profiles and {len(result['solutions'])} completions.")

if __name__ == '__main__':
    main()
