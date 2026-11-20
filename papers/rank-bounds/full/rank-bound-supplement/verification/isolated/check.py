"""Recheck all edgeless reservations and compare every outcome with the release."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

def main():
    here = Path(__file__).resolve().parent
    with TemporaryDirectory() as folder:
        actual = Path(folder)/'receipt.json'
        subprocess.run([sys.executable, '-B', str(here/'isolated209.py'), '--receipt', str(actual)], check=True)
        if json.loads(actual.read_text()) != json.loads((here/'isolated209.json').read_text()):
            raise ValueError('Per-reservation outcomes differ from release')

if __name__ == '__main__':
    main()
