"""Reproduce examples/sample-output.txt.

Run from the examples/ directory:  python3 ../pwtp/run_example.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from load import load_csv
from pwtp import BenchmarkNormaliser, PWTPParams, residue_pwtp, target_pwtp


def main() -> None:
    conformers = load_csv("sample-input.csv", pocketless=12)
    norm = BenchmarkNormaliser().fit(conformers)
    params = PWTPParams()
    scores = residue_pwtp(conformers, norm, params)

    print("PWTP reference implementation - sample output")
    print("input: examples/sample-input.csv  (8 pocket-bearing conformers + 12 pocketless = 20 total)")
    print()
    print(f"{'residue':>8} {'occupancy':>10} {'quality':>9} {'PWTP':>8}")
    for s in sorted(scores.values(), key=lambda s: -s.pwtp):
        print(f"{s.residue:>8} {s.occupancy:>10.4f} {s.quality:>9.4f} {s.pwtp:>8.4f}")
    print()
    print("target-level PWTP (top-5 residues):", round(target_pwtp(scores, k=5), 4))
    print()
    print("Read occupancy and quality as a pair: residues 41-44 line a high-quality")
    print("pocket in a minority of conformers, residues 80-82 line a poor pocket more")
    print("often. The product alone does not distinguish those two situations.")


if __name__ == "__main__":
    main()
