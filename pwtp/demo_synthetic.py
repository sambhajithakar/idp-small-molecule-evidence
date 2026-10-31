"""Synthetic sanity demonstration of the PWTP metric.

WHAT THIS IS: a check that the metric behaves as designed on inputs whose
ground truth we set ourselves. It belongs in Methods as a metric-validation
figure.

WHAT THIS IS NOT: a scientific result. No real protein, no real ensemble, no
real pocket detection. Nothing here says anything about whether any actual
disordered protein is druggable. Do not put these numbers in a Results
section.

Run from this directory:  python3 demo_synthetic.py
"""

from __future__ import annotations

import random

from pwtp import (
    BenchmarkNormaliser,
    Conformer,
    Pocket,
    PWTPParams,
    convergence_profile,
    residue_pwtp,
    target_pwtp,
)

SEED = 20261012
N_CONF = 4000
SITE = frozenset({41, 42, 43, 44, 45})  # the "binding site" residues
DECOY = frozenset({80, 81, 82})


def _pocket(rng: random.Random, residues: frozenset[int], scale: float) -> Pocket:
    """A pocket whose quality scales with `scale` in [0,1]."""
    return Pocket(
        volume=rng.gauss(120 * scale + 20, 10),
        hydrophobicity=max(0.01, min(1.0, rng.gauss(0.85 * scale + 0.05, 0.05))),
        druggability=max(0.01, min(1.0, rng.gauss(0.85 * scale + 0.05, 0.05))),
        residues=residues,
    )


def make_target(kind: str, rng: random.Random) -> list[Conformer]:
    """Three archetypes, matching the mechanism classes in the curation.

    pocket_like   a good pocket at the site in a minority of conformers
    remodelling   many poor, scattered pockets, nothing persistent at the site
    oligomer_like no single-chain pocket at all; engagement is intermolecular,
                  which a single-chain ensemble cannot represent
    """
    out: list[Conformer] = []
    for _ in range(N_CONF):
        pockets: list[Pocket] = []
        if kind == "pocket_like":
            if rng.random() < 0.08:                 # 8% occupancy, high quality
                pockets.append(_pocket(rng, SITE, 0.95))
            if rng.random() < 0.10:
                pockets.append(_pocket(rng, DECOY, 0.25))
        elif kind == "remodelling":
            if rng.random() < 0.45:                 # frequent but poor
                pockets.append(_pocket(rng, SITE, 0.18))
            if rng.random() < 0.45:
                pockets.append(_pocket(rng, DECOY, 0.18))
        elif kind == "oligomer_like":
            if rng.random() < 0.06:                 # only weak decoys
                pockets.append(_pocket(rng, DECOY, 0.20))
        else:
            raise ValueError(kind)
        out.append(Conformer(tuple(pockets)))
    return out


def main() -> None:
    rng = random.Random(SEED)
    targets = {k: make_target(k, rng) for k in ("pocket_like", "remodelling", "oligomer_like")}

    # Fit the normaliser across ALL targets, as the real pipeline must.
    norm = BenchmarkNormaliser().fit(c for ens in targets.values() for c in ens)
    print(f"normaliser ranges (fitted across all {len(targets)} targets):")
    for field, (lo, hi) in norm.ranges.items():
        print(f"  {field:<15} {lo:8.3f} .. {hi:8.3f}")

    k = len(SITE)
    print(f"\nper-target PWTP (top-{k} residues), and the site residues' decomposition:\n")
    header = f"{'target':<15}{'PWTP':>8}{'site occ':>10}{'site qual':>11}{'site PWTP':>11}"
    print(header)
    print("-" * len(header))

    for name, ens in targets.items():
        scores = residue_pwtp(ens, norm)
        tp = target_pwtp(scores, k=k)
        site = [scores[r] for r in sorted(SITE) if r in scores]
        if site:
            occ = sum(s.occupancy for s in site) / len(site)
            qual = sum(s.quality for s in site) / len(site)
            sp = sum(s.pwtp for s in site) / len(site)
            print(f"{name:<15}{tp:8.4f}{occ:10.3f}{qual:11.3f}{sp:11.4f}")
        else:
            print(f"{name:<15}{tp:8.4f}{'—':>10}{'—':>11}{'—':>11}  (no pocket ever at the site)")

    # The decomposition is the metric's reason to exist: show that a product
    # alone cannot distinguish rare-and-excellent from common-and-mediocre.
    print("\nwhy the pair (occupancy, quality) must be reported, not just the product:")
    rare_good = [
        Conformer((_pocket(rng, SITE, 1.0),)) if i < 100 else Conformer(())
        for i in range(1000)
    ]
    # scale 0.58 chosen by a scan so the two PWTP values land close together;
    # that near-equality is the point being demonstrated.
    common_mid = [
        Conformer((_pocket(rng, SITE, 0.58),)) if i < 460 else Conformer(())
        for i in range(1000)
    ]
    for label, ens in (("rare + excellent", rare_good), ("common + mediocre", common_mid)):
        s = residue_pwtp(ens, norm)[min(SITE)]
        print(f"  {label:<20} PWTP={s.pwtp:.4f}  occupancy={s.occupancy:.3f}  quality={s.quality:.3f}")
    print("  -> close PWTP, opposite practical meaning. A single scalar hides this.")

    # Convergence: occupancy is an ensemble property, so the estimate must plateau.
    print(f"\nconvergence of target PWTP for 'pocket_like' (k={k}):")
    prof = convergence_profile(targets["pocket_like"], norm, k=k, steps=8)
    final = prof[-1][1]
    for n, val in prof:
        drift = abs(val - final) / final if final else float("nan")
        print(f"  n={n:>5}  PWTP={val:.4f}  |rel. drift from full|={drift:6.2%}")

    # Exponent sensitivity, which the manuscript promises in Supporting Information.
    print("\nexponent sensitivity for 'pocket_like' (alpha,beta,gamma):")
    for params in (PWTPParams(1, 1, 1), PWTPParams(1, 0, 0), PWTPParams(0, 0, 1), PWTPParams(2, 1, 1)):
        s = residue_pwtp(targets["pocket_like"], norm, params)
        print(
            f"  ({params.alpha:g},{params.beta:g},{params.gamma:g})"
            f"  target PWTP={target_pwtp(s, k=k):.4f}"
        )

    print(
        "\nREMINDER: synthetic input, ground truth set by hand. "
        "Methods-validation only — not a Results number."
    )


if __name__ == "__main__":
    main()
