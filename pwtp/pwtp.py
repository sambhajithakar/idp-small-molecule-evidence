"""Persistence-weighted transient pocket (PWTP) score.

Reference implementation for the IDP-PocketBench manuscript.

The score collapses a conformational ensemble into a per-residue number that
separates two things a single druggability verdict confuses:

    PWTP(r) = occupancy(r) * conditional_quality(r)

where occupancy is the ensemble-weighted fraction of conformers in which the
residue lines a pocket, and conditional quality is the mean normalised pocket
score over just those conformers.

Report the pair, not only the product. A residue lining an excellent pocket in
2% of conformers and one lining a mediocre pocket in 40% can share a PWTP value
while meaning entirely different things for a drug-discovery decision.

Normalisation of volume / hydrophobicity / druggability must be computed across
the WHOLE benchmark, not per target, or targets stop being comparable. Use
`BenchmarkNormaliser` for that, and fit it once on all conformers of all
targets before scoring any of them.

No dependencies beyond the standard library and (optionally) numpy.
Pure-python by design so it is easy to audit, since the metric is a claimed
contribution of the paper.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

__all__ = [
    "Pocket",
    "Conformer",
    "BenchmarkNormaliser",
    "PWTPParams",
    "ResidueScore",
    "residue_pwtp",
    "target_pwtp",
    "convergence_profile",
]


# --------------------------------------------------------------------------
# data model
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Pocket:
    """One pocket detected in one conformer.

    volume: pocket volume in A^3, as reported by the pocket detector.
    hydrophobicity: hydrophobic-enclosure measure; detector-specific scale.
    druggability: detector's druggability score; detector-specific scale.
    residues: residue indices lining the pocket (1-based, matching the sequence).

    Raw detector scales are fine here. Normalisation happens later and across
    the whole benchmark.
    """

    volume: float
    hydrophobicity: float
    druggability: float
    residues: frozenset[int]

    def __post_init__(self) -> None:
        if self.volume < 0:
            raise ValueError("volume must be non-negative")
        if not self.residues:
            raise ValueError("a pocket must line at least one residue")


@dataclass(frozen=True)
class Conformer:
    """One ensemble member that passed the physical-plausibility filter.

    weight: statistical weight. Uniform for unbiased sampling; the reweighted
        value for biased sampling. Weights are normalised per ensemble by the
        scoring functions, so they need not sum to one on input.
    """

    pockets: tuple[Pocket, ...]
    weight: float = 1.0

    def __post_init__(self) -> None:
        if self.weight < 0:
            raise ValueError("conformer weight must be non-negative")


@dataclass
class PWTPParams:
    """Exponents emphasising each pocket property.

    Equal exponents are the default and the value reported in the main text.
    Any departure must be justified on a calibration subset held apart from the
    evaluation targets, with a sensitivity analysis.
    """

    alpha: float = 1.0  # volume
    beta: float = 1.0  # hydrophobicity
    gamma: float = 1.0  # druggability

    def __post_init__(self) -> None:
        for name in ("alpha", "beta", "gamma"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")


@dataclass
class ResidueScore:
    """PWTP decomposed. Always report occupancy and quality alongside pwtp."""

    residue: int
    pwtp: float
    occupancy: float
    quality: float


# --------------------------------------------------------------------------
# normalisation
# --------------------------------------------------------------------------


class BenchmarkNormaliser:
    """Min-max normaliser fitted across the whole benchmark.

    Fit once on every conformer of every target, then reuse the fitted
    instance for all scoring. Fitting per target would make targets
    incomparable, which defeats the purpose of the metric.

    Degenerate ranges (all values equal) map to 1.0 rather than dividing by
    zero: if a property does not vary across the benchmark it carries no
    discriminative information, and giving it a constant value of 1 makes it
    drop out of the product cleanly.
    """

    _FIELDS = ("volume", "hydrophobicity", "druggability")

    def __init__(self) -> None:
        self._lo: dict[str, float] = {}
        self._hi: dict[str, float] = {}
        self._fitted = False

    def fit(self, conformers: Iterable[Conformer]) -> "BenchmarkNormaliser":
        lo: dict[str, float] = {f: math.inf for f in self._FIELDS}
        hi: dict[str, float] = {f: -math.inf for f in self._FIELDS}
        seen = False
        for conf in conformers:
            for pocket in conf.pockets:
                seen = True
                for f in self._FIELDS:
                    v = float(getattr(pocket, f))
                    if v < lo[f]:
                        lo[f] = v
                    if v > hi[f]:
                        hi[f] = v
        if not seen:
            raise ValueError("cannot fit a normaliser with no pockets")
        self._lo, self._hi, self._fitted = lo, hi, True
        return self

    def __call__(self, pocket: Pocket) -> tuple[float, float, float]:
        if not self._fitted:
            raise RuntimeError("fit the normaliser on the benchmark first")
        out = []
        for f in self._FIELDS:
            lo, hi = self._lo[f], self._hi[f]
            v = float(getattr(pocket, f))
            if hi <= lo:
                out.append(1.0)
            else:
                # clamp so a value outside the fitted range (e.g. a target
                # added after fitting) cannot push a score above 1.
                out.append(min(1.0, max(0.0, (v - lo) / (hi - lo))))
        return tuple(out)  # type: ignore[return-value]

    @property
    def ranges(self) -> dict[str, tuple[float, float]]:
        if not self._fitted:
            raise RuntimeError("normaliser not fitted")
        return {f: (self._lo[f], self._hi[f]) for f in self._FIELDS}


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------


def _pocket_quality(
    pocket: Pocket,
    normaliser: BenchmarkNormaliser,
    params: PWTPParams,
) -> float:
    v, h, d = normaliser(pocket)
    return (v ** params.alpha) * (h ** params.beta) * (d ** params.gamma)


def residue_pwtp(
    ensemble: Sequence[Conformer],
    normaliser: BenchmarkNormaliser,
    params: PWTPParams | None = None,
) -> dict[int, ResidueScore]:
    """Per-residue PWTP for one target's ensemble.

    Where several pockets in the same conformer line the same residue, the
    highest-quality one is used: the residue's engagement potential in that
    conformer is set by the best pocket available to it, not by an average
    over pockets that a ligand would never occupy simultaneously.

    Residues never lining any pocket are absent from the result, not zero —
    the caller knows the sequence length and can fill zeros where it wants
    them. Returning only observed residues keeps the sequence-indexing
    responsibility with the caller rather than guessing it here.
    """
    params = params or PWTPParams()
    if not ensemble:
        raise ValueError("empty ensemble")

    total_weight = sum(c.weight for c in ensemble)
    if total_weight <= 0:
        raise ValueError("ensemble weights sum to zero")

    weighted_sum: dict[int, float] = defaultdict(float)
    occupied_weight: dict[int, float] = defaultdict(float)

    for conf in ensemble:
        w = conf.weight / total_weight
        if w == 0:
            continue
        best: dict[int, float] = {}
        for pocket in conf.pockets:
            q = _pocket_quality(pocket, normaliser, params)
            for r in pocket.residues:
                if q > best.get(r, -1.0):
                    best[r] = q
        for r, q in best.items():
            weighted_sum[r] += w * q
            occupied_weight[r] += w

    out: dict[int, ResidueScore] = {}
    for r, s in weighted_sum.items():
        occ = occupied_weight[r]
        # occ > 0 whenever r appears at all, so this division is safe.
        out[r] = ResidueScore(
            residue=r,
            pwtp=s,
            occupancy=occ,
            quality=s / occ,
        )
    return out


def target_pwtp(
    scores: Mapping[int, ResidueScore],
    k: int,
) -> float:
    """Target-level PWTP: mean over the top-k residues by PWTP.

    k is fixed across the benchmark so targets of different length stay
    comparable. If fewer than k residues line any pocket, the missing
    residues count as zero — a target with few pocket-lining residues should
    score low, and dropping to a smaller denominator would hide that.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    top = sorted((s.pwtp for s in scores.values()), reverse=True)[:k]
    return sum(top) / k


def convergence_profile(
    ensemble: Sequence[Conformer],
    normaliser: BenchmarkNormaliser,
    k: int,
    steps: int = 10,
    params: PWTPParams | None = None,
) -> list[tuple[int, float]]:
    """Target PWTP against a growing prefix of the ensemble.

    Occupancy is an ensemble property, so PWTP means nothing until it has
    plateaued. Plot this per target and exclude unconverged targets from
    headline numbers.

    Note this walks the ensemble in its given order. Shuffle beforehand, or
    pass trajectory frames in time order deliberately — a sorted ensemble
    would make the profile meaningless.
    """
    if steps <= 0:
        raise ValueError("steps must be positive")
    n = len(ensemble)
    if n == 0:
        raise ValueError("empty ensemble")

    profile: list[tuple[int, float]] = []
    for i in range(1, steps + 1):
        m = max(1, round(n * i / steps))
        prefix = ensemble[:m]
        try:
            scores = residue_pwtp(prefix, normaliser, params)
        except ValueError:
            # a prefix whose weights sum to zero carries no information
            continue
        profile.append((m, target_pwtp(scores, k)))
    return profile


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------


def _selftest() -> None:
    """Synthetic checks. Run from this directory: python3 pwtp.py

    These verify the arithmetic, not the science. The metric's scientific
    validity is what the benchmark in the manuscript is for.
    """
    P = Pocket

    # --- normalisation ---
    e_lo = Conformer((P(10.0, 0.0, 0.0, frozenset({1})),))
    e_hi = Conformer((P(110.0, 1.0, 1.0, frozenset({2})),))
    norm = BenchmarkNormaliser().fit([e_lo, e_hi])
    assert norm(P(10.0, 0.0, 0.0, frozenset({1}))) == (0.0, 0.0, 0.0)
    assert norm(P(110.0, 1.0, 1.0, frozenset({1}))) == (1.0, 1.0, 1.0)
    assert norm(P(60.0, 0.5, 0.5, frozenset({1}))) == (0.5, 0.5, 0.5)
    # out-of-range values clamp rather than exceeding 1
    assert norm(P(500.0, 9.0, 9.0, frozenset({1}))) == (1.0, 1.0, 1.0)

    # --- degenerate range maps to 1.0, so the property drops out ---
    flat = BenchmarkNormaliser().fit(
        [Conformer((P(50.0, 0.5, 0.5, frozenset({1})),))]
    )
    assert flat(P(50.0, 0.5, 0.5, frozenset({1}))) == (1.0, 1.0, 1.0)

    # --- occupancy: a perfect pocket in 1 of 4 conformers ---
    best = P(110.0, 1.0, 1.0, frozenset({7}))
    empty_like = P(10.0, 0.0, 0.0, frozenset({99}))
    ens = [Conformer((best,))] + [Conformer((empty_like,))] * 3
    scores = residue_pwtp(ens, norm)
    s7 = scores[7]
    assert abs(s7.occupancy - 0.25) < 1e-12, s7.occupancy
    assert abs(s7.quality - 1.0) < 1e-12, s7.quality
    assert abs(s7.pwtp - 0.25) < 1e-12, s7.pwtp

    # --- the decomposition is the whole point: equal pwtp, different meaning ---
    rare_good = [Conformer((P(110.0, 1.0, 1.0, frozenset({1})),))] + [
        Conformer((empty_like,))
    ] * 3
    a = residue_pwtp(rare_good, norm)[1]
    mid = P(60.0, 1.0, 1.0, frozenset({1}))  # quality 0.5 after normalisation
    common_mid = [Conformer((mid,))] * 2 + [Conformer((empty_like,))] * 2
    b = residue_pwtp(common_mid, norm)[1]
    assert abs(a.pwtp - b.pwtp) < 1e-12, (a.pwtp, b.pwtp)
    assert a.occupancy < b.occupancy and a.quality > b.quality

    # --- best pocket wins within a conformer, no averaging ---
    two = Conformer((P(110.0, 1.0, 1.0, frozenset({3})), P(10.0, 0.0, 0.0, frozenset({3}))))
    assert abs(residue_pwtp([two], norm)[3].pwtp - 1.0) < 1e-12

    # --- weights are honoured and normalised ---
    weighted = [
        Conformer((best,), weight=3.0),
        Conformer((empty_like,), weight=1.0),
    ]
    assert abs(residue_pwtp(weighted, norm)[7].occupancy - 0.75) < 1e-12
    # unnormalised weights give the same answer as normalised ones
    scaled = [
        Conformer((best,), weight=300.0),
        Conformer((empty_like,), weight=100.0),
    ]
    assert abs(
        residue_pwtp(weighted, norm)[7].pwtp - residue_pwtp(scaled, norm)[7].pwtp
    ) < 1e-12

    # --- zero-weight conformers are ignored, not counted ---
    with_zero = [
        Conformer((best,), weight=1.0),
        Conformer((empty_like,), weight=0.0),
    ]
    assert abs(residue_pwtp(with_zero, norm)[7].occupancy - 1.0) < 1e-12

    # --- exponents shift emphasis ---
    half = P(60.0, 0.5, 0.5, frozenset({5}))
    one = residue_pwtp([Conformer((half,))], norm, PWTPParams(1, 1, 1))[5].pwtp
    assert abs(one - 0.125) < 1e-12, one
    vol_only = residue_pwtp([Conformer((half,))], norm, PWTPParams(1, 0, 0))[5].pwtp
    assert abs(vol_only - 0.5) < 1e-12, vol_only

    # --- target-level score pads with zeros rather than shrinking k ---
    single = residue_pwtp([Conformer((best,))], norm)
    assert abs(target_pwtp(single, k=1) - 1.0) < 1e-12
    assert abs(target_pwtp(single, k=4) - 0.25) < 1e-12

    # --- convergence profile ---
    prof = convergence_profile(ens, norm, k=1, steps=4)
    assert [m for m, _ in prof] == [1, 2, 3, 4], prof
    assert abs(prof[-1][1] - 0.25) < 1e-12

    # --- input validation ---
    for bad in (
        lambda: P(-1.0, 0.0, 0.0, frozenset({1})),
        lambda: P(1.0, 0.0, 0.0, frozenset()),
        lambda: Conformer((best,), weight=-1.0),
        lambda: PWTPParams(-1, 1, 1),
        lambda: residue_pwtp([], norm),
        lambda: residue_pwtp([Conformer((best,), weight=0.0)], norm),
        lambda: target_pwtp(single, k=0),
        lambda: BenchmarkNormaliser().fit([]),
        lambda: BenchmarkNormaliser().fit([Conformer(())]),
    ):
        try:
            bad()
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError")

    try:
        BenchmarkNormaliser()(best)
    except RuntimeError:
        pass
    else:
        raise AssertionError("expected RuntimeError for unfitted normaliser")

    print("pwtp.py self-test: all checks passed")


if __name__ == "__main__":
    _selftest()
