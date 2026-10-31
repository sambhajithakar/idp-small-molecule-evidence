"""Load pocket-detection output into PWTP's data model.

Deliberately CSV-based rather than an fpocket/MDpocket parser: those output
formats differ by version and I could not verify one against a real run here,
so a guessed parser would be a silent source of error in a published metric.
Convert detector output to this schema yourself with a short script, check a
few rows by hand, and keep that converter in the repository.

Expected CSV schema, one row per (conformer, pocket):

    conformer_id    any string; rows sharing it belong to one conformer
    pocket_id       any string; unique within a conformer
    volume          float, pocket volume
    hydrophobicity  float, hydrophobic-enclosure measure
    druggability    float, detector druggability score
    residues        residue indices lining the pocket, space- or
                    semicolon-separated, 1-based
    weight          optional float, the conformer's statistical weight.
                    Must be identical on every row of a conformer. Omit the
                    column entirely for unbiased sampling.

Conformers with no detected pocket need no rows; pass their count to
`pocketless` so occupancy denominators stay correct. Getting this wrong
inflates occupancy, which is the single easiest way to overstate a target's
druggability.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable

from pwtp import Conformer, Pocket

__all__ = ["load_csv", "parse_residues"]

_SPLIT = re.compile(r"[;,\s]+")

_REQUIRED = {
    "conformer_id",
    "pocket_id",
    "volume",
    "hydrophobicity",
    "druggability",
    "residues",
}


def parse_residues(field: str) -> frozenset[int]:
    """Parse a residue list. Accepts space, comma or semicolon separators."""
    parts = [p for p in _SPLIT.split(field.strip()) if p]
    if not parts:
        raise ValueError("empty residue list")
    out = set()
    for p in parts:
        try:
            out.add(int(p))
        except ValueError as exc:
            raise ValueError(f"residue index {p!r} is not an integer") from exc
    return frozenset(out)


def load_csv(path: str | Path, pocketless: int = 0) -> list[Conformer]:
    """Read pocket rows into a list of Conformers.

    pocketless: how many plausibility-filtered conformers of this target had
        no pocket at all. They are appended as empty conformers so they
        contribute to the occupancy denominator. Pass the real number; a
        default of 0 is only correct when every conformer had a pocket.

    Row order is preserved so that `convergence_profile` sees frames in the
    order you supply them.
    """
    path = Path(path)
    grouped: dict[str, list[Pocket]] = {}
    weights: dict[str, float] = {}
    order: list[str] = []

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row")
        missing = _REQUIRED - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        has_weight = "weight" in reader.fieldnames

        for lineno, row in enumerate(reader, start=2):
            cid = (row["conformer_id"] or "").strip()
            if not cid:
                raise ValueError(f"{path}:{lineno} has an empty conformer_id")
            try:
                pocket = Pocket(
                    volume=float(row["volume"]),
                    hydrophobicity=float(row["hydrophobicity"]),
                    druggability=float(row["druggability"]),
                    residues=parse_residues(row["residues"]),
                )
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: {exc}") from exc

            if cid not in grouped:
                grouped[cid] = []
                order.append(cid)
            grouped[cid].append(pocket)

            if has_weight:
                raw = (row.get("weight") or "").strip()
                if raw:
                    w = float(raw)
                    prev = weights.get(cid)
                    if prev is not None and abs(prev - w) > 1e-12:
                        raise ValueError(
                            f"{path}:{lineno}: conformer {cid!r} has "
                            f"conflicting weights {prev} and {w}"
                        )
                    weights[cid] = w

    if not order and pocketless == 0:
        raise ValueError(f"{path} contains no rows")

    conformers = [
        Conformer(tuple(grouped[cid]), weight=weights.get(cid, 1.0))
        for cid in order
    ]
    if pocketless < 0:
        raise ValueError("pocketless must be non-negative")
    conformers.extend(Conformer((), weight=1.0) for _ in range(pocketless))
    return conformers


def _selftest() -> None:
    import tempfile

    assert parse_residues("1 2 3") == frozenset({1, 2, 3})
    assert parse_residues("4;5,6") == frozenset({4, 5, 6})
    assert parse_residues(" 7 ") == frozenset({7})

    rows = (
        "conformer_id,pocket_id,volume,hydrophobicity,druggability,residues,weight\n"
        "c1,p1,100,0.8,0.7,10 11 12,2.0\n"
        "c1,p2,50,0.4,0.3,20;21,2.0\n"
        "c2,p1,80,0.6,0.5,10 11,1.0\n"
    )
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "pockets.csv"
        p.write_text(rows)

        confs = load_csv(p, pocketless=2)
        assert len(confs) == 4, len(confs)
        assert len(confs[0].pockets) == 2
        assert confs[0].weight == 2.0
        assert confs[1].weight == 1.0
        assert confs[2].pockets == () and confs[3].pockets == ()
        assert confs[0].pockets[0].residues == frozenset({10, 11, 12})

        # occupancy denominator includes the pocketless conformers:
        # weights are 2,1,1,1 = 5; residue 10 appears in c1 (2) and c2 (1).
        from pwtp import BenchmarkNormaliser, residue_pwtp

        norm = BenchmarkNormaliser().fit(confs)
        occ = residue_pwtp(confs, norm)[10].occupancy
        assert abs(occ - 3 / 5) < 1e-12, occ

        # conflicting weights for one conformer are rejected
        bad = Path(d) / "bad.csv"
        bad.write_text(
            "conformer_id,pocket_id,volume,hydrophobicity,druggability,residues,weight\n"
            "c1,p1,100,0.8,0.7,10,2.0\n"
            "c1,p2,50,0.4,0.3,20,3.0\n"
        )
        try:
            load_csv(bad)
        except ValueError as exc:
            assert "conflicting weights" in str(exc), exc
        else:
            raise AssertionError("expected ValueError")

        # missing column is reported, not silently tolerated
        short = Path(d) / "short.csv"
        short.write_text("conformer_id,volume\nc1,10\n")
        try:
            load_csv(short)
        except ValueError as exc:
            assert "missing columns" in str(exc), exc
        else:
            raise AssertionError("expected ValueError")

        # a bad residue field names its line
        badres = Path(d) / "badres.csv"
        badres.write_text(
            "conformer_id,pocket_id,volume,hydrophobicity,druggability,residues\n"
            "c1,p1,100,0.8,0.7,abc\n"
        )
        try:
            load_csv(badres)
        except ValueError as exc:
            assert ":2" in str(exc), exc
        else:
            raise AssertionError("expected ValueError")

    print("load.py self-test: all checks passed")


if __name__ == "__main__":
    _selftest()
