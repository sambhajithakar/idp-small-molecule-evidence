# Evidence-graded assessment of small-molecule engagement of disordered proteins — data and code

Reference implementation, curated data and analysis scripts for:

> S. B. Thakar and P. N. Ghorpade, *AI-Driven Small-Molecule Discovery against
> Disordered Cancer Proteins: An Evidence-Graded Assessment* (submitted).

Everything here is **freely available for non-commercial use**. No new experimental
data were generated in this work; every measurement in the curated table was published
by another group and each row records the PMID or DOI of its source.

## What is in this repository

```
pwtp/                     the metric, its loader and the synthetic demonstration
  pwtp.py                 PWTP score: occupancy x conditional quality, per residue
  load.py                 CSV loader for pocket-detector output
  demo_synthetic.py       the synthetic validation reported in the paper
data/                     the curated evidence, as released
  01-evidence-tier-definitions.csv
  02-table1-graded-evidence-entries.csv      Table 1: every graded target-compound pair
  03-table2-tier-counts.csv
  04-mechanism-class-details.csv             the four mechanism classes
  05-table3-affinities.csv                   Table 3: every affinity located
  06-synthetic-validation-summary.csv
  07-indistinguishable-pwtp-pairs.csv        the pair that share a PWTP and mean different things
examples/
  sample-input.csv        a minimal pocket-detector export in the expected schema
  sample-output.txt       the exact output of the command below, for comparison
```

## Running it

Python 3.9 or later. No dependencies beyond the standard library — deliberately, so
that a metric claimed as a contribution can be audited by reading it.

```bash
# reproduce examples/sample-output.txt
cd examples && python3 ../pwtp/run_example.py

# the synthetic validation reported in the paper
cd pwtp && python3 demo_synthetic.py
```

`sample-input.csv` carries one row per (conformer, pocket), with the columns
`conformer_id, pocket_id, volume, hydrophobicity, druggability, residues, weight`.
Conformers in which the detector found no pocket contribute no rows; their count is
passed separately as `pocketless`, and getting that wrong inflates occupancy, which is
the easiest way to overstate a target's druggability.

## How to read the output

PWTP is reported as a **pair**, not as one number:

    PWTP(residue) = occupancy x conditional quality

`sample-output.txt` shows why. Residues 41–44 line a high-quality pocket in a minority
of conformers; residues 80–82 line a poor pocket more often. Two residues can share a
PWTP value and mean entirely different things for a discovery decision — one target
whose campaign must capture a rare state, one whose pockets are not worth pursuing.
`data/07-indistinguishable-pwtp-pairs.csv` is the measured instance of that.

**Scope.** PWTP measures apo-ensemble pocket persistence for the pocket-like mechanism
class only. Applied to an oligomer-stabilising or ensemble-remodelling target it returns
a low score *for the correct reason*, and reading that score as "not druggable" would be
wrong. It is not a predictor of binding, affinity or druggability, and the paper reports
a measured bound on what any score of this kind can resolve. Residue-level pocket
persistence over ensembles has been scored before; §2.4 of the paper states the prior
art rather than claiming the idea.

## Why this is code and data rather than a web server

The Journal of Molecular Biology asks authors who deposit code to say why the resource
is not offered as a web server. Three reasons, all of them properties of this particular
resource:

1. **A server would have to run the pocket detection, and that is the part we do not
   supply.** PWTP consumes a pocket-detector export, it does not produce one. The
   upstream step is a user's own choice of detector, version and parameters over their
   own ensemble — fpocket, MDpocket or another — and a server that fixed that choice
   would silently fix a methodological decision the paper argues should stay explicit.
2. **The input is a user's own ensemble, which is large and often unpublished.** The
   realistic unit of work is thousands of conformers of a target a group has not yet
   published. Uploading that to a third-party server is exactly what such groups avoid,
   and the computation is a few seconds of pure Python on their own machine.
3. **The resource is deliberately auditable rather than operable.** The metric is a
   claimed contribution, so the implementation is dependency-free and short enough to
   read end to end. A server would hide the one thing a reader of the paper should be
   able to check.

Who will download and apply it: groups running conformational ensembles of disordered
or partly disordered targets who want to separate ensemble occupancy from conditional
pocket quality before committing to a campaign; and readers checking the paper's
synthetic validation or re-using the curated evidence table, which is the larger part
of this deposit and is useful independently of the metric.

## Curated data

`data/02-table1-graded-evidence-entries.csv` and `data/05-table3-affinities.csv` are the
paper's Tables 1 and 3 as released. Each row is one **target–compound pair**, not one
target: evidence grade and mechanism can differ between two compounds of a series on one
protein, and the paper gives a case where they do. Every row carries its source
identifier, so each entry can be traced to and re-verified against the primary paper.
Bioactivity records examined in the database sweep are available from ChEMBL and
BindingDB under their own terms.

## Citing this

Cite the paper. If you use the code or the table directly, also cite the archived
release (Zenodo DOI, added on deposit).

## Licence

The code is released under the **MIT licence** (`LICENSE`). The curated data files in
`data/` are released under **CC BY 4.0** (`data/LICENSE-DATA`). Both are free for
non-commercial users, as the Journal of Molecular Biology requires of a deposited
resource. Every row of the curated tables records a measurement published by another
group; attribution under CC BY 4.0 does not replace citing those primary sources.
