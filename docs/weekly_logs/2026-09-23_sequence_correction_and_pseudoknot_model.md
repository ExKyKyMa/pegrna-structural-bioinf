# 2026-09-23 — Sequence correction, junction context, and the tevopreQ1 pseudoknot model

## 1. Sequence audit and correction

A full audit of the PegRNA3 reference against the canonical
spacer / scaffold / extension / 3'-motif architecture found **two
compensating errors** that cancelled in total length and therefore passed
every previous length check.

| Region | Was | Now | Coordinates |
|---|---|---|---|
| Spacer | ACAGGGUGUGGCAGAAGC (18 nt) | ACAGGGUGUGGCAGAAGCAGC (21 nt) | 1-21 |
| Scaffold fragment 1 | UUUUAGAGCUAGA (13 nt) | GUUUUAGAGCUAGA (14 nt) | 22-35 |
| Scaffold fragment 2 | unchanged (62 nt) | unchanged | 36-97 |
| Extension | GUGCGUCUUUCCUGGUGCUUCUGCCACA (28 nt) | GUCUUUCCUGGUGCUUCUGCCACA (24 nt) | 98-121 |
| 3' motif | unchanged (37 nt) | unchanged | 122-158 |

Cause: `AGCG` was missing after spacer position 18, and a duplicate `GUGC`
had been prepended to the extension. Total length stayed 158 nt in both
cases, so the two errors hid each other.

Consequence: all previously reported junction coordinates were wrong by
4 nt. Junctions had effectively been analysed at 31/32 and 93/94 instead
of the intended 35/36 and 97/98.

Every derived result was regenerated after the fix: MFE fold, secondary
structure drawings, dot plot, base-pair probability table and the
10,000-structure Boltzmann ensemble.

| Quantity | Before | After |
|---|---:|---:|
| MFE | -52.00 | -54.30 |
| Ensemble free energy | -54.91 | -57.45 |
| MFE frequency | — | 0.0060157 |
| Ensemble diversity | — | 22.59 |

Commit: `d079da4`.

## 2. Junction neighbourhood analysis

New script `scripts/analyze_junction_neighborhood.py` reports, for a
4-nt window on each side of every proposed linker junction, the MFE
partner and the ensemble probability that each nucleotide is paired.

A `--unpaired-spacer` flag was added, because in vivo the spacer (1-21)
is base-paired to target DNA and is therefore unavailable for
intramolecular pairing. Ignoring this adds 21.4 kcal/mol of spurious
predicted stability (MFE -54.30 free versus -32.90 constrained).

Commits: `6ad23c0`, `7637cd0`.

## 3. Design finding: the three-linker construct has no anti-over-extension junction

Shuto et al. show that M-MLV RT does not stop at the end of the RTT but
extends into the scaffold as far as U94, three nucleotides upstream of
the RTT 5' end, producing scaffold-derived 1-3 nt incorporations.

Blocking this requires a non-templated chemical junction at the
scaffold-to-extension boundary, i.e. at **97/98**. The two-linker
construct (35/36 and 97/98) has it. The three-linker construct
(35/36, 77/78, 121/122) does **not**: position 97/98 lies inside its
fragment 3 (78-121) as an ordinary phosphodiester bond.

A four-fragment alternative with junctions **35/36 + 77/78 + 97/98**,
leaving the 3' motif attached to the extension, keeps the protective
junction and is recommended for consideration.

## 4. Identification of the 3' motif

The 37-nt 3' motif of PegRNA3 is **tevopreQ1** (trimmed evopreQ1), a
prequeuosine1 class I riboswitch aptamer. The repository currently
labels it "Evopreq1"; renaming is pending and deliberately not done yet.

Published design context that matters for this project: Nelson et al.
place an **8-nt linker** between the PBS and the motif specifically to
prevent the motif from interfering with the rest of the pegRNA, and
wrote **pegLIT** to search for linker sequences that minimise base
pairing with the spacer, PBS, template and scaffold. PegRNA3 currently
has no nucleotide linker there at all, so a chemical junction at 121/122
occupies the slot that the published design reserves for a deliberately
non-interacting element.

Commit: `807d893`.

## 5. Correction of an earlier interpretation

An earlier conclusion stated that a linker at 121/122 would break a
large extension-to-motif stem. This was wrong. Folding the motif in
isolation and in full context gives the **identical** pair set
(122-143, 123-142, 124-141, 125-140, 129-136, 130-135); the only
cross-boundary pair is the single pair **121-144**. The junction removes
one spurious terminal pair and leaves the motif's own fold intact.

## 6. Folding with the pseudoknot enforced

tevopreQ1 folds as an H-type pseudoknot, and ViennaRNA cannot represent
pseudoknots. New script `scripts/fold_with_motif_pseudoknot.py` encodes
the consequence of the knot instead: nucleotides engaged in the
pseudoknot are marked unpaired, and everything else is folded around
them. Three candidate stem-2 registers were enumerated from sequence
complementarity; two of them contain the G15-C32 pair, consistent with
the published observation that the G15C point mutation abolishes
pseudoknot formation.

Results, spacer constrained unpaired:

| Model | MFE | Junction 121/122, P(paired) |
|---|---:|---|
| Motif free to misfold | -32.90 | 0.938 / 0.945 |
| Pseudoknot enforced, register A | -24.00 | 0.100 / 0.000 |
| register B | -24.00 | 0.100 / 0.000 |
| register C | -25.70 | 0.114 / 0.000 |

The predicted 8-bp duplex between scaffold 56-65 and the motif tail
147-156 disappears completely once the knot is enforced, and about
9 kcal/mol of apparent stability disappears with it. Junction 121/122
becomes the lowest-risk junction in the molecule, independently of the
assumed register. Junction 97/98 stays at P = 0.98 in every model.

Full analysis: `docs/designs/PegRNA3_pseudoknot_constrained_analysis.md`.
Commit: `7462f9a`.

## 7. Current risk ranking of candidate junctions

| Junction | Structural risk | Functional value |
|---|---|---|
| 35/36 | low, unpaired in every model, matches the sgRNA 34/35 positive control | segments the scaffold |
| 121/122 | low once the pseudoknot is modelled correctly | isolates the 3' motif, occupies the pegLIT slot |
| 77/78 | moderate, 3' side sits in a P ~ 0.8 helix | segments the scaffold |
| 97/98 | high, real SpCas9-recognised stem loop 83-86 : 94-97 at P = 0.98 | the only junction that can block RT over-extension |

## 8. Tooling notes

- ViennaRNA 2.7.2 installed as a Python module (`pip install ViennaRNA`).
  The command-line programs RNAfold and RNAsubopt are not included in
  the wheel; the Python API is used throughout.
- `fold_compound.hc_add_up_batch` does not exist in these bindings.
  Hard constraints must be added per position with
  `fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)`.
- `fc.pbacktrack(n)` returns nothing unless `md.uniq_ML = 1` is set
  before creating the fold compound.

## 9. Known limitations of the current model

1. No pseudoknot support in the energy model. Handled approximately, not
   solved. The residual predicted contact between the motif loop and
   scaffold 26-35 moves when the assumed register moves, which marks it
   as an artifact.
2. No protein. The real pegRNA is bound to SpCas9 over G21-C96, which
   enforces the scaffold fold rather than allowing it to be predicted.
3. No target DNA except as the `--unpaired-spacer` approximation.
4. No Mg2+ and no tertiary structure.
5. The extension (98-121) is not yet split into RTT and PBS, so no
   edit-specific analysis is possible yet.

## 10. Next step

Chemical definition and validation of the linkers with RDKit, then
comparison of the linker span against the phosphate-to-phosphate
distances measured at each junction in the cryo-EM structures
(PDB 8WUS, 8WUT, 8WUU, 8WUV, 8YGJ). See `docs/ROADMAP.md`.

## Sources

- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. *Nature* **631**, 224-231 (2024).
  https://doi.org/10.1038/s41586-024-07497-8
  PDB accession codes 8WUS, 8WUT, 8WUU, 8WUV, 8YGJ.
  Open-access version: https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
- Nelson, J. W. et al. Engineered pegRNAs improve prime editing
  efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
- Structure and function of preQ1 riboswitches.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4177978/
- pegLIT. https://github.com/sshen8/peglit and https://peglit.liugroup.us/
