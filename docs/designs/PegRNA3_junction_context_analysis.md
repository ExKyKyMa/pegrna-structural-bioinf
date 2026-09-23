# PegRNA3 junction base-pairing context

## Method

- Input: `data/raw/sequences/PegRNA3_unmodified_reference.fasta` (158 nt)
- Software: ViennaRNA 2.7.2 (Python API)
- Temperature: 37 °C
- MFE: -54.30 kcal/mol; ensemble free energy: -57.45 kcal/mol
- Script: `scripts/analyze_junction_neighborhood.py`
- Window: 4 nt on each side of every junction
- `P(paired)` = sum of all ensemble base-pair probabilities >= 0.01 for
  that nucleotide, i.e. the probability that it is paired at all.

## Two-linker construct (junctions 35/36 and 97/98)

| Pos | Base | Region | MFE partner | P(paired) |
|---:|---|---|---:|---:|
| 34 | G | Scaffold F1 | – | 0.011 |
| 35 | A | Scaffold F1 | – | 0.011 |
| 36 | A | Scaffold F2 | – | 0.000 |
| 37 | A | Scaffold F2 | – | 0.000 |
| 94 | G | Scaffold F2 | 86 | 0.997 |
| 95 | U | Scaffold F2 | 85 | 0.997 |
| 96 | G | Scaffold F2 | 84 | 0.997 |
| 97 | C | Scaffold F2 | 83 | 0.988 |
| 98 | G | Extension | 61 | 0.201 |
| 99 | U | Extension | – | 0.431 |
| 100 | C | Extension | 59 | 0.785 |
| 101 | U | Extension | 59/58 | 0.791 |

## Three-linker construct (junctions 35/36, 77/78 and 121/122)

| Pos | Base | Region | MFE partner | P(paired) |
|---:|---|---|---:|---:|
| 34 | G | Scaffold F1 | – | 0.011 |
| 35 | A | Scaffold F1 | – | 0.011 |
| 36 | A | Scaffold F2 | – | 0.000 |
| 37 | A | Scaffold F2 | – | 0.000 |
| 76 | A | Scaffold F2 | – | 0.029 |
| 77 | A | Scaffold F2 | – | 0.149 |
| 78 | A | Scaffold F2 | 73 | 0.775 |
| 79 | A | Scaffold F2 | 72 | 0.812 |
| 120 | C | Extension | 8 | 1.000 |
| 121 | A | Extension | 7 | 0.998 |
| 122 | C | Evopreq1 | 6 | 0.997 |
| 123 | G | Evopreq1 | – | 0.000 |

## Interpretation

### Junction 35/36 — low risk

Nucleotides 34-37 are essentially unpaired (P <= 0.011). The junction
lies in the apical loop of a scaffold hairpin whose closing pairs
32-39 (0.947) and 33-38 (0.912) are outside the linker site. This is
consistent with the experimentally tolerated sgRNA junction 34/35.

### Junction 77/78 — moderate risk

Positions 76-77 are weakly paired (P <= 0.15), but 78-79 participate in
a moderately probable local helix with 72-73 (P ~ 0.61-0.63). The
linker is placed on the 5′ shoulder of that helix, so the covalent
break itself falls in a low-pairing window while the adjacent helix is
only moderately populated in the ensemble. Acceptable, but the
structural cost is higher than at 35/36.

### Junction 97/98 — structurally the most demanding, functionally the most valuable

Positions 94-97 form the 3′ scaffold stem loop with 83-86
(P = 0.988-0.997). In the cryo-EM termination structure of Shuto et al.
(2024) this same stem loop (G82-C96 in their numbering) is directly
recognized by SpCas9: Y64/L99 contact it and R116 forms a hydrogen bond
with C96. Cutting the backbone at 97/98 therefore sits at the edge of a
protein-recognized, highly paired element.

The functional benefit is equally specific: 97/98 is the scaffold /
3′-extension boundary, i.e. exactly where M-MLV RT over-extends reverse
transcription by up to three nucleotides into the scaffold. A
non-natural triazole backbone at this position is the only junction in
either design that can physically terminate that read-through.

### Junction 121/122 — high apparent pairing, but largely a free-RNA artifact

Positions 118-122 pair with 6-10 at P ~ 1.0. This is the spacer pairing
intramolecularly with the PBS region, which cannot occur in the
assembled complex where the spacer is hybridised to the target DNA
strand. The pairing predicted here must therefore be re-evaluated with
the spacer constrained as unpaired before this junction is called risky.

## Consequence for construct design

The three-linker construct does **not** contain a junction at 97/98.
Its junctions are 35/36, 77/78 and 121/122, so position 97/98 lies
inside fragment 3 (78-121) as an ordinary phosphodiester bond. The
anti-over-extension effect is present only in the two-linker construct.

Recommended follow-up designs:

1. Re-fold with the spacer (1-21) constrained as unpaired, to remove the
   spacer-PBS artifact, before comparing junction risks.
2. Consider a four-fragment variant with junctions 35/36, 77/78 and
   97/98, which keeps the low-risk scaffold junction, the moderate
   scaffold junction and the functionally critical RT-termination
   junction.

## Source

Shuto, Y. et al. Structural basis for pegRNA-guided reverse
transcription by a prime editor. *Nature* **631**, 224-231 (2024).
https://doi.org/10.1038/s41586-024-07497-8
