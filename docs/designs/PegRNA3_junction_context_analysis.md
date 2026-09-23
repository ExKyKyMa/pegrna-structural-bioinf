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

## Spacer-constrained refolding

The free-RNA fold above lets the spacer (1-21) base-pair
intramolecularly with the PBS region. In the assembled prime editor the
spacer is hybridised to the target DNA strand, so this pairing cannot
occur. The analysis was therefore repeated with the spacer forced
unpaired:

```bash
python3 scripts/analyze_junction_neighborhood.py 35 77 121 \
    --label three_linker --unpaired-spacer
python3 scripts/analyze_junction_neighborhood.py 35 97 \
    --label two_linker --unpaired-spacer
```

Removing the artifact costs 21.4 kcal/mol of spurious stability:

| Model | MFE | Ensemble free energy |
|---|---:|---:|
| Free RNA | -54.30 | -57.45 |
| Spacer unpaired | -32.90 | -34.86 |

### Junction P(paired) in both models

| Junction | Position | Free RNA | Spacer unpaired | Dominant partner (constrained) |
|---|---:|---:|---:|---|
| 35/36 | 34 | 0.011 | 0.000 | – |
| 35/36 | 35 | 0.011 | 0.000 | – |
| 35/36 | 36 | 0.000 | 0.000 | – |
| 35/36 | 37 | 0.000 | 0.000 | – |
| 77/78 | 76 | 0.029 | 0.031 | 102 (0.031) |
| 77/78 | 77 | 0.149 | 0.167 | 103 (0.145) |
| 77/78 | 78 | 0.775 | 0.870 | 73 (0.686) |
| 77/78 | 79 | 0.812 | 0.911 | 72 (0.712) |
| 97/98 | 96 | 0.997 | 0.985 | 84 (0.974) |
| 97/98 | 97 | 0.988 | 0.976 | 83 (0.965) |
| 97/98 | 98 | 0.201 | 0.054 | 104 (0.022) |
| 121/122 | 120 | 1.000 | 0.215 | 108 (0.144) |
| 121/122 | 121 | 0.998 | 0.827 | 144 (0.806) |
| 121/122 | 122 | 0.997 | 0.938 | 143 (0.938) |
| 121/122 | 123 | 0.000 | 0.945 | 142 (0.945) |

### Revised conclusions

1. **35/36 is robustly unpaired.** P(paired) drops to 0.000 in the
   constrained model. The conclusion does not depend on the folding
   model, which agrees with the experimentally tolerated sgRNA junction
   34/35.

2. **97/98 is model-independent and genuinely paired.** The 83-86 /
   94-97 stem survives the constraint at P = 0.976-0.985, confirming it
   as a real structural element rather than a folding artifact. It
   corresponds to the stem loop G82-C96 that SpCas9 recognises directly
   in the cryo-EM termination structure. Note that the 3′ side of the
   junction becomes cleaner (98: 0.201 -> 0.054), so the linker strains
   the stem only from its 5′ side.

3. **121/122 was misdiagnosed as a pure artifact.** The spacer pairs
   (6-10) disappear as expected, but an extension-Evopreq1 helix takes
   their place: 121-144 (0.806), 122-143 (0.938), 123-142 (0.945). The
   junction is therefore paired in both models, only with different
   partners. A linker here breaks the extension-to-Evopreq1 stem and
   may compromise the protective function of the 3′ motif.

4. **77/78 becomes slightly worse under the constraint**
   (78: 0.775 -> 0.870; 79: 0.812 -> 0.911), because the spacer no
   longer competes for these partners. The break point itself (76-77)
   stays weakly paired, so the site remains acceptable.

### Caveat

Forcing the spacer unpaired is still only an approximation of the
protein-bound state. The PBS is also hybridised to the nicked target
DNA strand, and the scaffold is clamped by SpCas9. Secondary-structure
prediction should be used as a filter, not as a verdict: sites unpaired
in both models are safe, sites paired in both models carry real risk,
and sites that change between models require explicit 3D evaluation.
