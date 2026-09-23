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

## Reassessment of junction 121/122: the 3' motif should be structurally autonomous

The 3' motif of PegRNA3, `CGCGGUUCUAUCUAGUUACGCGUUAAACCAACUAGAA`
(122-158, 37 nt), is exactly **tevopreQ1** (trimmed evopreQ1), the
prequeuosine1 riboswitch aptamer used as the standard epegRNA 3' motif
(Nelson et al. 2022). Its purpose is to fold into an autonomous
structure that blocks 3'->5' exonucleases; any base-pairing with the
rest of the pegRNA is by definition unwanted.

### Published design intent supports insulating the motif

Nelson et al. inserted an 8-nt linker between the PBS and the motif
specifically to prevent interference, and developed **pegLIT** to
search for linker sequences that *minimize* base pairing between the
linker and the spacer, PBS, template and scaffold. In other words, the
position 121/122 is precisely the slot where the published design puts a
deliberately non-interacting spacer element.

A triazole linker is the limiting case of that idea: a non-natural
backbone cannot base-pair at all, so it is non-interfering by
construction rather than by sequence optimisation. PegRNA3 currently has
no nucleotide linker at all between the extension and tevopreQ1, so the
chemical junction fills a gap rather than creating one.

### The tevopreQ1 stem is unaffected by the junction

Folding 122-158 in isolation and in the full molecule (spacer
constrained unpaired) gives the **identical** pair set:

```
122-143, 123-142, 124-141, 125-140, 129-136, 130-135
```

Only a single cross-boundary pair exists: **121-144** (P = 0.806), an
extension nucleotide extending the motif stem by one base pair. The
earlier statement that a linker at 121/122 "breaks the
extension-to-Evopreq1 stem" was therefore wrong: the linker removes one
spurious terminal pair and leaves the motif's own fold intact. This
makes 121/122 a **low-risk, design-consistent** junction.

### The real interference is long-range and a linker cannot block it

In the full molecule the 3' tail of tevopreQ1 is sequestered by the
scaffold in an 8-bp duplex:

| Motif tail | Scaffold partner | Probability |
|---|---|---:|
| 147 A | 65 U | 0.643 |
| 148 A | 64 U | 0.743 |
| 149 C | 63 G | 0.762 |
| 152 A | 60 U | 0.821 |
| 153 C | 59 G | 0.844 |
| 154 U | 58 A | 0.844 |
| 155 A | 57 U | 0.843 |
| 156 G | 56 C | 0.835 |

Insulating the molecule at 121/122 (folding 1-121 and 122-158 as
separate chains) removes this duplex and costs 4.8 kcal/mol of
predicted stability (-32.90 vs -28.10 kcal/mol). However, a triazole
linker keeps the chain covalently continuous, so it **cannot** prevent
long-range pairing between positions 56-65 and 147-156. The insulated
fold is an upper bound on what a chemical junction can achieve, not a
prediction of it.

### Why this predicted interference is probably overestimated

tevopreQ1 is a **pseudoknotted** aptamer, and ViennaRNA cannot predict
pseudoknots. In the real fold the 3' tail participates in the pseudoknot
helix and is unavailable for intermolecular pairing. The scaffold-tail
duplex predicted above is therefore likely an artifact of a
pseudoknot-blind model, in the same way the spacer-PBS duplex was an
artifact of ignoring the target DNA.

### Revised risk ranking

| Junction | Risk | Basis |
|---|---|---|
| 35/36 | low | unpaired in both models; sgRNA positive control |
| 121/122 | low | breaks one spurious pair; matches pegLIT design slot |
| 77/78 | moderate | 5' side unpaired, 3' side in a P ~ 0.7 helix |
| 97/98 | high structurally, high value functionally | real SpCas9-recognised stem loop 83-86/94-97; only junction able to block RT over-extension |

### Next analytical step

Model tevopreQ1 with its pseudoknot enforced (for example as hard
constraints from the preQ1 aptamer secondary structure, or with a
pseudoknot-capable predictor), then re-fold the remainder of PegRNA3
with the motif's paired positions constrained. This removes the second
known artifact and gives the first junction risk estimate that is not
biased by a pseudoknot-blind model.

## Sources

- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. *Nature* **631**, 224-231 (2024).
  https://doi.org/10.1038/s41586-024-07497-8
- Nelson, J. W. et al. Engineered pegRNAs improve prime editing
  efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
- pegLIT source code. https://github.com/sshen8/peglit
