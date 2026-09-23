# PegRNA3 folding with the tevopreQ1 pseudoknot enforced

## Why this analysis exists

ViennaRNA's energy model only admits nested (non-crossing) base pairs, so it
cannot represent a pseudoknot. The 3' motif of PegRNA3 (positions 122-158) is
**tevopreQ1**, a prequeuosine1 class I riboswitch aptamer that folds as an
H-type pseudoknot. Every unconstrained fold of PegRNA3 therefore mis-folds the
motif and then lets its "leftover" nucleotides pair with the scaffold, which
inflates the apparent structural risk of the 121/122 junction.

This analysis removes that bias. Script: `scripts/fold_with_motif_pseudoknot.py`.

## Method

A pseudoknot cannot be written as a hard constraint, but its *consequence* can:
nucleotides already engaged in the pseudoknot are unavailable for pairing with
the rest of the molecule. Those positions are marked unpaired
(`fc.hc_add_up(...)`) and ViennaRNA folds everything else around them. The
spacer (1-21) is optionally constrained the same way, because in vivo it is
base-paired to target DNA.

### Pseudoknot model

Stem 1 is unambiguous, taken from the isolated fold of the motif:

| Motif-local | Global | Pairs |
|---|---|---|
| 1-4 : 19-22 | 122-125 : 140-143 | CGCG / CGCG |

Stem 2 (the pseudoknot helix, formed between the stem-1 loop and the 3' tail)
is not published as an explicit diagram, so three candidate registers were
enumerated from sequence complementarity and all three were tested:

| Register | Motif-local | Global | Loop 1 | Contains G15 |
|---|---|---|---|---|
| A (default) | 11-16 : 31-36 | 132-137 : 152-157 | 6 nt | yes |
| B | 12-17 : 30-35 | 133-138 : 151-156 | 7 nt | yes |
| C | 6-10 : 33-37 | 127-131 : 154-158 | 1 nt | no |

Registers A and B contain the **G15-C32** pair. Nelson et al. report that the
single point mutation **G15C abolishes pseudoknot formation** and the entire
epegRNA benefit, which is direct evidence that G15 is base-paired in the knot;
registers A and B are therefore preferred and A is the script default.
Register C matches the preQ1-I consensus loop-1 length (2-3 nt) better but
leaves G15 unpaired, so it is retained only as a sensitivity test.

## Results

Predicted stability, spacer constrained unpaired in all cases:

| Model | MFE | Ensemble FE |
|---|---:|---:|
| Motif left free to misfold | -32.90 | -34.86 |
| Pseudoknot enforced, register A | **-24.00** | -26.03 |
| Pseudoknot enforced, register B | -24.00 | -26.02 |
| Pseudoknot enforced, register C | -25.70 | -27.92 |

About 9 kcal/mol of the earlier value came from the motif pairing with the
scaffold instead of with itself. That is not real stability.

### The scaffold-to-3'-tail duplex disappears

In the unconstrained fold the motif's 3' tail was captured by the scaffold in
an 8-bp duplex (56-65 with 147-156, P = 0.64-0.84). With the pseudoknot
enforced that duplex is gone entirely:

| Motif region (register A) | Cross-boundary probability mass |
|---|---:|
| loop 2, 138-139 | 0.001 |
| loop 3, 144-151 | 0.368 |
| 3' nucleotide 158 | 0.000 |
| stem nucleotides 152-157 | 0 by construction |

Scaffold positions 56-65 now pair locally within the scaffold
(64-53 at P = 0.65, 65-52 at P = 0.70) instead of reaching into the motif.
The predicted interference reported in the previous document was an artifact of
a pseudoknot-blind model.

### Junction risk under the corrected model

P(paired), spacer constrained unpaired, register A:

| Junction | 5' side | 3' side | Change vs motif-free model |
|---|---|---|---|
| 35/36 | 0.051 | 0.023 | unchanged, still safe |
| 77/78 | 0.511 | 0.764 | slightly worse |
| 97/98 | 0.980 | 0.133 | unchanged, stem loop is real |
| **121/122** | **0.100** | **0.000** | **large improvement** |

Junction 121/122 becomes the cleanest cut site in the molecule. This holds for
all three registers (P(121) = 0.100, 0.100, 0.114), so the conclusion does not
depend on the register choice. The user hypothesis that the motif should be
structurally autonomous, and that a non-pairing chemical linker at 121/122 is
consistent with that role rather than harmful, is supported.

### One residual interaction, and why it is probably also an artifact

Whichever register is assumed, roughly 6 nt of pseudoknot *loop* remain formally
free, and ViennaRNA pairs them with the scaffold repeat region:

| Register | Free loop nt | Predicted partner | P |
|---|---|---|---|
| A | 126-131 | scaffold 26-31 | 0.72-0.83 |
| C | 132-139 | scaffold 28-35 | 0.83-0.999 |

The predicted partner moves when the assumed register moves, which is the
signature of a model artifact rather than a physical interaction. In a real
H-type pseudoknot the loops are short, cross the helical grooves and are
sterically locked by stacking and A-minor contacts; they cannot open up into a
6-bp duplex with a distant sequence without dismantling the knot. Secondary
structure prediction has no way to know this.

Note also that under register C this residual pairing lands on positions 32-35
and would corrupt junction 35/36, whereas under register A junction 35/36 stays
clean. Junction 35/36 should still be considered safe, since it is unpaired in
the motif-free model, in register A, and in the sgRNA positive control.

## Conclusions

1. The motif is tevopreQ1 and it folds as a pseudoknot. Any PegRNA3 fold that
   does not account for this overestimates motif-to-pegRNA interference.
2. With the pseudoknot enforced, junction **121/122 is the lowest-risk junction
   in the whole molecule** (P = 0.10 / 0.00), not a high-risk one.
3. Junction **97/98 remains high-risk and structurally real** in every model:
   the SpCas9-recognised stem loop 83-86 : 94-97 stays at P = 0.98.
4. Junction **35/36 remains safe**; junction **77/78 is moderate** and gets
   slightly worse once the motif is folded correctly.
5. The remaining predicted motif-to-scaffold contact is register-dependent and
   therefore not trustworthy. Resolving it requires a pseudoknot-capable
   predictor or an experimental probe (SHAPE, DMS-MaPseq), not more RNAfold.

## Reproduce

```bash
python3 scripts/fold_with_motif_pseudoknot.py --register A --unpaired-spacer
python3 scripts/fold_with_motif_pseudoknot.py --register C --unpaired-spacer
```

Outputs are written to `results/secondary_structure/pseudoknot/`.

## Sources

- Nelson, J. W. et al. Engineered pegRNAs improve prime editing efficiency.
  *Nat. Biotechnol.* **40**, 402-410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
  (evopreQ1 is a modified preQ1 riboswitch aptamer forming a pseudoknot;
  G15C disrupts pseudoknot formation; tevopreQ1 is evopreQ1 trimmed by 5 nt;
  an 8-nt pegLIT linker is placed between the PBS and the motif.)
- Structure and function of preQ1 riboswitches.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4177978/
  (preQ1-I aptamers form H-type pseudoknots: a P1 hairpin followed by an
  A-rich tail whose 3' end pairs with the centre of the P1 loop; stem 1 is
  5 bp, stem 2 is 3-4 bp, loop 1 is 2-3 nt, loop 3 is 8-9 nt.)
- Shuto, Y. et al. Structural basis for pegRNA-guided reverse transcription by
  a prime editor. *Nature* **631**, 224-231 (2024).
  https://doi.org/10.1038/s41586-024-07497-8
- pegLIT source code. https://github.com/sshen8/peglit
