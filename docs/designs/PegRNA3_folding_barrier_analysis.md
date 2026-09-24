# PegRNA3: secondary structure, Folding Barrier, and spacer accessibility

Status: computational analysis, 2026-09-24. Every number in this document
is reproduced by `scripts/folding_barrier_analysis.py`; the tables in
`results/folding_barrier/` are the machine-readable source. Nothing here
has been tested experimentally.

This analysis was run out of interest, alongside the planned cell-based
comparison of the linker constructs. It turned up one result that is not
about the linkers at all and that is worth acting on.

---

## Summary

1. The unmodified PegRNA3 has a **Folding Barrier of 32.8 kcal/mol**
   (relaxed definition) or **39.6 kcal/mol** (strict definition).
   In the reference study every effective guide sat below 10 kcal/mol
   and every defective guide above it.
2. The plain sgRNA control — the same spacer and scaffold with the 3'
   extension removed — has a barrier of **12.0 kcal/mol**. The 3'
   extension therefore adds **+20.8 kcal/mol**. The high number is a
   property of this specific extension, not an artefact of the pegRNA
   being long or of the method.
3. The cause is a single, strong, intramolecular duplex: the **spacer
   pairs with the PBS**, ΔG = **−23.4 kcal/mol**, with pairing
   probability ≈ 1.00 across 12 consecutive positions. Mean spacer
   accessibility collapses from **0.355** in the sgRNA to **0.197** in
   PegRNA3; seed accessibility at positions 19–21 collapses from
   **0.917** to **0.339**.
4. This is a known and published failure mode of prime editing, called
   auto-inhibition. It predicts low, not zero, editing efficiency.
5. The linkers **lower** the barrier rather than raising it. The
   three-linker construct drops it by 6.6 kcal/mol, the four-linker
   construct by 10.5 kcal/mol. This is a small predicted benefit on top
   of a large underlying problem.

---

## 1. Why a barrier rather than a stability

The obvious way to ask whether a guide RNA is well folded is to compute
its free energy: the more negative, the more stable. That turns out to be
the wrong question, and the reason matters for how these numbers are read.

A guide RNA has two relevant structures. The **ground state** is what the
free RNA folds into on its own, its minimum free energy (MFE) structure.
The **active state** is the conformation Cas9 requires, which is known
experimentally from the cryo-EM structure of the prime editor. These are
not the same structure. A very stable ground state is bad, not good,
because the molecule has to leave it.

Zalatan and co-workers formalised this as the **Folding Barrier**: the
energy of the highest point on the path from the ground state to the
active state. Working on CRISPRa guide RNAs, they found it correlates with
activity at Spearman r_s = 0.8, explaining about 80% of the variation, and
that the separation was clean — every effective guide had a barrier below
roughly 10 kcal/mol, every defective guide above it. Folding Energy, the
plain energy gap, gave r_s = 0.7. Binding Energy alone was insufficient:
33% of the guides with the most favourable Binding Energy were still
defective. General-purpose guide design scores did far worse on the same
dataset: Azimuth r_s = 0.22, Doench 2016 r_s = 0.02, Moreno-Mateos
r_s = 0.09 ([Nature Communications 15, 6341,
2024](https://www.nature.com/articles/s41467-024-50528-1)).

The barrier is computed with the `findpath` heuristic in ViennaRNA: it
searches direct refolding paths between two structures and reports the
saddle point ([ViennaRNA Package 2.0, Algorithms for Molecular Biology 6,
26, 2011](https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26)).
Because it is a heuristic over a limited number of paths, **every barrier
reported here is an upper bound** — the true barrier can only be lower.

### The threshold does not transfer directly

The 10 kcal/mol threshold was calibrated on ~100-nt bacterial CRISPRa
scRNAs. PegRNA3 is 158 nt. Longer RNA has more structure and a rougher
landscape, so a naive comparison would be unfair. That is why this
analysis always computes a control.

---

## 2. Defining the active state

The active state is built from the base pairs actually observed in the
cryo-EM structure of the prime editor ([Shuto et al., Nature 631,
224–231, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/); PDB
8WUS, 8WUT, 8WUU, 8WUV, 8YGJ), stored in
`data/raw/designs/article_pegrna_bound_rna_pairs.tsv`.

The published pegRNA has a 20-nt spacer and PegRNA3 has 21, so every
position shifts by +1. The script verifies that offset at run time by
locating the scaffold start in both sequences rather than hard-coding it.

Of the 24 observed pairs, 22 are usable. Two are dropped, and the reasons
are recorded rather than hidden:

| Pair | Position | Bases | Why it is dropped |
|---|---|---|---|
| P13 | 56–59 | C–G | encloses a 2-nt hairpin loop; the nearest-neighbour model forbids loops under 3 nt, and the conformation exists only because the protein strains the backbone |
| P14 | 69–82 | A–G | non-canonical pair, the energy model has no parameters for it |

Neither is a failure of the comparison. Both are cases where a protein
enforces geometry that free RNA cannot adopt.

Two definitions of the active state are reported throughout, because the
choice changes the number and it would be dishonest to pick one silently:

- **Strict** — only the cryo-EM pairs; every other nucleotide unpaired.
  This overstates the barrier, because the real bound molecule is not
  bare outside the resolved core.
- **Relaxed** — the cryo-EM pairs enforced, the spacer held unpaired
  because Cas9 holds it extended for DNA interrogation, and everything
  else allowed to fold. More physical, and the one used for the headline
  numbers.

---

## 3. Results

### The control comes first

| | plain sgRNA (1–97) | PegRNA3 (1–158) |
|---|---|---|
| length | 97 nt | 158 nt |
| ground state energy | −24.3 | −54.3 |
| active state energy, relaxed | −14.9 | −27.0 |
| Folding Energy, relaxed | 9.4 | 27.3 |
| **Folding Barrier, relaxed** | **12.0** | **32.8** |
| Folding Barrier, strict | 12.0 | 39.6 |
| mean spacer accessibility | 0.355 | **0.197** |
| seed accessibility, 19–21 | 0.917 | **0.339** |

All energies in kcal/mol.

The sgRNA control sits at 12.0 kcal/mol, close to the 10 kcal/mol
threshold, which is what a working guide should look like and which
suggests the method is behaving sensibly on a 97-nt RNA. Adding the 3'
extension raises the barrier by **+20.8 kcal/mol** and drops seed
accessibility from 0.917 to 0.339.

That is a large, specific effect, and it is not explained by length alone:
truncating PegRNA3 at position 121, removing only the tevopreQ1 motif and
keeping the extension, still gives a barrier of 27.4 kcal/mol. The
extension, not the motif, carries most of the cost.

### The linker variants

| Variant | Junctions | Ground | Barrier, strict | Barrier, relaxed | Δ relaxed | Spacer acc. |
|---|---|---|---|---|---|---|
| unmodified | — | −54.3 | 39.6 | 32.8 | — | 0.197 |
| two_linker | 35/36, 97/98 | −50.2 | 37.0 | 29.4 | −3.4 | 0.197 |
| three_linker | 35/36, 77/78, 121/122 | −47.5 | 34.3 | **26.2** | **−6.6** | 0.313 |
| four_linker | 35/36, 77/78, 97/98, 121/122 | −43.6 | 30.4 | **22.3** | **−10.5** | 0.313 |
| 35/36 only | 35/36 | −54.3 | 39.6 | 32.8 | 0.0 | 0.197 |
| 77/78 only | 77/78 | −53.9 | 40.7 | 33.9 | +1.1 | 0.197 |
| 97/98 only | 97/98 | −50.2 | 37.0 | 29.4 | −3.4 | 0.197 |
| 121/122 only | 121/122 | −47.9 | 33.2 | 26.6 | −6.2 | 0.313 |

Reading these:

- **35/36 is thermodynamically silent.** Zero change under both
  definitions. Consistent with the wet-lab result that this junction does
  not reduce CRISPR/Cas9 efficiency, and consistent with the low pairing
  probability previously computed for these positions.
- **77/78 is the only junction that raises the barrier**, by
  +1.1 kcal/mol. Small, but it is the one junction retained for
  convenience rather than validated, and it is the one where the
  independent pairing analysis gave the highest occupancy
  (P = 0.511 / 0.764). Two independent calculations flag the same
  junction. 57/58 remains the registered fallback.
- **121/122 gives the largest single-junction improvement**, −6.2, by
  decoupling the tevopreQ1 motif from the rest of the molecule.
- The improvements are not additive, because they act on overlapping
  structure.

The mechanism of the improvement is straightforward: forcing a nucleotide
unpaired destabilises the ground state more than it destabilises the
active state, so the molecule has less to unfold. This is a prediction
about **assembly kinetics**, not about catalysis, and it is unverified.

### Caveat on how the linker is modelled

A triazole linker is represented here as a hard constraint that the two
flanking nucleotides cannot pair. That captures the loss of a
phosphodiester bond and the inability of the triazole to stack into a
helix. It does not capture the linker's own geometry, its flexibility, or
any steric effect on the protein. The direction of the effect is probably
right; the magnitudes should not be quoted as precise.

---

## 4. What is actually wrong: the spacer is sequestered by the PBS

Decomposing where the spacer's pairing probability goes:

| Partner region | Total pairing probability |
|---|---|
| **PBS (110–121)** | **11.995** |
| motif (122–158) | 2.883 |
| RTT (98–109) | 1.956 |
| scaffold (22–97) | 0.027 |
| spacer (self) | 0.000 |

Out of 21 spacer nucleotides, essentially 12 full units of pairing
probability go to the PBS. Position by position, spacer 7–18 pair with
PBS 121–110 at probability 0.998 to 1.000, and spacer 19 and 21 pair into
the RTT at 0.98 and 0.96. Almost the entire spacer is locked up.

The duplex energies:

- spacer : PBS — **ΔG = −23.4 kcal/mol**
- spacer : RTT + PBS — **ΔG = −25.2 kcal/mol**

### This is expected by construction, and it is documented

The logic is unavoidable given how a pegRNA is designed. The PBS is
complementary to the nicked target DNA strand. That strand is
complementary to the protospacer. The protospacer is the spacer sequence.
So the spacer and the PBS are complementary to each other by construction
— and an RNA:RNA duplex is more stable than the RNA:DNA duplex the PBS is
supposed to form with the DNA. The intramolecular reaction also wins on
effective concentration.

Ponnienselvan and co-workers described exactly this and named it
auto-inhibition. Their findings ([Nucleic Acids Research 51, 6966–6980,
2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10359601/)):

- pegRNAs with a standard 13- or 14-nt PBS gave **negligible or no
  detectable in vitro Cas9 cleavage** of the target.
- Adding an oligonucleotide complementary to the PBS–RTT region
  **restored** cleavage, proving the inhibition was the intramolecular
  duplex and not something else.
- Shortening the PBS restored activity. For synthetic, chemically
  end-protected pegRNAs and for epegRNAs carrying a 3' pseudoknot, the
  optimal PBS was as short as **7 nt**, against 13–14 nt for
  plasmid-expressed pegRNAs.
- Complementarity can extend into the first three RTT nucleotides when
  those match the target — which is what positions 19 and 21 pairing into
  the RTT show here.
- Designing the PBS to a predicted PBS:target melting temperature near
  37 °C gave 49.3% editing with PE2 and 73% with PE3 at one site.

The relevance is direct. PegRNA3 has a **12-nt PBS**, is **chemically
modified**, and carries a **3' pseudoknot (tevopreQ1)** — all three
features that in that study shifted the optimum toward a shorter PBS.

Two further papers point the same way. Yao and colleagues found that a
**more** stable 3' extension gave **worse** editing: ΔG going from −9.6 to
−12.5 kcal/mol and stems from 3 to 6 dropped efficiency, while adding CGC
to reduce stem count restored editing from 4.68% to 7.95% ([FEBS Letters,
2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12421704/)). Wong and
colleagues, on ordinary sgRNAs, found that the single most discriminating
structural feature between functional and non-functional guides was the
accessibility of spacer positions 18–20, the 3' seed — non-functional
guides sequestered them in an extended stem-loop, with mean self-folding
ΔG of −3.1 versus −1.9 kcal/mol, P = 6.7 × 10⁻¹¹ ([Genome Biology 16, 218,
2015](https://link.springer.com/article/10.1186/s13059-015-0784-0)).
PegRNA3's seed accessibility of 0.339, against 0.917 for the same spacer
without the extension, falls on the wrong side of that distinction.

There is also a kinetic argument that this may be less severe in reality
than in the calculation. Single-molecule FRET shows Cas9 actively shifts
sgRNA into a compact folded state, increasing the folded population by
113% and eliminating a slow-folding fraction that otherwise takes 15.8 s
to resolve ([Journal of Physical Chemistry B, 2022](https://pubs.acs.org/doi/10.1021/acs.jpcb.2c05428)).
Protein binding can pull an RNA out of a trap that a free-RNA calculation
says it should be stuck in. That is a reason for uncertainty, not a reason
to dismiss the result — the free state still determines how often the
complex assembles at all.

### What shortening the PBS would predict

A design calculation only. The present construct is already synthesised.

| PBS length | spacer:PBS ΔG | mean spacer accessibility |
|---|---|---|
| 12 nt (current) | −23.4 | 0.197 |
| 11 nt | −21.6 | 0.252 |
| 10 nt | −19.1 | 0.267 |
| 9 nt | −17.3 | 0.297 |
| 8 nt | −14.0 | 0.234 |
| 7 nt | −10.8 | 0.297 |
| 6 nt | −8.3 | 0.358 |

Accessibility improves but does not recover to the sgRNA level of 0.355
until the PBS is down to 6 nt, and it does not improve monotonically —
shortening the PBS changes which alternative structures compete, so the
trend is bumpy. The published optimum of 7 nt for chemically modified,
pseudoknot-carrying pegRNAs sits in a reasonable part of this range, but
a shorter PBS also primes reverse transcription less efficiently, so
there is a real trade-off that thermodynamics alone does not resolve.

---

## 5. What this does and does not say

**Supported by the calculation:**

- PegRNA3's free state is dominated by a spacer:PBS duplex that leaves
  the spacer largely inaccessible.
- The 3' extension, not length in general and not the tevopreQ1 motif,
  is responsible.
- This matches a documented, experimentally demonstrated inhibition
  mechanism in prime editing.
- The linkers reduce the barrier rather than increasing it; 77/78 is the
  single exception and it is the unvalidated junction.

**Not supported, and should not be claimed:**

- That PegRNA3 will not work. Auto-inhibition predicts reduced
  efficiency, not zero. The magnitude is not predictable from these
  numbers.
- That the 10 kcal/mol threshold applies to a 158-nt pegRNA. It was
  calibrated on a different RNA class in a different organism.
- That the linkers will improve editing in cells. The calculation
  describes free-RNA folding thermodynamics, not catalysis, not RT
  processivity, not nuclear delivery.
- Any precise magnitude for the linker effects. The triazole model is
  crude.

**Useful consequence for the planned experiment.** If the cell-based
comparison shows low absolute efficiency for all constructs including the
unmodified control, that is consistent with auto-inhibition and is a
property of the pegRNA design, not of the linker chemistry. The unmodified
control is what separates those two explanations, which is a good reason
to keep it in the experiment.

A falsifiable prediction worth recording before the experiment is run:
ranked by predicted assembly, four_linker > three_linker ≈ 121/122 >
two_linker ≈ 97/98 > unmodified ≈ 35/36 > 77/78. If the measured ranking
is unrelated to this, the folding-barrier model does not describe this
system, which is itself a result.

---

## Reproducing

```bash
python3 scripts/folding_barrier_analysis.py
```

Outputs:

- `results/folding_barrier/PegRNA3_folding_barrier.tsv`
- `results/folding_barrier/PegRNA3_spacer_pairing_partners.tsv`
- `results/folding_barrier/PegRNA3_pbs_length_scan.tsv`

Inputs, which are the source of truth for every sequence:

- `data/raw/sequences/PegRNA3_unmodified_reference.fasta`
- `data/raw/sequences/article_pegrna_reference.fasta`
- `data/raw/designs/article_pegrna_bound_rna_pairs.tsv`

Computed with ViennaRNA 2.7.2, Turner 2004 parameters, 37 °C, findpath
search width 400.

---

## References

- Zalatan, J. G. et al. Guide RNA structure design enables combinatorial
  CRISPRa programs for biosynthetic profiling. *Nature Communications*
  **15**, 6341 (2024).
  https://www.nature.com/articles/s41467-024-50528-1
- Ponnienselvan, K., Liu, P., Nyalile, T., Oikemus, S., Maitland, S. A.,
  Lawson, N. D., Luban, J. & Wolfe, S. A. Reducing the inherent
  auto-inhibitory interaction within the pegRNA enhances prime editing
  efficiency. *Nucleic Acids Research* **51**, 6966–6980 (2023).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10359601/
- Wong, N., Liu, W. & Wang, X. WU-CRISPR: characteristics of functional
  guide RNAs for the CRISPR/Cas9 system. *Genome Biology* **16**, 218
  (2015).
  https://link.springer.com/article/10.1186/s13059-015-0784-0
- Yao, X. et al. Secondary structure of the pegRNA 3' extension affects
  prime editing efficiency. *FEBS Letters* (2025).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12421704/
- Okafor, I. C. & Ha, T. Single molecule FRET analysis of CRISPR Cas9
  single guide RNA folding dynamics. *Journal of Physical Chemistry B*
  (2022).
  https://pubs.acs.org/doi/10.1021/acs.jpcb.2c05428
- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. *Nature* **631**, 224–231 (2024).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
- Nelson, J. W. et al. Engineered pegRNAs improve prime editing
  efficiency. *Nature Biotechnology* **40**, 402–410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
- Lorenz, R. et al. ViennaRNA Package 2.0. *Algorithms for Molecular
  Biology* **6**, 26 (2011).
  https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26
