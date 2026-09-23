# PegRNA3: target site, intended edit, and the RTT / PBS boundary

Status: verified computationally 2026-09-23. The construct itself is a
synthesis order from a collaborating biology group; the target and edit
below were reconstructed from the sequence and then checked against
public reference sequences. No data file has been changed yet.

## Reference sequences used

- **NM_000345.4** — *Homo sapiens* SNCA transcript variant 1, mRNA.
  https://www.ncbi.nlm.nih.gov/nuccore/NM_000345.4
- **NG_011851.1** — SNCA RefSeqGene (genomic).
  https://www.ncbi.nlm.nih.gov/nuccore/NG_011851.1

CDS numbering below is `c.` numbering of NM_000345.4, with c.1 = A of the
initiator ATG.

## Target site

| Item | Value |
|---|---|
| Gene | SNCA (alpha-synuclein) |
| Protospacer | c.69-c.89, `ACAGGGTGTGGCAGAAGCAGC` (21 nt) |
| Strand | coding / sense strand |
| PAM | c.90-c.92, `AGG` |
| Contiguous in genomic DNA | yes, protospacer+PAM found in NG_011851.1 |
| SpCas9 nick (PAM strand) | between **c.86 and c.87** |

The protospacer matches PegRNA3 positions 1-21 exactly, so the spacer is
correct. The PAM is a canonical `NGG` and the protospacer plus PAM is
contiguous in genomic DNA, so the site does not straddle an intron.

## Intended edit

Reverse transcription was reconstructed nucleotide by nucleotide from the
3' extension. The primer 3'-OH is c.86; RT then reads the pegRNA
3' to 5' and writes the complement.

| | Sequence (c.87-c.98) |
|---|---|
| genomic template | `AGCAGGAAAGAC` |
| synthesised 3' flap | `ACCAGGAAAGAC` |

Exactly one position differs:

| | Value |
|---|---|
| Edit | **c.88G>C** |
| Codon 30 | `GCA` (Ala) to `CCA` (Pro) |
| Protein | **p.Ala30Pro (A30P)** |
| Direction | the construct **installs** A30P, it does not correct it |
| Edit distance from the nick | +2 |
| 3' homology after the edit | 10 nt (c.89-c.98) |

This is therefore a disease-model construct: A30P is a pathogenic
missense variant associated with familial Parkinson's disease. If the
biological intent was to *correct* A30P back to wild type, the construct
as ordered does the opposite and the RTT must be changed. This should be
confirmed with the requesting group before any further modelling.

## The RTT / PBS boundary is 12 + 12, not 11 + 13

The extension occupies positions 98-121, `GUCUUUCCUGGUGCUUCUGCCACA`
(24 nt). An earlier working note recorded the split as RTT = 11 nt
(98-108) and PBS = 13 nt (109-121). That is off by one nucleotide.

A 13-nt PBS would have to anneal to c.75-c.87, which requires a nick
between c.87 and c.88, only 2 nt from the PAM. SpCas9 nicks 3 nt from
the PAM. Enumerating every possible PBS length against the canonical
nick position gives exactly one consistent answer:

| PBS length | Required PBS for a nick at c.86/c.87 | Matches construct |
|---:|---|---|
| 11 | `GCUUCUGCCAC` | no |
| **12** | **`GCUUCUGCCACA`** | **yes** |
| 13 | `GCUUCUGCCACAC` | no |

Correct region map:

| Region | Positions | Length | Sequence | Pairs with / templates |
|---|---|---:|---|---|
| Spacer | 1-21 | 21 | `ACAGGGUGUGGCAGAAGCAGC` | protospacer c.69-c.89 |
| Scaffold | 22-97 | 76 | canonical SpCas9 scaffold + `GGUGC` | bound by SpCas9 over G21-C96 |
| **RTT** | **98-109** | **12** | **`GUCUUUCCUGGU`** | templates c.87-c.98 |
| **PBS** | **110-121** | **12** | **`GCUUCUGCCACA`** | anneals c.75-c.86 |
| tevopreQ1 | 122-158 | 37 | `CGCGGUUCUAUCUAGUUACGCGUUAAACCAACUAGAA` | 3' protective pseudoknot |

The reason for the confusion: position 109 (`U`) base-pairs with nothing
in the primer duplex. It is the **first templating nucleotide**, the one
that directs incorporation of A opposite c.87, so it belongs to the RTT.
Both splits encode the same c.88G>C edit, because only the boundary
label was wrong, not the sequence. The construct is correct; the
annotation was not.

Note that the total is unchanged at 24 nt, which is why no length check
caught this. This is the same failure mode as the two compensating
sequence errors found earlier, and the same remedy applies: the
validator must reconstruct positions against the reference, not sum
lengths.

## Design assessment

| Parameter | Value | Assessment |
|---|---|---|
| PBS length | 12 nt | within the usual 8-15 nt range, GC 58% |
| Edit position | +2 from the nick | favourable, PE is most efficient for edits close to the nick |
| 3' homology after the edit | 10 nt | at the low end; longer 3' homology usually improves flap resolution and is worth raising with the requesting group |
| PAM disruption | none | the PAM is untouched |
| Re-nicking of the edited allele | reduced | the edit sits at protospacer position 20 of 21, i.e. 2 nt from the PAM; PAM-proximal mismatches strongly impair SpCas9 re-binding, so the edited allele is poorly re-cut. This is a favourable property of the site. |
| First nucleotide of the 3' extension | `G` (position 98) | correct; a `C` here is known to mispair with scaffold G81 |

## Experimentally validated junction positions

The requesting work already established that fragment ligation of a
plain sgRNA at positions **34/35** and **57/58** does not reduce
CRISPR/Cas9 activity. Because PegRNA3 shares the spacer and scaffold with
that sgRNA over positions 1-97, those coordinates transfer directly.

Comparison under the current best model (tevopreQ1 pseudoknot enforced,
spacer constrained unpaired), P(paired) for the two nucleotides on each
side of the cut:

| Junction | 5' side | 3' side | Status |
|---|---|---|---|
| 34/35 | 0.073 | 0.051 / 0.023 | **experimentally validated** |
| 35/36 | 0.051 | 0.023 / 0.011 | currently chosen |
| **57/58** | **0.059** | **0.140 / 0.117** | **experimentally validated** |
| 58/59 | 0.140 | 0.117 / 0.088 | — |
| 77/78 | 0.511 | 0.764 / 0.811 | currently chosen |
| 97/98 | 0.980 | 0.133 / 0.079 | blocks RT over-extension |
| 121/122 | 0.100 | 0.000 / 0.000 | isolates tevopreQ1 |

Two consequences:

1. **34/35 and 35/36 are equivalent** by prediction. Since 34/35 is
   experimentally validated and 35/36 is not, there is no reason to
   prefer 35/36. Using the validated coordinate costs nothing and gains
   an experimental precedent.
2. **57/58 is far better than 77/78** and is experimentally validated.
   Junction 77/78 sits with its 3' side inside a helix at P = 0.76-0.81,
   whereas 57/58 has P = 0.06 on the 5' side and 0.12-0.14 on the 3'
   side. Replacing 77/78 with 57/58 removes the only moderate-risk
   scaffold junction in the design and replaces it with a position that
   has already been shown not to impair Cas9.

Recommended junction set for the multi-fragment construct, pending
confirmation with the requesting group:

```
34/35  (validated)  +  57/58  (validated)  +  97/98  (blocks RT over-extension)
```

with 121/122 as an optional fourth junction if the 3' motif is to be
synthesised separately.

## Open items

- Confirm with the requesting group whether A30P installation is
  intended, or whether correction of A30P or of A53T (c.157G>A) was
  meant. c.157 is in codon 53 and is not reachable from this protospacer.
- Confirm whether the 3' homology of 10 nt is deliberate.
- Once confirmed, update `data/raw/designs/PegRNA3_fragments.tsv` and the
  analysis scripts to the 12 + 12 boundary and re-run everything.

## Sources

- SNCA mRNA NM_000345.4 and RefSeqGene NG_011851.1, NCBI Nucleotide.
  https://www.ncbi.nlm.nih.gov/nuccore/NM_000345.4
- Anzalone, A. V. et al. Search-and-replace genome editing without
  double-strand breaks or donor DNA. *Nature* **576**, 149-157 (2019).
  https://doi.org/10.1038/s41586-019-1711-4
- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. *Nature* **631**, 224-231 (2024).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
- Nelson, J. W. et al. Engineered pegRNAs improve prime editing
  efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
