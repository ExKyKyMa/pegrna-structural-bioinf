# PegRNA3

## Construct architecture

### Chemically assembled construct

```text
Spacer–Scaffold fragment 1–Linker 1–Scaffold fragment 2–Linker 2–Extension–Evopreq1
```

### Unmodified reference

```text
Spacer–Scaffold fragment 1–Scaffold fragment 2–Extension–Evopreq1
```

## Components

| Component | Sequence or description | Role |
|---|---|---|
| Spacer | `ACAGGGUGUGGCAGAAGC` | Target recognition |
| Scaffold fragment 1 | `UUUUAGAGCUAGA` | First scaffold fragment |
| Linker 1 | 3′-azide–propargylamine–C(O)– | Chemical junction 1 |
| Scaffold fragment 2 | `AAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC` | Second scaffold fragment |
| Linker 2 | 3′-amino–azidobutanoic-acid–propargylamine–C(O)– | Chemical junction 2 |
| Extension | `GUGCGUCUUUCCUGGUGCUUCUGCCACA` | PBS/RTT-containing region; boundaries unknown |
| Evopreq1 | `CGCGGUUCUAUCUAGUUACGCGUUAAACCAACUAGAA` | Protective 3′ motif against exonuclease degradation |

## PBS and RTT annotation

The exact PBS and RTT boundaries are currently unknown because the
original target sequence and pegRNA design parameters are unavailable.

For the initial structural analysis, the full extension is treated as
one region.

## Structural hypothesis

Linker 1 and Linker 2 are positioned in non-base-paired and
solvent-exposed regions of RNA hairpins.

The expected structure is:

- preserved scaffold stem geometry;
- preserved Extension structure;
- preserved Evopreq1 structure;
- localized structural differences only near the linker junctions;
- no major steric overlap with the prime-editing complex.

## Computational comparison

### Reference system

```text
S–F1–F2–E–V
```

### Modified system

```text
S–F1–L1–F2–L2–E–V
```

The nucleotide-level FASTA representations are identical because the
chemical linkers are not represented by standard RNA letters.

The chemically explicit difference will be introduced during 3D model
construction and topology preparation.

## Analysis plan

1. Validate the annotated sequence lengths.
2. Verify that all regions occur in the full sequence.
3. Predict the reference secondary structure.
4. Use a documented surrogate representation for linker junctions.
5. Compare secondary-structure ensembles.
6. Build local 3D models of both linker junctions.
7. Compare stem geometry and Evopreq1 structure.
8. Prepare MD only after chemical connectivity is verified.

## Initial folding observation

The nucleotide-level reference model predicts stable base-pairing near
both intended linker junctions.

For Junction 1, high-probability pairs include:

- 26–37: 0.8297
- 27–36: 0.8576
- 28–35: 0.8548
- 29–34: 0.8231

For Junction 2, the scaffold stem remains highly probable near the
junction:

- 79–93: 0.9962
- 80–92: 0.9987
- 81–91: 0.9984
- 82–90: 0.9988
- 83–89: 0.9986
- 84–88: 0.8739

Moderate-probability interactions between Scaffold_fragment_2 and the
Extension were also detected:

- 76–97: 0.7226
- 77–96: 0.7159

## Interpretation

The initial continuous-sequence model does not support the assumption
that both linker junctions are located in completely unpaired regions.

This result must be interpreted cautiously because standard nucleotide
folding does not explicitly represent the chemical linkers or their
connectivity.

## Revised computational plan

1. Use the continuous sequence as a baseline.
2. Build constrained surrogate models with linker junctions treated as
   disconnected from standard nucleotide pairing.
3. Build explicit 3D models of both chemical linker junctions.
4. Compare local geometry and structural compatibility.
5. Perform MD only after explicit chemical connectivity has been validated.
