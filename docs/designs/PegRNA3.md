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