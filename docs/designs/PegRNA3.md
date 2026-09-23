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
| Spacer | `ACAGGGUGUGGCAGAAGCAGC` | Target recognition |
| Scaffold fragment 1 | `GUUUUAGAGCUAGA` | First scaffold fragment |
| Linker 1 | 3′-azide–propargylamine–C(O)– | Chemical junction 1 |
| Scaffold fragment 2 | `AAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC` | Second scaffold fragment |
| Linker 2 | 3′-amino–azidobutanoic-acid–propargylamine–C(O)– | Chemical junction 2 |
| Extension | `GUCUUUCCUGGUGCUUCUGCCACA` | RTT (98–109) + PBS (110–121) |
| Evopreq1 | `CGCGGUUCUAUCUAGUUACGCGUUAAACCAACUAGAA` | Protective 3′ motif against exonuclease degradation |

## PBS and RTT annotation

Resolved 2026-09-23. The boundary is 12 + 12:

| Region | Positions | Length | Sequence | Role |
|---|---|---:|---|---|
| RTT | 98–109 | 12 | `GUCUUUCCUGGU` | templates SNCA c.87–c.98 |
| PBS | 110–121 | 12 | `GCUUCUGCCACA` | anneals to SNCA c.75–c.86 |

An earlier working note recorded 11 + 13 (RTT 98–108, PBS 109–121). That
is off by one: position 109 (`U`) pairs with nothing in the primer
duplex, it is the first templating nucleotide. A 13-nt PBS would require
an SpCas9 nick 2 nt from the PAM instead of the canonical 3 nt. Both
splits encode the same edit, so the construct is correct and only the
annotation was wrong.

The target site and the intended edit (SNCA, protospacer c.69–c.89, PAM
`AGG` at c.90–c.92, nick between c.86 and c.87, edit c.88G>C =
p.Ala30Pro) are derived and verified in
[PegRNA3_target_and_edit_definition.md](PegRNA3_target_and_edit_definition.md)
and recorded in `data/raw/designs/PegRNA3_regions.tsv` and
`data/raw/designs/PegRNA3_target_edit.tsv`.

Run `python3 scripts/check_pegrna_design.py` to verify that every region
matches the reference FASTA at its stated coordinates.

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

The nucleotide-level reference model (158 nt, MFE = -54.30 kcal/mol)
predicts base-pairing near both intended linker junctions.

Region coordinates in the full sequence:

| Region | Start | End | Length |
|---|---:|---:|---:|
| Spacer | 1 | 21 | 21 |
| Scaffold fragment 1 | 22 | 35 | 14 |
| Scaffold fragment 2 | 36 | 97 | 62 |
| Extension | 98 | 121 | 24 |
| Evopreq1 | 122 | 158 | 37 |

Junction 1 is located between positions 35 and 36.
Junction 2 is located between positions 97 and 98.

For Junction 1, high-probability pairs include:

- 32-39: 0.9469
- 33-38: 0.9119

For Junction 2, the scaffold stem remains highly probable near the
junction:

- 83-97: 0.9877
- 84-96: 0.9968
- 85-95: 0.9966
- 86-94: 0.9970

Moderate-probability interactions between Scaffold_fragment_1 and the
Extension were also detected:

- 27-101: 0.3658
- 28-100: 0.3762

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
