# PegRNA3 sequence validation

## Validation result

The PegRNA3 nucleotide sequence was successfully validated.

| Region | Length | Start | End |
|---|---:|---:|---:|
| Spacer | 21 | 1 | 21 |
| Scaffold fragment 1 | 14 | 22 | 35 |
| Scaffold fragment 2 | 62 | 36 | 97 |
| Extension | 24 | 98 | 121 |
| Evopreq1 | 37 | 122 | 158 |

## Total length

- Total assembled sequence length: 158 nt

## Architecture

```text
S(1–21)–F1(22–35)–F2(36–97)–E(98–121)–V(122–158)
```

## Interpretation

The nucleotide-level reference sequence is internally consistent.
The chemical linker junctions are not represented as nucleotide symbols
and will be introduced later in the explicit 3D model.

## Current limitation

The exact PBS and RTT boundaries remain unknown. The full extension
region is therefore treated as one annotated region for the initial
structural analysis.

## Next step

Run RNAfold for the unmodified nucleotide-level reference sequence.

## RNAfold output

- Sequence length: 158 nt
- MFE structure energy: -54.30 kcal/mol
- Ensemble free energy: -57.45 kcal/mol
- Centroid structure energy: -51.30 kcal/mol
- Centroid-to-ensemble base-pair distance: 15.60

## Interpretation

The MFE structure is not expected to represent the entire structural
ensemble. The ensemble free energy is more favorable than the MFE
structure energy because it includes contributions from alternative
secondary structures.

The centroid structure represents a structure that is close to the
Boltzmann-weighted ensemble by base-pair distance. The value d=15.60
is a base-pair distance measure, not a physical distance in angstroms
or nanometers.

These global values are baseline descriptors. Local analysis of the
scaffold, linker junctions, extension and Evopreq1 is still required.