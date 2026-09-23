# PegRNA3 base-pair probability analysis

## Input

- Sequence: PegRNA3 unmodified reference
- Length: 158 nt
- Software: ViennaRNA Python API
- Temperature: 37 °C
- Probability cutoff: 0.01

## Regions

| Region | Start | End |
|---|---:|---:|
| Spacer | 1 | 21 |
| Scaffold fragment 1 | 22 | 35 |
| Scaffold fragment 2 | 36 | 97 |
| Extension | 98 | 121 |
| Evopreq1 | 122 | 158 |

## Interpretation rules

- Probability >= 0.80: high-confidence predicted pair
- Probability 0.50–0.80: moderate-confidence pair
- Probability 0.10–0.50: alternative or weakly supported pair
- Probability < 0.10: not treated as a stable pair

## Results

- MFE: -54.30 kcal/mol
- Ensemble free energy: -57.45 kcal/mol
- MFE frequency: 0.0060157
- Ensemble diversity: 22.59

## Main observations

- High-confidence (p >= 0.80) pairs by region pair:
  - Spacer–Extension: 14
  - Evopreq1–Evopreq1: 7
  - Scaffold fragment 2–Scaffold fragment 2: 6
  - Scaffold fragment 1–Scaffold fragment 2: 4
  - Spacer–Evopreq1: 2
- Linker junction 1: positions 35/36. Nucleotide 35 itself is essentially
  unpaired, but the flanking pairs 32-39 (0.9469) and 33-38 (0.9119)
  close a hairpin across the junction.
- Linker junction 2: positions 97/98. Nucleotide 97 is paired with
  position 83 (0.9877) inside the scaffold stem; the first extension
  nucleotides (98-101) show only moderate pairing (<= 0.38).

## Next step

Visualize base-pair probabilities and annotate the secondary structure
by functional region.