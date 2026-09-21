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
| Spacer | 1 | 18 |
| Scaffold fragment 1 | 19 | 31 |
| Scaffold fragment 2 | 32 | 93 |
| Extension | 94 | 121 |
| Evopreq1 | 122 | 158 |

## Interpretation rules

- Probability >= 0.80: high-confidence predicted pair
- Probability 0.50–0.80: moderate-confidence pair
- Probability 0.10–0.50: alternative or weakly supported pair
- Probability < 0.10: not treated as a stable pair

## Results

- MFE: -52.00 kcal/mol
- Ensemble free energy: -54.91 kcal/mol
- MFE frequency: 0.0089435
- Ensemble diversity: 21.92

## Main observations

- Scaffold high-confidence pairs: pending
- Evopreq1 high-confidence pairs: pending
- Extension–Evopreq1 pairs: pending
- Potentially unpaired linker junction 1: positions 31/32
- Potentially unpaired linker junction 2: positions 93/94

## Next step

Visualize base-pair probabilities and annotate the secondary structure
by functional region.