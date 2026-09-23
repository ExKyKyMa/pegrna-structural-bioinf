# PegRNA3 secondary-structure ensemble sampling

## Method

- Software: ViennaRNA RNAsubopt
- Sampling method: Boltzmann stochastic backtracking
- Number of sampled structures: 10,000
- Input: PegRNA3 nucleotide-level reference sequence
- Architecture: S-F1-F2-E-V

## Results

- Most frequent sampled structure: MFE structure
- Frequency of most frequent structure: 0.0056
- MFE frequency from RNAfold: 0.0060157
- Ensemble diversity from RNAfold: 22.59
- MFE–ensemble free-energy gap: 3.15 kcal/mol

## Interpretation

The MFE structure is also the most frequently sampled individual
structure, but its estimated ensemble frequency is only approximately
0.6%.

The first ten sampled structures account for approximately 3.76% of
the sampled ensemble. Therefore, no single secondary structure
dominates the ensemble.

The major structural framework appears to be partially conserved,
while local alternative base-pairing patterns occur in the internal
and terminal regions.

## Consequence for further modeling

The MFE structure will be retained as one reference conformer, but
the subsequent analysis will use ensemble-derived base-pair
probabilities and representative structural clusters rather than
assuming that the MFE structure is the only relevant conformation.

## Next step

Calculate and visualize base-pair probabilities by region and identify
structural clusters relevant to the two linker junctions.