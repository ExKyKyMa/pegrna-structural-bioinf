# PegRNA3 structure definition

## Completed

- Defined the construct architecture as S–F1–L1–F2–L2–E–V.
- Separated Evopreq1 from the pegRNA extension.
- Created a nucleotide-level unmodified reference sequence.
- Created a nucleotide-level modified reference sequence.
- Created a fragment table.
- Created a linker junction table.
- Marked PBS and RTT boundaries as unknown.
- Preserved the planned positions of Linker 1 and Linker 2.

## Important modeling decision

The modified and unmodified FASTA files are identical at the
nucleotide level. The chemical difference will be represented later
using explicit linker atoms and modified connectivity in 3D models.

## Current uncertainty

The exact fragment boundaries and reaction connectivity need to be
validated before building explicit linker structures.

## Next step

Run an automated sequence validation script to calculate fragment
lengths, total sequence length and region positions.