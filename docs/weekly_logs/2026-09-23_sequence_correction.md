# Correction report — PegRNA3 reference sequence (2026-09-23)

## Summary

An automated comparison between the README synthesis design and the
`data/` + `results/` files revealed **one systematic error** affecting every
nucleotide-level reference file, the RNAfold input and all downstream
secondary-structure results.

The error is a **4-nt block that was deleted from the spacer/scaffold boundary
and re-inserted at the scaffold/extension boundary**. Total length stayed 158 nt,
so a length check could not detect it.

## The error

| | Spacer | Scaffold fragment 1 | Extension |
|---|---|---|---|
| Files (wrong) | `ACAGGGUGUGGCAGAAGC` (18 nt) | `UUUUAGAGCUAGA` (13 nt) | `GUGCGUCUUUCC...` (28 nt) |
| README / correct | `ACAGGGUGUGGCAGAAGCAGC` (21 nt) | `GUUUUAGAGCUAGA` (14 nt) | `GUCUUUCC...` (24 nt) |

- `AGC` was dropped after position 18 and the scaffold-opening `G` was dropped,
  i.e. the block `AGCG` (canonical positions 19–22) is missing at the 5' end.
- The same four letters reappear as a spurious `GUGC` at the start of the
  extension (files positions 94–97).

## Why the README version is correct

1. **Canonical SpCas9 scaffold.** The corrected sequence contains the full
   76-nt SpCas9 sgRNA scaffold `GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCC
   GUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC` exactly at positions 22–97.
   The version in the files **does not contain it** — the scaffold is
   truncated at its 5' end, which would break the repeat:anti-repeat duplex
   that binds Cas9.
2. **Fragment reassembly.** Both synthesis sets from the README
   (PegRNA31/32/33 and PegRNA41/42/43/44) reassemble byte-for-byte into the
   corrected 158-nt sequence. They do **not** reassemble into the file version.
3. **PBS complementarity.** In the corrected sequence the 3' end of the
   extension (`UGCUUCUGCCACA`) is the reverse complement of spacer positions
   9–21 — exactly what a primer binding site must be. The spurious `GUGC`
   in the file version sits inside the RT template and has no such rationale.
4. **Junction coordinates.** The corrected numbering reproduces the agreed
   linker positions 35/36, 77/78, 97/98, 121/122 as clean fragment boundaries.

## Corrected region map (158 nt)

| Region | Positions | Length |
|---|---|---|
| Spacer | 1–21 | 21 |
| Scaffold fragment 1 | 22–35 | 14 |
| Scaffold fragment 2 | 36–97 | 62 |
| RTT (RT template) | 98–108 | 11 |
| PBS | 109–121 | 13 |
| evopreQ1 motif | 122–158 | 37 |

Previously PBS and RTT were marked `unknown`; they are now assigned from
PBS/spacer reverse complementarity and should still be confirmed against the
original synthesis design.

## Impact on existing results — must be recomputed

All secondary-structure work used the wrong sequence and is **invalid**:

- `PegRNA3_unmodified_reference.fold` (MFE −52.00, ensemble −54.91, centroid −48.80)
- `PegRNA3_10000_samples.txt` and the sampling statistics
- `PegRNA3_base_pair_probabilities.tsv` — its region labels are also wrong:
  e.g. pairs reported as "Spacer 11 G — Extension 117 C" are shifted by the
  4-nt displacement and the "Extension" assignment reflects the corrupted
  boundary.
- `PegRNA3_unmodified_reference_ss.svg` / `_ss.ps` / `_dp.ps`

## Other issues found

- `components.tsv` still contains only placeholders; no SMILES have been added,
  although linker SMILES and xTB geometries exist elsewhere in the project.
- `pegrna_designs.tsv` has `full_sequence_length = unknown` — now resolvable (158).
- `PegRNA3_fragments.tsv` had `length = unknown` for every fragment.
- The fragment table described only the two-linker architecture; the
  three-linker (4-fragment) design from the README was not represented.
- The modified and unmodified FASTA files were byte-identical, which is the
  intended modeling decision, but there was no machine-readable record of
  linker positions tying them together.

## Corrected files

`PegRNA3_full_sequence.txt`, `PegRNA3_unmodified_reference.fasta`,
`PegRNA3_MFE_input.txt`, `PegRNA3_regions.tsv`, `PegRNA3_fragments.tsv`,
`PegRNA3_junctions.tsv`, `validate_pegrna3.py`.

`validate_pegrna3.py` re-checks length, alphabet, scaffold identity, all region
and fragment coordinates, reassembly of both constructs, FASTA/MFE consistency
and PBS complementarity. Current status: 34/34 checks passed.
