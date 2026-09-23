"""Analyse the base-pairing context around chemical linker junctions.

Usage:
    python scripts/analyze_junction_neighborhood.py 35 77 121
    python scripts/analyze_junction_neighborhood.py 35 97 --label two_linker
    python scripts/analyze_junction_neighborhood.py 35 77 121 --unpaired-spacer

A junction given as N means "the linker is inserted between nucleotide N
and nucleotide N+1" in the full-length PegRNA3 numbering (1-based).

The script:
  1. reads the single source of truth FASTA;
  2. folds the continuous (unmodified) sequence with ViennaRNA;
  3. reports, for every junction, the MFE pairing state and the
     ensemble base-pair probabilities of the nucleotides in a window
     around the junction;
  4. writes a TSV report per construct.
"""

import argparse
from pathlib import Path

import RNA
import pandas as pd


FASTA = Path("data/raw/sequences/PegRNA3_unmodified_reference.fasta")
OUTDIR = Path("results/secondary_structure/junctions")

REGIONS = [
    ("Spacer", 1, 21),
    ("Scaffold_fragment_1", 22, 35),
    ("Scaffold_fragment_2", 36, 97),
    ("Extension", 98, 121),
    ("Evopreq1", 122, 158),
]

WINDOW = 4  # nucleotides inspected on each side of a junction

SPACER_RANGE = (1, 21)  # positions forced unpaired by --unpaired-spacer


def read_fasta(path):
    lines = path.read_text().splitlines()
    seq = "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    )
    return seq.upper().replace("T", "U")


def region_of(position):
    for name, start, end in REGIONS:
        if start <= position <= end:
            return name
    return "Unknown"


def pair_table(structure):
    """Return {position: partner} for a dot-bracket string (1-based)."""
    stack, pairs = [], {}
    for index, symbol in enumerate(structure, start=1):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            opener = stack.pop()
            pairs[opener] = index
            pairs[index] = opener
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("junctions", nargs="+", type=int)
    parser.add_argument("--label", default=None)
    parser.add_argument(
        "--unpaired-spacer",
        action="store_true",
        help=(
            "Force the spacer (positions "
            f"{SPACER_RANGE[0]}-{SPACER_RANGE[1]}) to stay unpaired. "
            "In the assembled prime editor the spacer is hybridised to "
            "the target DNA strand, so its intramolecular pairing with "
            "the PBS region seen in free-RNA folding is an artifact."
        ),
    )
    args = parser.parse_args()

    junctions = sorted(args.junctions)
    label = args.label or f"{len(junctions)}_linker"

    if args.unpaired_spacer:
        label = f"{label}_unpaired_spacer"

    sequence = read_fasta(FASTA)
    print(f"Sequence: {len(sequence)} nt from {FASTA}")

    model = RNA.md()
    model.temperature = 37.0
    fold_compound = RNA.fold_compound(sequence, model)

    if args.unpaired_spacer:
        for position in range(SPACER_RANGE[0], SPACER_RANGE[1] + 1):
            fold_compound.hc_add_up(
                position,
                RNA.CONSTRAINT_CONTEXT_ALL_LOOPS,
            )
        print(
            f"Hard constraint: positions "
            f"{SPACER_RANGE[0]}-{SPACER_RANGE[1]} forced unpaired"
        )
    else:
        print("No constraints: free-RNA folding")

    mfe_structure, mfe = fold_compound.mfe()
    fold_compound.exp_params_rescale(mfe)
    _, ensemble_energy = fold_compound.pf()

    probabilities = fold_compound.bpp()
    pairs = pair_table(mfe_structure)

    print(f"MFE: {mfe:.2f} kcal/mol")
    print(f"Ensemble free energy: {ensemble_energy:.2f} kcal/mol")
    print(f"MFE structure:\n{mfe_structure}")

    rows = []

    for junction in junctions:
        low = max(1, junction - WINDOW + 1)
        high = min(len(sequence), junction + WINDOW)

        print(f"\n=== Junction {junction}/{junction + 1} "
              f"({region_of(junction)} -> {region_of(junction + 1)}) ===")
        print(f"local sequence {low}-{high}: {sequence[low - 1:high]}")

        for position in range(low, high + 1):
            partners = []
            for other in range(1, len(sequence) + 1):
                if other == position:
                    continue
                i, j = min(position, other), max(position, other)
                probability = probabilities[i][j]
                if probability >= 0.01:
                    partners.append((other, probability))

            partners.sort(key=lambda item: -item[1])
            paired_probability = sum(p for _, p in partners)
            mfe_partner = pairs.get(position, 0)

            rows.append({
                "construct": label,
                "spacer_constrained": bool(args.unpaired_spacer),
                "junction": f"{junction}/{junction + 1}",
                "position": position,
                "base": sequence[position - 1],
                "region": region_of(position),
                "side": "left" if position <= junction else "right",
                "mfe_partner": mfe_partner,
                "mfe_paired": bool(mfe_partner),
                "total_pairing_probability": round(paired_probability, 4),
                "top_partners": "; ".join(
                    f"{other}:{probability:.3f}"
                    for other, probability in partners[:3]
                ),
            })

            print(
                f"  {position:>3} {sequence[position - 1]} "
                f"{region_of(position):<20} "
                f"MFE_partner={mfe_partner:<4} "
                f"P(paired)={paired_probability:.3f}  "
                + "; ".join(
                    f"{other}:{probability:.3f}"
                    for other, probability in partners[:3]
                )
            )

    output = pd.DataFrame(rows)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out_file = OUTDIR / f"PegRNA3_{label}_junction_context.tsv"
    output.to_csv(out_file, sep="\t", index=False)
    print(f"\nSaved: {out_file}")


if __name__ == "__main__":
    main()
