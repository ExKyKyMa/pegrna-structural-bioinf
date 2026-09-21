from pathlib import Path

import RNA
import pandas as pd


INPUT_FASTA = Path(
    "data/raw/sequences/PegRNA3_unmodified_reference.fasta"
)

OUTPUT_TABLE = Path(
    "results/secondary_structure/PegRNA3_base_pair_probabilities.tsv"
)

REGIONS = [
    ("Spacer", 1, 18),
    ("Scaffold_fragment_1", 19, 31),
    ("Scaffold_fragment_2", 32, 93),
    ("Extension", 94, 121),
    ("Evopreq1", 122, 158),
]

JUNCTION_REGIONS = {
    "Junction_1": (25, 38),
    "Junction_2": (87, 100),
}


def read_fasta(path):
    lines = path.read_text().splitlines()

    sequence = "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    )

    return sequence.replace("T", "U").upper()


def region_for_position(position):
    for name, start, end in REGIONS:
        if start <= position <= end:
            return name

    return "Unknown"


def main():
    sequence = read_fasta(INPUT_FASTA)

    print(f"Sequence length: {len(sequence)} nt")

    if len(sequence) != 158:
        raise ValueError(
            f"Expected 158 nt, but found {len(sequence)} nt"
        )

    model = RNA.md()
    model.temperature = 37.0

    fold_compound = RNA.fold_compound(sequence, model)

    mfe_structure, mfe = fold_compound.mfe()

    fold_compound.exp_params_rescale(mfe)
    fold_compound.pf()

    probabilities = fold_compound.bpp()

    rows = []

    for position_i in range(1, len(sequence) + 1):
        for position_j in range(
            position_i + 1,
            len(sequence) + 1,
        ):
            probability = probabilities[position_i][position_j]

            if probability < 0.01:
                continue

            rows.append({
                "position_i": position_i,
                "base_i": sequence[position_i - 1],
                "region_i": region_for_position(position_i),
                "position_j": position_j,
                "base_j": sequence[position_j - 1],
                "region_j": region_for_position(position_j),
                "pair_probability": round(probability, 6),
                "distance_nt": position_j - position_i,
            })

    output = pd.DataFrame(rows)

    if not output.empty:
        output = output.sort_values(
            "pair_probability",
            ascending=False,
        )

    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT_TABLE, sep="\t", index=False)

    print(f"MFE: {mfe:.2f} kcal/mol")
    print(f"Pairs with probability >= 0.01: {len(output)}")
    print(f"Saved: {OUTPUT_TABLE}")

    print("\nTop 20 base pairs:")
    print(output.head(20).to_string(index=False))

    print("\nHigh-confidence pairs by region:")

    high_confidence = output[
        output["pair_probability"] >= 0.80
    ]

    if high_confidence.empty:
        print("No pairs with probability >= 0.80")
    else:
        region_counts = (
            high_confidence
            .groupby(["region_i", "region_j"])
            .size()
            .sort_values(ascending=False)
        )

        print(region_counts)

    print("\nPairs near linker junctions:")

    for junction_name, (start, end) in JUNCTION_REGIONS.items():
        nearby = output[
            (
                output["position_i"].between(start, end)
                | output["position_j"].between(start, end)
            )
        ].sort_values(
            "pair_probability",
            ascending=False,
        )

        print()
        print(
            f"{junction_name}: "
            f"positions {start}-{end}"
        )

        if nearby.empty:
            print("No pairs with probability >= 0.01")
        else:
            print(nearby.to_string(index=False))


if __name__ == "__main__":
    main()