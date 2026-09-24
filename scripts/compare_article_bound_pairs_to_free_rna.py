from pathlib import Path

import pandas as pd
import RNA


FASTA_FILE = Path(
    "data/raw/sequences/article_pegrna_reference.fasta"
)

BOUND_PAIRS_FILE = Path(
    "data/raw/designs/article_pegrna_bound_rna_pairs.tsv"
)

OUTPUT_FILE = Path(
    "results/secondary_structure/article_pegrna_bound_pairs_in_free_ensemble.tsv"
)


def read_fasta(path):
    lines = path.read_text().splitlines()

    return "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    ).upper().replace("T", "U")


def main():
    sequence = read_fasta(FASTA_FILE)

    print(f"Article pegRNA length: {len(sequence)} nt")

    pairs = pd.read_csv(
        BOUND_PAIRS_FILE,
        sep="\t",
        dtype=str,
    )

    model = RNA.md()
    model.temperature = 37.0

    fold_compound = RNA.fold_compound(sequence, model)

    mfe_structure, mfe = fold_compound.mfe()

    fold_compound.exp_params_rescale(mfe)
    ensemble_free_energy = fold_compound.pf()[1]

    probabilities = fold_compound.bpp()

    results = []

    for _, row in pairs.iterrows():
        position_i = int(row["position_i"])
        position_j = int(row["position_j"])

        probability = probabilities[position_i][position_j]

        results.append({
            "pair_id": row["pair_id"],
            "position_i": position_i,
            "base_i": sequence[position_i - 1],
            "position_j": position_j,
            "base_j": sequence[position_j - 1],
            "element": row["region_or_element"],
            "free_rna_pair_probability": round(probability, 6),
            "bound_state_pair": True,
        })

    output = pd.DataFrame(results)
    output = output.sort_values(
        ["element", "position_i"],
    )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.to_csv(
        OUTPUT_FILE,
        sep="\t",
        index=False,
    )

    print(f"MFE: {mfe:.2f} kcal/mol")
    print(f"Ensemble free energy: {ensemble_free_energy:.2f} kcal/mol")
    print()
    print("Protein-bound RNA-RNA pairs in free-RNA ensemble:")
    print(output.to_string(index=False))

    print()
    print("Mean probability by structural element:")

    summary = (
        output
        .groupby("element")["free_rna_pair_probability"]
        .agg(["mean", "min", "max", "count"])
        .round(4)
    )

    print(summary)


if __name__ == "__main__":
    main()