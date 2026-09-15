from pathlib import Path

import pandas as pd


INPUT_FILE = Path("data/raw/designs/PegRNA3_fragments.tsv")
OUTPUT_FILE = Path("data/processed/PegRNA3_sequence_validation.tsv")


EXPECTED_ORDER = [
    "Spacer",
    "Scaffold_fragment_1",
    "Scaffold_fragment_2",
    "Extension",
    "Evopreq1",
]

ALLOWED_SYMBOLS = set("AUGCN")


def clean_sequence(sequence):
    return "".join(str(sequence).split()).upper().replace("T", "U")


def find_column(df, expected_name):
    if expected_name in df.columns:
        return expected_name

    normalized = {
        column.lower().replace(" ", "_"): column
        for column in df.columns
    }

    normalized_name = expected_name.lower().replace(" ", "_")

    if normalized_name in normalized:
        return normalized[normalized_name]

    raise KeyError(
        f"Column '{expected_name}' was not found. "
        f"Available columns: {list(df.columns)}"
    )


def main():
    df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str).fillna("")

    fragment_column = find_column(df, "fragment_name")
    sequence_column = find_column(df, "sequence")

    rows = []
    assembled_parts = []

    for expected_region in EXPECTED_ORDER:
        matching_rows = df[df[fragment_column] == expected_region]

        if len(matching_rows) == 0:
            rows.append({
                "region": expected_region,
                "status": "missing_from_table",
                "length_nt": "",
                "sequence": "",
                "invalid_symbols": "",
                "start_nt": "",
                "end_nt": "",
            })
            continue

        if len(matching_rows) > 1:
            print(f"Warning: multiple rows found for {expected_region}")

        row = matching_rows.iloc[0]
        sequence = clean_sequence(row[sequence_column])
        invalid_symbols = sorted(set(sequence) - ALLOWED_SYMBOLS)

        start_nt = sum(len(part) for part in assembled_parts) + 1
        end_nt = start_nt + len(sequence) - 1

        rows.append({
            "region": expected_region,
            "status": "valid" if not invalid_symbols else "invalid_symbols",
            "length_nt": len(sequence),
            "sequence": sequence,
            "invalid_symbols": ",".join(invalid_symbols),
            "start_nt": start_nt,
            "end_nt": end_nt,
        })

        assembled_parts.append(sequence)

    assembled_sequence = "".join(assembled_parts)

    output_df = pd.DataFrame(rows)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("Region validation:")
    print(output_df.to_string(index=False))

    print()
    print(f"Total assembled sequence length: {len(assembled_sequence)} nt")
    print(f"Saved validation table to: {OUTPUT_FILE}")

    print()
    print("Assembled sequence:")
    print(assembled_sequence)


if __name__ == "__main__":
    main()