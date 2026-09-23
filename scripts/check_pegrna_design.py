from pathlib import Path

import pandas as pd


INPUT_FILE = Path("data/raw/designs/PegRNA3_fragments.tsv")
REGIONS_FILE = Path("data/raw/designs/PegRNA3_regions.tsv")
REFERENCE_FASTA = Path("data/raw/sequences/PegRNA3_unmodified_reference.fasta")
OUTPUT_FILE = Path("data/processed/PegRNA3_sequence_validation.tsv")


EXPECTED_ORDER = [
    "Spacer",
    "Scaffold_fragment_1",
    "Scaffold_fragment_2",
    "Extension",
    "Evopreq1",
]

ALLOWED_SYMBOLS = set("AUGCN")


def read_reference():
    """Read the single source of truth for the PegRNA3 sequence."""
    lines = REFERENCE_FASTA.read_text().splitlines()
    return clean_sequence("".join(l for l in lines if l and not l.startswith(">")))


def check_against_reference(assembled):
    """Compare the assembled fragments with the reference character by character.

    Summing fragment lengths is not enough: two compensating errors that keep
    the total length unchanged pass a length check but change every downstream
    coordinate. This comparison is the check that would have caught them.
    """
    reference = read_reference()
    if assembled == reference:
        print(f"Reference check: OK, {len(assembled)} nt identical to "
              f"{REFERENCE_FASTA}")
        return True

    print("Reference check: FAILED")
    print(f"  assembled length {len(assembled)} nt, "
          f"reference length {len(reference)} nt")
    for i, (a, b) in enumerate(zip(assembled, reference), start=1):
        if a != b:
            print(f"  first difference at position {i}: "
                  f"assembled {a}, reference {b}")
            lo, hi = max(0, i - 11), i + 10
            print(f"  assembled {lo + 1}-{hi}: {assembled[lo:hi]}")
            print(f"  reference {lo + 1}-{hi}: {reference[lo:hi]}")
            break
    return False


def check_regions():
    """Check that every functional region matches the reference at its stated
    coordinates. This is what verifies the RTT / PBS boundary."""
    if not REGIONS_FILE.exists():
        print(f"Region check: skipped, {REGIONS_FILE} not found")
        return True

    reference = read_reference()
    regions = pd.read_csv(REGIONS_FILE, sep="\t", dtype=str).fillna("")
    ok = True
    print("Region check:")
    for _, row in regions.iterrows():
        start, end = int(row["start_nt"]), int(row["end_nt"])
        expected = clean_sequence(row["sequence"])
        actual = reference[start - 1:end]
        stated_length = int(row["length_nt"])
        problems = []
        if actual != expected:
            problems.append(f"sequence mismatch, reference has {actual}")
        if len(expected) != stated_length:
            problems.append(f"stated length {stated_length} != {len(expected)}")
        if end - start + 1 != stated_length:
            problems.append("coordinates do not match stated length")
        status = "OK" if not problems else "FAILED: " + "; ".join(problems)
        print(f"  {row['region']:<20} {start:>4}-{end:<4} "
              f"{stated_length:>3} nt  {status}")
        ok = ok and not problems

    covered = sorted(
        (int(r["start_nt"]), int(r["end_nt"])) for _, r in regions.iterrows()
    )
    position = 1
    for start, end in covered:
        if start != position:
            print(f"  gap or overlap in coverage at position {position}")
            ok = False
        position = end + 1
    if position - 1 != len(reference):
        print(f"  regions cover {position - 1} nt, reference is "
              f"{len(reference)} nt")
        ok = False
    return ok


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

    print()
    reference_ok = check_against_reference(assembled_sequence)
    print()
    regions_ok = check_regions()

    if not (reference_ok and regions_ok):
        raise SystemExit(1)


if __name__ == "__main__":
    main()