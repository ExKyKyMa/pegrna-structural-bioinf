from collections import Counter
from pathlib import Path


INPUT_FILE = Path(
    "results/secondary_structure/sampling/PegRNA3_10000_samples.txt"
)


def is_dot_bracket(line):
    allowed = set(".()[]{}<>,|")
    structure = line.split()[0]

    return (
        len(structure) > 0
        and all(symbol in allowed for symbol in structure)
        and any(symbol in structure for symbol in "()")
    )


def main():
    structures = []

    with INPUT_FILE.open() as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                continue

            if is_dot_bracket(line):
                structure = line.split()[0]
                structures.append(structure)

    counts = Counter(structures)

    if not counts:
        print("No structures were found.")
        print("Please inspect the RNAsubopt output format.")
        return

    total = len(structures)

    print(f"Structures read: {total}")
    print()
    print("Most frequent structures:")

    for structure, count in counts.most_common(10):
        frequency = count / total
        print(f"{frequency:.4f}\t{count}\t{structure}")


if __name__ == "__main__":
    main()