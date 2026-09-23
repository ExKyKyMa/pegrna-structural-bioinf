#!/usr/bin/env python3
"""Fold PegRNA3 while treating the 3' tevopreQ1 motif as a folded pseudoknot.

ViennaRNA's energy model is restricted to nested (non-crossing) base pairs, so
a pseudoknot cannot be written down as a hard constraint directly. What we can
do instead is encode its *consequence*: nucleotides that are already engaged in
the pseudoknot are unavailable for pairing with the rest of the molecule. We
therefore mark those positions as unpaired and let ViennaRNA fold everything
else around them.

Two additional artifacts are handled the same way:
  * the spacer (1-21) is bound to target DNA in vivo (--unpaired-spacer),
  * nothing else in the motif is constrained, so genuine single-stranded loop
    nucleotides of the pseudoknot may still pair if the model wants them to.

Reference structure of tevopreQ1 (local numbering 1-37 = global 122-158):
  H-type pseudoknot of the preQ1 class I riboswitch aptamer family.
  stem 1  1-4  :  19-22
  stem 2 (register A, default) 11-16 : 31-36
  stem 2 (register B)          12-17 : 30-35
  stem 2 (register C)           6-10 : 33-37   (shortest loop 1, closest to the
                                                preQ1-I consensus topology, but
                                                does not involve G15)
Both registers contain the G15-C32 pair, which is consistent with the
published observation that the G15C point mutation abolishes pseudoknot
formation and the epegRNA benefit (Nelson et al. 2022).

Sources:
  Nelson, J. W. et al. Nat. Biotechnol. 40, 402-410 (2022).
    https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
  Kang, M. et al. / Roth, A. et al., reviewed in
    Structure and function of preQ1 riboswitches.
    https://pmc.ncbi.nlm.nih.gov/articles/PMC4177978/
"""
import argparse
import csv
import pathlib

import RNA

REPO = pathlib.Path(__file__).resolve().parents[1]
FASTA = REPO / "data/raw/sequences/PegRNA3_unmodified_reference.fasta"
OUTDIR = REPO / "results/secondary_structure/pseudoknot"

SPACER = (1, 21)
MOTIF_OFFSET = 121  # local motif position 1 == global position 122

STEM1 = [(1, 22), (2, 21), (3, 20), (4, 19)]
REGISTERS = {
    "A": [(11, 36), (12, 35), (13, 34), (14, 33), (15, 32), (16, 31)],
    "B": [(12, 35), (13, 34), (14, 33), (15, 32), (16, 31), (17, 30)],
    "C": [(6, 37), (7, 36), (8, 35), (9, 34), (10, 33)],
}
JUNCTIONS = {"35/36": 35, "77/78": 77, "97/98": 97, "121/122": 121}
WINDOW = 4


def read_sequence() -> str:
    return "".join(
        line.strip() for line in FASTA.read_text().splitlines()
        if line and not line.startswith(">")
    )


def pair_table(structure: str) -> dict:
    stack, pairs = [], {}
    for pos, char in enumerate(structure, start=1):
        if char == "(":
            stack.append(pos)
        elif char == ")":
            opened = stack.pop()
            pairs[opened] = pos
            pairs[pos] = opened
    return pairs


def fold(sequence: str, unpaired: set):
    md = RNA.md()
    md.temperature = 37.0
    fc = RNA.fold_compound(sequence, md)
    for pos in sorted(unpaired):
        fc.hc_add_up(pos, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    structure, energy = fc.mfe()
    fc.exp_params_rescale(energy)
    fc.pf()
    return fc, structure, energy


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--register", choices=sorted(REGISTERS), default="A",
                    help="pseudoknot stem 2 register (default: A)")
    ap.add_argument("--unpaired-spacer", action="store_true",
                    help="also force the spacer (1-21) unpaired, as it is bound "
                         "to target DNA in vivo")
    args = ap.parse_args()

    sequence = read_sequence()
    motif_pairs = STEM1 + REGISTERS[args.register]
    committed = sorted(
        {p + MOTIF_OFFSET for pair in motif_pairs for p in pair}
    )

    unpaired = set(committed)
    label = f"register{args.register}"
    if args.unpaired_spacer:
        unpaired |= set(range(SPACER[0], SPACER[1] + 1))
        label += "_unpaired_spacer"

    fc, structure, energy = fold(sequence, unpaired)
    ens = fc.pf()[1]
    pairs = pair_table(structure)
    bpp = fc.bpp()

    print(f"length                 : {len(sequence)} nt")
    print(f"pseudoknot register    : {args.register}")
    print(f"motif nt held in knot  : {len(committed)} "
          f"({committed[0]}-{committed[-1]})")
    print(f"spacer forced unpaired : {args.unpaired_spacer}")
    print(f"MFE                    : {energy:8.2f} kcal/mol")
    print(f"ensemble free energy    : {ens:8.2f} kcal/mol")
    print()
    print(sequence)
    print(structure)
    print()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, cut in JUNCTIONS.items():
        lo = max(1, cut - WINDOW + 1)
        hi = min(len(sequence), cut + WINDOW)
        for pos in range(lo, hi + 1):
            total = sum(
                bpp[min(pos, other)][max(pos, other)]
                for other in range(1, len(sequence) + 1) if other != pos
            )
            partner = pairs.get(pos)
            rows.append({
                "junction": name,
                "position": pos,
                "nucleotide": sequence[pos - 1],
                "side": "5prime" if pos <= cut else "3prime",
                "mfe_partner": partner if partner else "",
                "p_paired": round(total, 4),
                "pseudoknot_register": args.register,
                "spacer_constrained": args.unpaired_spacer,
            })

    out = OUTDIR / f"PegRNA3_pseudoknot_{label}_junction_context.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for name in JUNCTIONS:
        sel = [r for r in rows if r["junction"] == name]
        summary = " ".join(f"{r['position']}{r['nucleotide']}:{r['p_paired']:.3f}"
                           for r in sel)
        print(f"junction {name:>8} | {summary}")
    print(f"\nwritten: {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
