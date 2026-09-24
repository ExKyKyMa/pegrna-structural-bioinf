"""Where to place destabilising modifications in PegRNA3.

The idea
--------
PegRNA3 has a high Folding Barrier because the spacer is sequestered in a
strong intramolecular duplex with the PBS. A modification that locally
prevents or weakens base pairing would destabilise that duplex, lower the
barrier, and make the RNA easier for Cas9 to load. The question is where
to put such a modification so that it helps folding without breaking
anything the molecule needs in order to work.

This script answers that question systematically:

  1. It defines a protected set: every position the molecule cannot afford
     to lose. Modifying these would be self-defeating.
  2. It scans every remaining position one at a time, recomputing the
     Folding Barrier and the spacer accessibility.
  3. It searches combinations inside the safest zone and reports the best.
  4. It checks whether the effect combines with the triazole linkers.

How a modification is modelled
------------------------------
As a hard constraint that the position cannot pair. This is exact for an
abasic site, which has no nucleobase at all, and approximate for an
unlocked nucleic acid (UNA) residue, which still pairs but is strongly
destabilised: +4.0 to +6.6 kcal/mol for an internal UNA, with effects
that are 90-98% additive across multiple substitutions
(https://pmc.ncbi.nlm.nih.gov/articles/PMC2965255/). Treating UNA as
fully non-pairing is therefore an upper bound on its effect, and the
script says so rather than pretending otherwise.

Run from the repository root:

    python3 scripts/design_destabilising_modifications.py

Sources
-------
- Pasternak, A. & Wengel, J. Thermodynamics of RNA duplexes modified with
  unlocked nucleic acid nucleotides. Nucleic Acids Res. 38, 6697-6706
  (2010). https://pmc.ncbi.nlm.nih.gov/articles/PMC2965255/
- Yin, H. et al. Structure-guided chemical modification of guide RNA
  enables potent non-viral in vivo genome editing. Nat. Biotechnol. 35,
  1179-1187 (2017). https://pmc.ncbi.nlm.nih.gov/articles/PMC5901668/
- Ponnienselvan, K. et al. Reducing the inherent auto-inhibitory
  interaction within the pegRNA enhances prime editing efficiency.
  Nucleic Acids Res. 51, 6966-6980 (2023).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10359601/
- Zalatan, J. G. et al. Guide RNA structure design enables combinatorial
  CRISPRa programs. Nat. Commun. 15, 6341 (2024).
  https://www.nature.com/articles/s41467-024-50528-1
- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. Nature 631, 224-231 (2024).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
"""

from itertools import combinations
from pathlib import Path

import pandas as pd
import RNA

ROOT = Path(__file__).resolve().parent.parent

REFERENCE_FASTA = ROOT / "data/raw/sequences/PegRNA3_unmodified_reference.fasta"
ARTICLE_FASTA = ROOT / "data/raw/sequences/article_pegrna_reference.fasta"
BOUND_PAIRS_TSV = ROOT / "data/raw/designs/article_pegrna_bound_rna_pairs.tsv"
OUTPUT_DIR = ROOT / "results/modification_design"

SCAFFOLD_ANCHOR = "GUUUUAGAGCUAGA"

SPACER_RANGE = (1, 21)
SCAFFOLD_RANGE = (22, 97)
RTT_RANGE = (98, 109)
PBS_RANGE = (110, 121)
MOTIF_RANGE = (122, 158)
SEED_RANGE = (19, 21)

# The priming region: the 5' end of the PBS, immediately adjacent to the
# RTT. This is where the 3'-OH of the nicked DNA strand sits and where
# reverse transcription starts. A modification here would attack the
# chemistry the construct exists to perform.
PRIMING_REGION = (110, 115)

# The distal end of the PBS. This is the part that PBS-shortening removes,
# and shortening is a published, functionally tolerated intervention, so
# weakening it is the conservative place to start.
DISTAL_PBS = (116, 121)

LINKER_JUNCTIONS = [35, 36, 77, 78, 121, 122]

SEARCH_WIDTH = 400
CANONICAL_PAIRS = {"AU", "UA", "GC", "CG", "GU", "UG"}
MAX_COMBINATION_SIZE = 4


def read_fasta(path):
    return "".join(
        line.strip() for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ).upper().replace("T", "U")


def model():
    md = RNA.md()
    md.temperature = 37.0
    return md


def region_of(position):
    for name, (start, end) in [
        ("spacer", SPACER_RANGE), ("scaffold", SCAFFOLD_RANGE),
        ("RTT", RTT_RANGE), ("PBS", PBS_RANGE), ("motif", MOTIF_RANGE),
    ]:
        if start <= position <= end:
            return name
    return "outside"


def load_cryoem(sequence, article_sequence):
    """Return (usable pairs, every contacted position)."""
    offset = (sequence.index(SCAFFOLD_ANCHOR)
              - article_sequence.index(SCAFFOLD_ANCHOR))
    table = pd.read_csv(BOUND_PAIRS_TSV, sep="\t", dtype=str)

    usable, contacted = [], set()
    for _, row in table.iterrows():
        i = int(row["position_i"]) + offset
        j = int(row["position_j"]) + offset
        contacted.update((i, j))
        if (sequence[i - 1] + sequence[j - 1] in CANONICAL_PAIRS
                and j - i - 1 >= 3):
            usable.append((i, j))

    return usable, contacted


def protected_positions(contacted):
    """Positions that must not be modified, each with a stated reason.

    Being explicit about this is the point of the exercise. A modification
    that lowers the folding barrier by breaking the spacer or the priming
    site would look good in the calculation and fail in the experiment.
    """
    reasons = {}

    for position in contacted:
        reasons[position] = "contacts Cas9 in the cryo-EM structure"

    for position in range(SPACER_RANGE[0], SPACER_RANGE[1] + 1):
        reasons[position] = ("spacer, must base pair with the DNA target "
                             "for recognition")

    for position in range(RTT_RANGE[0], RTT_RANGE[1] + 1):
        reasons[position] = "RTT, is the template that is copied"

    for position in range(PRIMING_REGION[0], PRIMING_REGION[1] + 1):
        reasons[position] = ("PBS priming region, the DNA 3'-OH anneals "
                             "here and RT initiates")

    return reasons


def metrics(sequence, usable_pairs, blocked=()):
    """Folding Barrier and accessibility with the given positions blocked."""
    blocked = list(blocked)
    n = len(sequence)

    ground = RNA.fold_compound(sequence, model())
    for position in blocked:
        ground.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    ground_structure, ground_energy = ground.mfe()

    active = RNA.fold_compound(sequence, model())
    for i, j in usable_pairs:
        active.hc_add_bp(i, j,
                         RNA.CONSTRAINT_CONTEXT_ALL_LOOPS
                         | RNA.CONSTRAINT_CONTEXT_ENFORCE)
    for position in range(SPACER_RANGE[0], SPACER_RANGE[1] + 1):
        active.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    for position in blocked:
        active.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    active_structure, _ = active.mfe()

    saddle = RNA.fold_compound(sequence, model()).path_findpath_saddle(
        ground_structure, active_structure, SEARCH_WIDTH) / 100.0

    partition = RNA.fold_compound(sequence, model())
    for position in blocked:
        partition.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    _, rescale = partition.mfe()
    partition.exp_params_rescale(rescale)
    partition.pf()
    bpp = partition.bpp()

    def accessibility(start, end):
        values = [
            1.0 - sum(bpp[min(p, o)][max(p, o)]
                      for o in range(1, n + 1) if o != p)
            for p in range(start, end + 1)
        ]
        return sum(values) / len(values)

    return {
        "barrier": round(saddle - ground_energy, 2),
        "ground_energy": round(ground_energy, 2),
        "spacer_accessibility": round(accessibility(*SPACER_RANGE), 4),
        "seed_accessibility": round(accessibility(*SEED_RANGE), 4),
    }


def main():
    sequence = read_fasta(REFERENCE_FASTA)
    article_sequence = read_fasta(ARTICLE_FASTA)
    usable_pairs, contacted = load_cryoem(sequence, article_sequence)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    baseline = metrics(sequence, usable_pairs)
    sgrna = sequence[:SCAFFOLD_RANGE[1]]
    sgrna_pairs, _ = load_cryoem(sgrna, article_sequence)
    sgrna_reference = metrics(sgrna, sgrna_pairs)

    print(f"PegRNA3, {len(sequence)} nt")
    print(f"unmodified : barrier {baseline['barrier']}, "
          f"spacer accessibility {baseline['spacer_accessibility']}")
    print(f"sgRNA target: barrier {sgrna_reference['barrier']}, "
          f"spacer accessibility {sgrna_reference['spacer_accessibility']}")
    print()

    # ------------------------------------------------ 1. what is off limits
    reasons = protected_positions(contacted)
    protected = pd.DataFrame(
        [{"position": p, "base": sequence[p - 1], "region": region_of(p),
          "reason": reasons[p]} for p in sorted(reasons)])
    protected.to_csv(OUTPUT_DIR / "protected_positions.tsv",
                     sep="\t", index=False)

    free = [p for p in range(1, len(sequence) + 1) if p not in reasons]
    print("=" * 74)
    print("1. Positions that must not be modified")
    print("=" * 74)
    print(f"  protected : {len(reasons)} of {len(sequence)}")
    print(f"  available : {len(free)}")
    for region in ("spacer", "scaffold", "RTT", "PBS", "motif"):
        count = sum(1 for p in free if region_of(p) == region)
        print(f"    {region:<9} {count:3d} available")
    print()

    # ---------------------------------------------- 2. single-position scan
    print("=" * 74)
    print("2. Single-position scan over every unprotected position")
    print("=" * 74)

    rows = []
    for position in free:
        result = metrics(sequence, usable_pairs, [position])
        rows.append({
            "position": position,
            "base": sequence[position - 1],
            "region": region_of(position),
            **result,
            "delta_barrier": round(result["barrier"] - baseline["barrier"], 2),
        })

    scan = pd.DataFrame(rows).sort_values("delta_barrier")
    scan.to_csv(OUTPUT_DIR / "single_position_scan.tsv", sep="\t", index=False)
    print(scan.head(15).to_string(index=False))
    print()

    # ------------------------------------------- 3. combinations, safe zone
    print("=" * 74)
    print(f"3. Combinations within the distal PBS "
          f"({DISTAL_PBS[0]}-{DISTAL_PBS[1]})")
    print("=" * 74)
    print("This zone is chosen because shortening the PBS from its 3' end")
    print("is a published, functionally tolerated intervention, so weakening")
    print("the same region is the conservative version of the same idea.")
    print()

    candidates = list(range(DISTAL_PBS[0], DISTAL_PBS[1] + 1))
    rows = []
    for size in range(1, MAX_COMBINATION_SIZE + 1):
        for combo in combinations(candidates, size):
            result = metrics(sequence, usable_pairs, combo)
            rows.append({
                "n_modifications": size,
                "positions": ", ".join(str(p) for p in combo),
                "bases": "".join(sequence[p - 1] for p in combo),
                **result,
                "delta_barrier": round(
                    result["barrier"] - baseline["barrier"], 2),
            })

    combos = pd.DataFrame(rows).sort_values(
        ["delta_barrier", "n_modifications"])
    combos.to_csv(OUTPUT_DIR / "distal_pbs_combinations.tsv",
                  sep="\t", index=False)
    print(combos.head(10).to_string(index=False))
    print()

    # ------------------------------------------ 4. combining with the linkers
    print("=" * 74)
    print("4. Does it combine with the triazole linkers?")
    print("=" * 74)

    best = combos.iloc[0]
    best_positions = [int(p) for p in best["positions"].split(", ")]

    designs = [
        ("unmodified", []),
        ("linkers only", LINKER_JUNCTIONS),
        (f"modifications only ({best['positions']})", best_positions),
        ("modifications + linkers",
         sorted(set(best_positions + LINKER_JUNCTIONS))),
    ]

    rows = []
    for name, blocked in designs:
        rows.append({"design": name, **metrics(sequence, usable_pairs,
                                               blocked)})
    rows.append({"design": "plain sgRNA, for reference", **sgrna_reference})

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_DIR / "design_comparison.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print()

    print("Caveats that belong with these numbers:")
    print("  - A blocked position is an exact model of an abasic site and an")
    print("    upper bound for UNA, which destabilises by 4.0-6.6 kcal/mol")
    print("    internally rather than abolishing pairing entirely.")
    print("  - findpath is a heuristic, so every barrier is an upper bound.")
    print("  - None of this is validated experimentally. It ranks designs;")
    print("    it does not predict editing efficiency.")
    print()
    print(f"Files written to {OUTPUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
