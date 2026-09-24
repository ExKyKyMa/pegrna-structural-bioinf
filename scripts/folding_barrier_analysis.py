"""Folding Barrier and spacer accessibility analysis for PegRNA3.

Background
----------
Two different structures matter for a guide RNA.

The *ground state* is what the free RNA folds into on its own: the minimum
free energy (MFE) structure. The *active state* is the conformation the RNA
must adopt to be bound by Cas9, which is known experimentally from the
cryo-EM structure of the prime editor.

These are not the same structure. The question that predicts activity is not
"how stable is the RNA" but "how hard is it for the RNA to get from the
ground state into the active state". Zalatan and co-workers called that
quantity the Folding Barrier: the height of the highest energy point on the
refolding path. For CRISPRa they found Spearman r_s = 0.8 against activity,
explaining roughly 80% of the variation, with every effective guide below
about 10 kcal/mol and every defective guide above it. Thermodynamic
stability alone gave only r_s = 0.7, and general-purpose guide design
scores did far worse (Azimuth r_s = 0.22, Doench 2016 r_s = 0.02).

    Nature Communications 15, 6341 (2024)
    https://www.nature.com/articles/s41467-024-50528-1

That threshold was calibrated on ~100-nt bacterial CRISPRa scRNAs, not on a
158-nt pegRNA, so the absolute number should not be transferred without a
control. This script therefore always computes the plain sgRNA (spacer plus
scaffold, no 3' extension) alongside PegRNA3, so that the effect of the
extension can be separated from the effect of length and of the method.

What is computed
----------------
  1. the MFE ground state of the free RNA;
  2. the active state, built from the base pairs observed in the cryo-EM
     structure of the prime editor, under two definitions;
  3. the Folding Energy, the energy gap between them;
  4. the Folding Barrier, the saddle point on the refolding path,
     computed with ViennaRNA's findpath heuristic;
  5. spacer and seed accessibility from the partition function;
  6. the spacer:PBS intramolecular duplex, which is a documented
     auto-inhibition mechanism in prime editing, and how it responds to
     shortening the PBS.

Run from the repository root:

    python3 scripts/folding_barrier_analysis.py

Sources
-------
- Zalatan, J. G. et al. Guide RNA structure design enables combinatorial
  CRISPRa programs for biosynthetic profiling. Nat. Commun. 15, 6341 (2024).
  https://www.nature.com/articles/s41467-024-50528-1
- Wong, N., Liu, W. & Wang, X. WU-CRISPR: characteristics of functional
  guide RNAs for the CRISPR/Cas9 system. Genome Biol. 16, 218 (2015).
  https://link.springer.com/article/10.1186/s13059-015-0784-0
- Ponnienselvan, K. et al. Reducing the inherent auto-inhibitory interaction
  within the pegRNA enhances prime editing efficiency. Nucleic Acids Res.
  51, 6966-6980 (2023).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10359601/
- Shuto, Y. et al. Structural basis for pegRNA-guided reverse transcription
  by a prime editor. Nature 631, 224-231 (2024).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
- Lorenz, R. et al. ViennaRNA Package 2.0. Algorithms Mol. Biol. 6, 26
  (2011). https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26
"""

from pathlib import Path

import pandas as pd
import RNA

ROOT = Path(__file__).resolve().parent.parent

REFERENCE_FASTA = ROOT / "data/raw/sequences/PegRNA3_unmodified_reference.fasta"
BOUND_PAIRS_TSV = ROOT / "data/raw/designs/article_pegrna_bound_rna_pairs.tsv"
OUTPUT_DIR = ROOT / "results/folding_barrier"

# The cryo-EM pegRNA of Shuto et al. has a 20-nt spacer, PegRNA3 has 21, so
# every position in the published pair list shifts by one. The offset is
# verified at run time by locating the scaffold start in both sequences
# rather than being trusted as a constant.
SCAFFOLD_ANCHOR = "GUUUUAGAGCUAGA"

# Region boundaries in PegRNA3, from data/raw/designs/PegRNA3_regions.tsv.
SPACER_RANGE = (1, 21)
SCAFFOLD_RANGE = (22, 97)
RTT_RANGE = (98, 109)
PBS_RANGE = (110, 121)
MOTIF_RANGE = (122, 158)

# Seed region: the three nucleotides at the 3' end of the spacer. Wong et al.
# found their accessibility to be the single most discriminating structural
# feature between functional and non-functional guides.
SEED_RANGE = (19, 21)

# Junction positions where a triazole linker replaces the phosphodiester
# bond. The triazole cannot participate in base pairing, so the two flanking
# nucleotides are forced unpaired.
JUNCTION_SETS = {
    "two_linker": [(35, 36), (97, 98)],
    "three_linker": [(35, 36), (77, 78), (121, 122)],
    "four_linker": [(35, 36), (77, 78), (97, 98), (121, 122)],
    "linker_35_36_only": [(35, 36)],
    "linker_77_78_only": [(77, 78)],
    "linker_97_98_only": [(97, 98)],
    "linker_121_122_only": [(121, 122)],
}

# findpath search width. Higher explores more refolding paths and gives a
# tighter upper bound on the true barrier, at the cost of run time.
SEARCH_WIDTH = 400

CANONICAL_PAIRS = {"AU", "UA", "GC", "CG", "GU", "UG"}


# --------------------------------------------------------------- utilities

def read_fasta(path):
    lines = path.read_text().splitlines()
    return "".join(
        line.strip() for line in lines
        if line.strip() and not line.startswith(">")
    ).upper().replace("T", "U")


def model():
    md = RNA.md()
    md.temperature = 37.0
    return md


def pairs_to_structure(pairs, length):
    """Turn a list of (i, j) pairs into dot-bracket notation."""
    structure = ["."] * length
    for i, j in pairs:
        structure[i - 1] = "("
        structure[j - 1] = ")"
    return "".join(structure)


def region_of(position):
    for name, (start, end) in [
        ("spacer", SPACER_RANGE), ("scaffold", SCAFFOLD_RANGE),
        ("RTT", RTT_RANGE), ("PBS", PBS_RANGE), ("motif", MOTIF_RANGE),
    ]:
        if start <= position <= end:
            return name
    return "outside"


# ----------------------------------------------------------- active state

def load_active_pairs(sequence, article_sequence):
    """Build the active (Cas9-bound) pair list from the cryo-EM data.

    Two classes of observed pair cannot be represented by the
    nearest-neighbour energy model and are dropped, with the reason
    recorded so the exclusion is visible rather than silent:

      - non-canonical pairs, for which the model has no parameters;
      - pairs enclosing a hairpin loop shorter than 3 nt, which is
        sterically impossible in free RNA and exists in the structure only
        because the protein holds the backbone in a strained conformation.
    """
    offset = (sequence.index(SCAFFOLD_ANCHOR)
              - article_sequence.index(SCAFFOLD_ANCHOR))

    table = pd.read_csv(BOUND_PAIRS_TSV, sep="\t", dtype=str)
    kept, dropped = [], []

    for _, row in table.iterrows():
        i = int(row["position_i"]) + offset
        j = int(row["position_j"]) + offset

        if j > len(sequence):
            dropped.append((row["pair_id"], i, j, "--",
                            "outside this construct"))
            continue

        bases = sequence[i - 1] + sequence[j - 1]

        if bases not in CANONICAL_PAIRS:
            dropped.append((row["pair_id"], i, j, bases,
                            "non-canonical pair, no energy parameters"))
        elif j - i - 1 < 3:
            dropped.append((row["pair_id"], i, j, bases,
                            f"hairpin loop of {j - i - 1} nt, minimum is 3"))
        else:
            kept.append((i, j))

    return kept, dropped, offset


def build_relaxed_active(sequence, active_pairs, forced_unpaired=()):
    """Active state with the cryo-EM pairs enforced and the rest free.

    The strict definition leaves every nucleotide not seen pairing in the
    cryo-EM map unpaired. That is not what the bound molecule looks like:
    the map resolves the protein-contacted core, while the 3' extension is
    flexible and still forms structure. Enforcing the observed pairs,
    holding the spacer unpaired because Cas9 holds it extended for DNA
    interrogation, and letting everything else fold is the more physical
    model. Both are reported because the choice changes the number.
    """
    fc = RNA.fold_compound(sequence, model())

    for i, j in active_pairs:
        fc.hc_add_bp(i, j,
                     RNA.CONSTRAINT_CONTEXT_ALL_LOOPS
                     | RNA.CONSTRAINT_CONTEXT_ENFORCE)

    spacer_end = min(SPACER_RANGE[1], len(sequence))
    for position in range(SPACER_RANGE[0], spacer_end + 1):
        fc.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    for position in forced_unpaired:
        if position <= len(sequence):
            fc.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    structure, _ = fc.mfe()
    return structure


# ------------------------------------------------------------ core maths

def fold(sequence, forced_unpaired=()):
    """Fold with the given positions forced unpaired. Returns (struct, e)."""
    fc = RNA.fold_compound(sequence, model())

    for position in forced_unpaired:
        if position <= len(sequence):
            fc.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    return fc.mfe()


def saddle_energy(sequence, start_structure, end_structure):
    """Saddle-point energy on the direct refolding path, in kcal/mol.

    path_findpath_saddle returns energy in units of 0.01 kcal/mol. The
    result is an upper bound on the true barrier: findpath is a heuristic
    that explores a limited number of paths and may miss a lower one.
    """
    fc = RNA.fold_compound(sequence, model())
    saddle = fc.path_findpath_saddle(start_structure, end_structure,
                                     SEARCH_WIDTH)
    return saddle / 100.0


def unpaired_probabilities(sequence, forced_unpaired=()):
    """P(unpaired) for every position, from the partition function."""
    fc = RNA.fold_compound(sequence, model())

    for position in forced_unpaired:
        if position <= len(sequence):
            fc.hc_add_up(position, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    _, mfe_energy = fc.mfe()
    fc.exp_params_rescale(mfe_energy)
    fc.pf()
    bpp = fc.bpp()

    n = len(sequence)
    return {
        position: 1.0 - sum(
            bpp[min(position, other)][max(position, other)]
            for other in range(1, n + 1) if other != position
        )
        for position in range(1, n + 1)
    }


def mean_accessibility(probabilities, start, end):
    values = [probabilities[p] for p in range(start, end + 1)
              if p in probabilities]
    return sum(values) / len(values) if values else float("nan")


# ------------------------------------------------------------- variants

def unpaired_positions(junctions):
    """Positions a triazole linker forces to be unpaired.

    Nothing else is constrained. The ground state has to be the structure
    the free RNA actually adopts, so the spacer and the 3' motif are left
    free here even though both are constrained elsewhere in this repository
    for different questions.
    """
    positions = set()
    for left, right in junctions:
        positions.update((left, right))
    return sorted(positions)


def analyse(name, sequence, junctions, active_pairs, strict_active):
    """Compute every quantity for one variant, under both definitions."""
    forced = unpaired_positions(junctions)

    ground_structure, ground_energy = fold(sequence, forced)
    relaxed_active = build_relaxed_active(sequence, active_pairs, forced)

    free = RNA.fold_compound(sequence, model())
    strict_energy = free.eval_structure(strict_active)
    relaxed_energy = free.eval_structure(relaxed_active)

    strict_saddle = saddle_energy(sequence, ground_structure, strict_active)
    relaxed_saddle = saddle_energy(sequence, ground_structure, relaxed_active)

    accessibility = unpaired_probabilities(sequence, forced)

    return {
        "variant": name,
        "junctions": ", ".join(f"{a}/{b}" for a, b in junctions) or "none",
        "ground_energy": round(ground_energy, 2),
        "active_energy_strict": round(strict_energy, 2),
        "active_energy_relaxed": round(relaxed_energy, 2),
        "folding_energy_strict": round(strict_energy - ground_energy, 2),
        "folding_energy_relaxed": round(relaxed_energy - ground_energy, 2),
        "folding_barrier_strict": round(strict_saddle - ground_energy, 2),
        "folding_barrier_relaxed": round(relaxed_saddle - ground_energy, 2),
        "spacer_accessibility": round(
            mean_accessibility(accessibility, *SPACER_RANGE), 4),
        "seed_accessibility": round(
            mean_accessibility(accessibility, *SEED_RANGE), 4),
    }


# --------------------------------------------------------- spacer vs PBS

def spacer_pbs_analysis(sequence):
    """Quantify the intramolecular spacer:PBS duplex.

    The pegRNA PBS is by construction complementary to the nicked target
    DNA strand, which is itself complementary to the protospacer, which is
    the spacer sequence. The spacer and the PBS are therefore complementary
    to each other, and an RNA:RNA duplex is more stable than the RNA:DNA
    duplex the PBS is supposed to form with the DNA. Ponnienselvan et al.
    showed this auto-inhibition prevents R-loop formation and that
    shortening the PBS relieves it.
    """
    spacer = sequence[SPACER_RANGE[0] - 1:SPACER_RANGE[1]]
    rtt = sequence[RTT_RANGE[0] - 1:RTT_RANGE[1]]
    pbs = sequence[PBS_RANGE[0] - 1:PBS_RANGE[1]]

    duplex_pbs = RNA.duplexfold(spacer, pbs)
    duplex_rtt_pbs = RNA.duplexfold(spacer, rtt + pbs)

    rows = []
    full_length = PBS_RANGE[1] - PBS_RANGE[0] + 1

    for keep in range(full_length, 5, -1):
        # Trim the PBS from its 3' end, keeping RTT and the 3' motif.
        trimmed = (sequence[:PBS_RANGE[0] - 1 + keep]
                   + sequence[PBS_RANGE[1]:])
        probabilities = unpaired_probabilities(trimmed)
        duplex = RNA.duplexfold(spacer, pbs[:keep])
        _, energy = fold(trimmed)

        rows.append({
            "pbs_length": keep,
            "spacer_pbs_duplex_energy": round(duplex.energy, 2),
            "ground_energy": round(energy, 2),
            "spacer_accessibility": round(
                mean_accessibility(probabilities, *SPACER_RANGE), 4),
            "seed_accessibility": round(
                mean_accessibility(probabilities, *SEED_RANGE), 4),
        })

    return duplex_pbs, duplex_rtt_pbs, pd.DataFrame(rows)


def spacer_partners(sequence, threshold=0.05):
    """Which positions the spacer actually pairs with, and how strongly."""
    fc = RNA.fold_compound(sequence, model())
    _, mfe_energy = fc.mfe()
    fc.exp_params_rescale(mfe_energy)
    fc.pf()
    bpp = fc.bpp()

    n = len(sequence)
    rows = []
    mass = {}

    for position in range(SPACER_RANGE[0], SPACER_RANGE[1] + 1):
        for other in range(1, n + 1):
            if other == position:
                continue
            probability = bpp[min(position, other)][max(position, other)]
            region = region_of(other)
            mass[region] = mass.get(region, 0.0) + probability
            if probability > threshold:
                rows.append({
                    "spacer_position": position,
                    "spacer_base": sequence[position - 1],
                    "partner_position": other,
                    "partner_base": sequence[other - 1],
                    "partner_region": region,
                    "probability": round(probability, 4),
                })

    return pd.DataFrame(rows), mass


# ------------------------------------------------------------------ main

def main():
    sequence = read_fasta(REFERENCE_FASTA)
    n = len(sequence)

    article_sequence = read_fasta(
        ROOT / "data/raw/sequences/article_pegrna_reference.fasta")

    print(f"PegRNA3, {n} nt")
    print()

    active_pairs, dropped, offset = load_active_pairs(sequence,
                                                      article_sequence)

    print("Active state from the cryo-EM observed base pairs")
    print(f"  numbering offset, article to repository : +{offset}")
    print(f"  pairs in the published list             : "
          f"{len(active_pairs) + len(dropped)}")
    print(f"  usable in the nearest-neighbour model   : {len(active_pairs)}")

    if dropped:
        print("  excluded:")
        for pair_id, i, j, bases, reason in dropped:
            print(f"    {pair_id}  {i}-{j}  {bases}  {reason}")

    strict_active = pairs_to_structure(active_pairs, n)
    print()
    print("  " + strict_active)
    print()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------- 1. the sgRNA control
    print("=" * 78)
    print("1. Control: plain sgRNA, spacer plus scaffold, no 3' extension")
    print("=" * 78)

    sgrna = sequence[:SCAFFOLD_RANGE[1]]
    sgrna_pairs, _, _ = load_active_pairs(sgrna, article_sequence)
    sgrna_strict = pairs_to_structure(sgrna_pairs, len(sgrna))
    sgrna_row = analyse("sgRNA_control", sgrna, [], sgrna_pairs, sgrna_strict)

    print(f"  length                 : {len(sgrna)} nt")
    print(f"  ground state energy    : {sgrna_row['ground_energy']}")
    print(f"  Folding Barrier, relaxed: "
          f"{sgrna_row['folding_barrier_relaxed']}")
    print(f"  spacer accessibility   : {sgrna_row['spacer_accessibility']}")
    print(f"  seed accessibility     : {sgrna_row['seed_accessibility']}")
    print()

    # -------------------------------------------- 2. PegRNA3 and variants
    print("=" * 78)
    print("2. PegRNA3 and its linker variants")
    print("=" * 78)

    rows = [sgrna_row,
            analyse("unmodified", sequence, [], active_pairs, strict_active)]

    for name, junctions in JUNCTION_SETS.items():
        rows.append(analyse(name, sequence, junctions, active_pairs,
                            strict_active))

    results = pd.DataFrame(rows)
    barrier_file = OUTPUT_DIR / "PegRNA3_folding_barrier.tsv"
    results.to_csv(barrier_file, sep="\t", index=False)

    print(results.to_string(index=False))
    print()

    reference = results.loc[results.variant == "unmodified"].iloc[0]

    print("Change in Folding Barrier relative to unmodified PegRNA3:")
    print(f"  {'variant':<22} {'strict':>8} {'relaxed':>9}")
    for _, row in results.iterrows():
        if row.variant in ("unmodified", "sgRNA_control"):
            continue
        print(f"  {row.variant:<22} "
              f"{row.folding_barrier_strict - reference.folding_barrier_strict:+8.2f} "
              f"{row.folding_barrier_relaxed - reference.folding_barrier_relaxed:+9.2f}")
    print()

    extension_cost = (reference.folding_barrier_relaxed
                      - sgrna_row["folding_barrier_relaxed"])
    print(f"Cost of the 3' extension: the barrier rises by "
          f"{extension_cost:+.2f} kcal/mol going from the plain sgRNA "
          f"to PegRNA3.")
    print()

    # --------------------------------------------- 3. what blocks the spacer
    print("=" * 78)
    print("3. What sequesters the spacer")
    print("=" * 78)

    partners, mass = spacer_partners(sequence)
    partners_file = OUTPUT_DIR / "PegRNA3_spacer_pairing_partners.tsv"
    partners.to_csv(partners_file, sep="\t", index=False)

    print("Total spacer pairing probability, summed by partner region:")
    for region, value in sorted(mass.items(), key=lambda item: -item[1]):
        print(f"  {region:<10} {value:7.3f}")
    print()

    print("Strong spacer contacts, P > 0.05:")
    print(partners.to_string(index=False))
    print()

    # ------------------------------------------------ 4. the spacer:PBS duplex
    print("=" * 78)
    print("4. The spacer:PBS auto-inhibitory duplex")
    print("=" * 78)

    duplex_pbs, duplex_rtt_pbs, scan = spacer_pbs_analysis(sequence)
    scan_file = OUTPUT_DIR / "PegRNA3_pbs_length_scan.tsv"
    scan.to_csv(scan_file, sep="\t", index=False)

    print(f"  spacer : PBS          dG = {duplex_pbs.energy:6.2f} kcal/mol")
    print(f"  spacer : RTT + PBS    dG = {duplex_rtt_pbs.energy:6.2f} "
          f"kcal/mol")
    print()
    print("Effect of shortening the PBS from its 3' end. This is a design")
    print("calculation only; the current construct is already synthesised.")
    print(scan.to_string(index=False))
    print()

    print("Files written:")
    for path in (barrier_file, partners_file, scan_file):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
