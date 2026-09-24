"""Validate every PegRNA3 sequence file against each other and against the
real SNCA reference.

Run from the repository root:

    python3 scripts/validate_pegrna3.py

Exit code 0 means every check passed. Exit code 1 means at least one failed,
which is what lets this script be used as an automatic gate before a commit.

There are two kinds of check here and the difference matters.

Internal consistency checks ask "do the files in this repository agree with
each other". They catch copy-paste damage, but they cannot catch an error
that was made consistently everywhere.

External checks ask "does the construct agree with the real SNCA sequence".
These are the ones that can catch a genuinely wrong design. They use a
stored window of NM_000345.4 so that the script runs offline and so that the
reference itself is version controlled.
"""

from pathlib import Path
import sys

# Paths are resolved relative to the repository root, which is the parent of
# the scripts/ directory. This makes the script work no matter which
# directory it is called from.
ROOT = Path(__file__).resolve().parent.parent

REFERENCE_FASTA = ROOT / "data/raw/sequences/PegRNA3_unmodified_reference.fasta"
MODIFIED_FASTA = ROOT / "data/raw/sequences/PegRNA3_modified_nucleotide_reference.fasta"
REGIONS_TSV = ROOT / "data/raw/designs/PegRNA3_regions.tsv"
FRAGMENTS_TSV = ROOT / "data/raw/designs/PegRNA3_fragments.tsv"
TARGET_TSV = ROOT / "data/raw/designs/PegRNA3_target_edit.tsv"
MFE_INPUT = ROOT / "results/secondary_structure/PegRNA3_MFE_input.txt"
SNCA_FASTA = ROOT / "data/reference/SNCA_NM_000345.4_CDS_c40_c140.fasta"

# The stored SNCA window starts at this CDS coordinate.
SNCA_WINDOW_FIRST = 40

# The canonical SpCas9 sgRNA scaffold, written 5' to 3'. Hard-coded on
# purpose: this is an external constant, not something derived from our own
# files, so it can catch an error that was propagated through all of them.
CANONICAL_SCAFFOLD = (
    "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"
)

EXPECTED_LENGTH = 158
ALLOWED = set("AUGC")

COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}

checks = []


def check(name, ok, detail=""):
    """Record one check. Printing happens at the end so the report is tidy."""
    checks.append((name, bool(ok), detail))


def to_rna(text):
    return "".join(text.split()).upper().replace("T", "U")


def read_fasta(path):
    lines = path.read_text().splitlines()
    return to_rna("".join(l for l in lines if l and not l.startswith(">")))


def read_tsv(path):
    """Minimal TSV reader: returns a list of dicts keyed by column name."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def reverse_complement(rna):
    return "".join(COMPLEMENT[b] for b in reversed(rna))


def snca(start, end):
    """Return SNCA CDS positions start..end (c. numbering, inclusive) as RNA."""
    window = read_fasta(SNCA_FASTA)
    lo = start - SNCA_WINDOW_FIRST
    hi = end - SNCA_WINDOW_FIRST + 1
    if lo < 0 or hi > len(window):
        raise ValueError(
            f"c.{start}-c.{end} is outside the stored window "
            f"c.{SNCA_WINDOW_FIRST}-c.{SNCA_WINDOW_FIRST + len(window) - 1}"
        )
    return window[lo:hi]


def main():
    full = read_fasta(REFERENCE_FASTA)

    # ---------------------------------------------------------------- basics
    check("full length is 158 nt", len(full) == EXPECTED_LENGTH,
          f"got {len(full)}")
    check("alphabet is AUGC only", set(full) <= ALLOWED,
          f"unexpected symbols: {sorted(set(full) - ALLOWED)}")

    # ------------------------------------------- canonical scaffold, 22 to 97
    # An external constant, so this check is independent of our own tables.
    check("canonical SpCas9 scaffold at 22-97",
          full[21:97] == CANONICAL_SCAFFOLD,
          "scaffold does not match the canonical sequence")

    # --------------------------------------------------- regions versus FASTA
    regions = {r["region"]: r for r in read_tsv(REGIONS_TSV)}

    for name, row in regions.items():
        start, end = int(row["start_nt"]), int(row["end_nt"])
        stated_length = int(row["length_nt"])
        stored = to_rna(row["sequence"])
        check(f"region {name}: sequence matches reference at {start}-{end}",
              full[start - 1:end] == stored,
              f"reference has {full[start - 1:end]}, table has {stored}")
        check(f"region {name}: length is self-consistent",
              len(stored) == stated_length == end - start + 1,
              f"sequence {len(stored)} nt, stated {stated_length}, "
              f"coordinates span {end - start + 1}")

    # Regions must tile the whole molecule with no gap and no overlap.
    spans = sorted((int(r["start_nt"]), int(r["end_nt"]))
                   for r in regions.values())
    position = 1
    tiling_ok = True
    for start, end in spans:
        if start != position:
            tiling_ok = False
            break
        position = end + 1
    check("regions tile positions 1-158 without gaps or overlaps",
          tiling_ok and position - 1 == len(full),
          f"coverage ends at {position - 1}")

    # ------------------------------------------------- fragments versus FASTA
    # PegRNA3_fragments.tsv has no coordinate columns, so positions are
    # derived from the order column and then checked against the reference.
    fragments = sorted(read_tsv(FRAGMENTS_TSV), key=lambda r: int(r["order"]))
    cursor = 1
    for row in fragments:
        stored = to_rna(row["sequence"])
        end = cursor + len(stored) - 1
        check(f"fragment {row['fragment_name']}: matches reference "
              f"at {cursor}-{end}",
              full[cursor - 1:end] == stored,
              "fragment does not match the reference at its implied position")
        cursor = end + 1
    check("fragments concatenate to the full sequence",
          cursor - 1 == len(full),
          f"fragments total {cursor - 1} nt, reference is {len(full)} nt")

    # ------------------------------------------------------- derived files
    check("modified-nucleotide FASTA matches the reference",
          read_fasta(MODIFIED_FASTA) == full,
          "modified reference differs from the unmodified reference")

    mfe_lines = [l for l in MFE_INPUT.read_text().splitlines() if l.strip()]
    check("RNAfold input file matches the reference",
          any(to_rna(l) == full for l in mfe_lines),
          "no line in the RNAfold input equals the reference sequence")

    # =====================================================================
    # External checks against the real SNCA sequence.
    #
    # These are the checks that can detect a wrong design rather than an
    # inconsistent one.
    # =====================================================================
    if not SNCA_FASTA.exists():
        check("SNCA reference is present", False, f"{SNCA_FASTA} not found")
    else:
        spacer = to_rna(regions["Spacer"]["sequence"])
        rtt = to_rna(regions["RTT"]["sequence"])
        pbs = to_rna(regions["PBS"]["sequence"])

        # The spacer equals the protospacer on the coding strand.
        check("spacer equals SNCA protospacer c.69-c.89",
              spacer == snca(69, 89),
              f"SNCA has {snca(69, 89)}")

        # SpCas9 requires NGG immediately 3' of the protospacer.
        pam = snca(90, 92)
        check("PAM at c.90-c.92 is NGG", pam[1:] == "GG", f"PAM is {pam}")

        # The PBS anneals to the nicked strand, so it must be the reverse
        # complement of the genomic sequence immediately 5' of the nick.
        #
        # Note this is a check against the genome, not against the spacer.
        # Checking the PBS against the spacer would be wrong: the PBS pairs
        # with target DNA, and it only resembles the spacer region because
        # the two overlap on the same locus.
        check("PBS is the reverse complement of SNCA c.75-c.86",
              reverse_complement(pbs) == snca(75, 86),
              f"PBS reverse complement is {reverse_complement(pbs)}, "
              f"SNCA c.75-c.86 is {snca(75, 86)}")

        # Reverse transcription copies the RTT, so the newly written flap is
        # the reverse complement of the RTT. Compare it with the original
        # genomic sequence over the same span to recover the intended edit.
        flap = reverse_complement(rtt)
        original = snca(87, 87 + len(flap) - 1)
        differences = [
            (87 + i, o, n)
            for i, (o, n) in enumerate(zip(original, flap))
            if o != n
        ]
        check("reverse transcription introduces exactly one change",
              len(differences) == 1,
              f"found {len(differences)}: " + ", ".join(
                  f"c.{p} {o}>{n}" for p, o, n in differences))

        if len(differences) == 1:
            position, old, new = differences[0]
            check("the single change is c.88G>C",
                  (position, old, new) == (88, "G", "C"),
                  f"found c.{position}{old}>{new}")

            # Codon 30 spans c.88-c.90. Confirm the protein consequence
            # rather than trusting the note written in the table.
            codon_before = snca(88, 90).replace("U", "T")
            codon_after = (new + snca(89, 90)).replace("U", "T")
            check("codon 30 changes GCA (Ala) to CCA (Pro)",
                  codon_before == "GCA" and codon_after == "CCA",
                  f"{codon_before} -> {codon_after}")

        # Cross-check the free-text claims recorded in the target table, so
        # that the table cannot silently drift away from the sequences.
        if TARGET_TSV.exists():
            target = {r["field"]: r["value"] for r in read_tsv(TARGET_TSV)}
            check("target table records edit c.88G>C",
                  target.get("edit", "").strip() == "c.88G>C",
                  f"table says {target.get('edit')!r}")
            check("target table records p.Ala30Pro",
                  "Ala30Pro" in target.get("protein_consequence", ""),
                  f"table says {target.get('protein_consequence')!r}")

    # ------------------------------------------------------------- reporting
    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        line = ("PASS  " if ok else "FAIL  ") + name.ljust(width)
        if not ok and detail:
            line += f"   <- {detail}"
        print(line.rstrip())

    passed = sum(1 for _, ok, _ in checks if ok)
    print(f"\n{passed}/{len(checks)} checks passed")
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
