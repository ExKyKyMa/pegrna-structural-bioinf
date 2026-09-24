"""Validate PegRNA3 sequence files. Run:  python validate_pegrna3.py"""
from pathlib import Path
import sys

HERE = Path(__file__).parent
FULL = (HERE / "PegRNA3_full_sequence.txt").read_text().strip()
ALLOWED = set("AUGC")
errors, checks = [], []

def check(name, ok, detail=""):
    checks.append((name, ok, detail))
    if not ok:
        errors.append(f"{name}: {detail}")

check("full length == 158", len(FULL) == 158, f"got {len(FULL)}")
check("alphabet AUGC only", set(FULL) <= ALLOWED, f"bad: {sorted(set(FULL)-ALLOWED)}")

SCAFFOLD = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"
check("canonical SpCas9 scaffold at 22-97", FULL[21:97] == SCAFFOLD, "scaffold mismatch")

def read_tsv(p):
    lines = (HERE / p).read_text().strip().split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]

for r in read_tsv("PegRNA3_regions.tsv"):
    a, b = int(r["start_nt"]), int(r["end_nt"])
    check(f"region {r['region']} coords", FULL[a-1:b] == r["sequence"], "slice != sequence")
    check(f"region {r['region']} length", len(r["sequence"]) == int(r["length_nt"]) == b-a+1, "length mismatch")

frags = {}
for r in read_tsv("PegRNA3_fragments.tsv"):
    a, b = int(r["start_nt"]), int(r["end_nt"])
    check(f"{r['construct']}/{r['fragment']} coords", FULL[a-1:b] == r["sequence"], "slice != sequence")
    check(f"{r['construct']}/{r['fragment']} length", len(r["sequence"]) == int(r["length_nt"]), "length mismatch")
    frags.setdefault(r["construct"], []).append((int(r["order"]), r["sequence"]))
for c, fs in frags.items():
    joined = "".join(s for _, s in sorted(fs))
    check(f"{c} reassembles to full length", joined == FULL, f"got {len(joined)} nt")

fa = (HERE / "PegRNA3_unmodified_reference.fasta").read_text().strip().split("\n")
check("FASTA == full sequence", "".join(fa[1:]) == FULL, "fasta mismatch")
mfe = (HERE / "PegRNA3_MFE_input.txt").read_text().strip().split("\n")
check("MFE input == full sequence", mfe[1] == FULL, "mfe input mismatch")

comp = {"A":"U","U":"A","G":"C","C":"G"}
pbs = next(r for r in read_tsv("PegRNA3_regions.tsv") if r["region"] == "PBS")["sequence"]
rc = "".join(comp[c] for c in reversed(pbs))
check("PBS complementary to spacer 3' end", rc in FULL[:21], f"revcomp {rc} not in spacer")

for name, ok, detail in checks:
    print(("PASS  " if ok else "FAIL  ") + name + ("" if ok else f"  <- {detail}"))
print(f"\n{sum(1 for _,ok,_ in checks if ok)}/{len(checks)} checks passed")
sys.exit(1 if errors else 0)
