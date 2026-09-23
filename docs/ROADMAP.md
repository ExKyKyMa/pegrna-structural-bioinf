# Roadmap: from sequence to molecular dynamics to thesis defence

Status date: 2026-09-23. Assumed defence: spring 2027.

## The claim the thesis has to defend

A master's thesis is judged on one defensible claim, not on the number of
calculations performed. The claim this project is positioned to make is:

> A pegRNA can be assembled from separately synthesised fragments joined
> by non-natural triazole linkers. Computational analysis identifies
> which junction positions are structurally permissible, predicts which
> one additionally suppresses reverse-transcriptase over-extension into
> the scaffold, and specifies the linker length required to span each
> junction without distorting the RNA.

Everything below exists to support or falsify that sentence. Any result
that does neither is optional.

Three falsifiable predictions should be stated early and tested:

| # | Prediction | Test |
|---|---|---|
| P1 | Junctions 35/36 and 121/122 do not perturb the fold; 97/98 does | base-pair persistence in MD, versus unmodified control |
| P2 | A triazole linker at 97/98 blocks scaffold-derived incorporations | RT extension assay, or RT footprint in the modelled complex |
| P3 | The linker must span the measured P-to-P distance at each junction; a linker that is too short strains the backbone | RDKit conformer end-to-end distribution versus cryo-EM distances, then MD |

## Stage 0 — close the data layer (1-2 weeks)

Small, cheap, and blocking everything downstream.

- [ ] Rename `Evopreq1` to `tevopreQ1` across data files and docs, with a
      note recording the old name.
- [x] **Split the extension (98-121) into RTT and PBS.** Done 2026-09-23:
      RTT = 98-109 (12 nt), PBS = 110-121 (12 nt). The 11+13 split
      recorded earlier was off by one; position 109 is the first
      templating nucleotide, not part of the primer duplex.
- [x] Record the target genomic site and the intended edit in
      `data/raw/designs/`. Done: `PegRNA3_regions.tsv` and
      `PegRNA3_target_edit.tsv`. SNCA protospacer c.69-c.89, PAM AGG at
      c.90-c.92, nick c.86/c.87, edit c.88G>C = p.A30P (installed, not
      corrected — still to be confirmed with the requesting group).
- [ ] Add a GitHub Actions workflow running `scripts/check_pegrna_design.py`
      on every push, so a sequence error like the one found this week
      cannot survive a commit again.
- [ ] Extend the validator to check that fragments concatenate back to
      the reference exactly, not only that lengths sum correctly. The
      two compensating errors passed a length check.

## Stage 1 — finish the secondary-structure argument (2-3 weeks)

- [ ] Run **pegLIT** on this construct to obtain the nucleotide linker
      the published method would recommend at 121/122. This gives a
      reference point: the chemical linker can then be presented as a
      non-pairing limiting case of an established design principle,
      which is far stronger than presenting it as an invention.
- [ ] Resolve the pseudoknot register with a pseudoknot-capable
      predictor (IPknot, pKiss, HotKnots) instead of the three-register
      sensitivity test now in the repository.
- [ ] Positive control: repeat the whole junction analysis on the plain
      sgRNA at positions 34/35 and 57/58, where fragment ligation is
      already known experimentally not to reduce activity. If the method
      passes the control it is credible for the untested junctions.
- [ ] Write the limitations section now, while the artifacts are fresh.
      The honest treatment of the spacer artifact and the pseudoknot
      artifact is one of the strongest parts of this project and should
      be presented as method, not as apology.

## Stage 2 — chemistry with RDKit (3-4 weeks)

RDKit is the right tool here, but for a narrower job than it may appear.
RDKit is a cheminformatics toolkit: it handles molecular graphs,
reactions, descriptors and conformer generation. It does **not** do
quantum chemistry and it does **not** do molecular dynamics. Its role in
this project is to turn the prose descriptions in
`data/raw/designs/PegRNA3_linkers.tsv` into exact, machine-checkable
molecules, and then to answer one geometric question.

### 2.1 Validate the component structures

`data/raw/components.tsv` is still a template with empty SMILES columns,
and `docs/components/*.md` files carry the note "SMILES status:
unverified. The structure should be checked computationally with RDKit."

A script `scripts/validate_components.py` should, for every component:

- parse and sanitise the SMILES, failing loudly if it does not parse;
- compute the molecular formula and exact mass and compare them with the
  values written in the component's markdown file;
- emit the canonical SMILES and the InChIKey, so the identity is
  citable and diffable;
- write everything back into `data/raw/components.tsv` with a
  `smiles_status` of `verified` or `mismatch`.

This converts ten documentation files with unverified formulas into a
validated small-molecule inventory. It is unglamorous and it is exactly
the kind of thing a committee checks.

### 2.2 Build the linker *products*, not the reagents

The two linkers are currently described only as text:

| Linker | Junction | Description |
|---|---|---|
| Linker_1 | scaffold F1 to scaffold F2 | 3'-azide - propargylamine - C(O)- |
| Linker_2 | scaffold F2 to extension | 3'-amino - azidobutanoic acid - propargylamine - C(O)- |

Encode the copper-catalysed azide-alkyne cycloaddition as a reaction
SMARTS and let RDKit generate the 1,4-triazole product from the
reagents. The output is an explicit atom-by-atom structure of the
junction, which is what any later modelling step needs. Write it to
`data/processed/linkers/Linker_1_product.smi` and
`Linker_2_product.smi`, with a unit test asserting the expected formula.

Also build a **capped junction fragment**: the triazole linker plus one
nucleotide on each side, with the RNA cut points capped by methyl
groups. This capped fragment, not the bare linker, is the molecule that
later gets quantum-chemical treatment and force-field parameters.

### 2.3 Answer the geometric question

This is the step that connects chemistry to structural biology and it is
the heart of the RDKit work.

1. Generate conformers of each linker product with ETKDGv3, optimise
   with MMFF94s, and cluster.
2. For each conformer, measure the distance between the two atoms that
   attach to the flanking nucleotides. This gives the **distribution of
   spans** the linker can adopt, and the number of rotatable bonds tells
   you how much entropy is paid to reach the extremes.
3. Download PDB **8WUS** (termination state) and measure the actual
   phosphorus-to-phosphorus distance across each candidate junction:
   35/36, 77/78, 97/98, 121/122. A normal RNA P-to-P step is about 6 Å;
   a junction in a loop can be longer.
4. Compare. A linker whose accessible span comfortably contains the
   measured distance is geometrically permissible. One that must sit at
   the tail of its distribution will strain the backbone, and that
   strain is what MD will show.

The deliverable is one table: junction, measured distance, linker span
range, verdict. That table is a thesis figure.

### 2.4 Refine geometry and charges

RDKit conformers from MMFF are adequate for measuring spans but not for
force-field parameterisation. Take the best conformer of the capped
junction fragment into **xTB** (GFN2-xTB) for geometry optimisation, and
compute the torsion profile of the two or three dihedrals around the
triazole. That profile is the reference the force field must reproduce.

## Stage 3 — everything that must happen before MD (4-6 weeks)

This stage is usually underestimated. The list below is the actual
blocking work.

### 3.1 Decide what system is being simulated

Two options, and they answer different questions.

| System | Answers | Cost |
|---|---|---|
| Free modified pegRNA in solution | does the linker perturb the fold, does the motif stay autonomous | ~60-80k atoms, cheap |
| pegRNA inside the SpCas9 - M-MLV RT complex, from 8WUS | does the linker survive protein contacts, does 97/98 block RT | ~300-400k atoms, expensive |

Recommendation: do the free RNA first for all four junctions, then the
complex for only the best and the worst junction plus the unmodified
control. Trying to do the complex for everything will consume the whole
year.

### 3.2 Obtain a 3D starting structure

- For the complex: use the deposited coordinates (8WUS, 8WUT, 8WUU,
  8WUV, 8YGJ). Note that the deposited pegRNA is not identical to
  PegRNA3, so the spacer, RTT and PBS must be mutated in silico and the
  3' motif, which is disordered or absent in cryo-EM, must be modelled.
- For the free RNA: predict 3D from the secondary structure already
  computed, using RNAComposer, FARFAR2, trRosettaRNA, AlphaFold3 or
  Boltz. Generate several models and keep the ensemble; do not pretend a
  single predicted RNA 3D model is correct.
- Restrain the prediction with the secondary structure derived in
  Stage 1, including the pseudoknot, rather than folding blind.

### 3.3 Splice the linker into the 3D model

Delete the phosphodiester group at the junction, place the linker
conformer whose span matches the measured distance, join the atoms,
and locally minimise while restraining the rest of the molecule. Do this
with a script, not by hand, so it is reproducible and so all four
junctions are built the same way.

### 3.4 Force-field parameters — the real blocker

No standard RNA force field contains triazole-linker parameters. This is
the single most likely place for the project to stall, so plan for it.

- RNA: AMBER OL3 (chi-OL3) or CHARMM36. OL3 is the usual choice for RNA.
- Linker: GAFF2 with AM1-BCC charges, or OpenFF Sage. Generate with
  `antechamber` / `acpype`, or the OpenFF toolkit.
- The boundary is the difficult part. The linker becomes a custom
  residue, and the bonded terms crossing from the linker into the
  adjacent nucleotides must be defined by hand. Document every one of
  them; a committee will ask.
- Validate the parameters before using them: scan the key dihedrals with
  the force field and compare against the xTB or DFT profile from Stage
  2.4. If they disagree, refit. An unvalidated custom parameter set is
  the easiest thing for a reviewer to reject.

### 3.5 System preparation

- Water: OPC or TIP3P. OPC behaves better for RNA.
- Ions: neutralising K+ plus excess KCl, and explicit **Mg2+**. RNA
  folding is Mg-dependent and standard Mg parameters are known to be
  poor; use a corrected set and say which one and why.
- Box: at least 12 Å padding, ideally more for an extended pegRNA.
- Protocol: minimise, heat under restraints, equilibrate in NPT with
  restraints, release gradually, then production in NPT.

### 3.6 Sampling plan, decided in advance

- At least **three independent replicas** per system with different
  velocity seeds. One trajectory is an anecdote.
- 500 ns to 1 microsecond per replica for the free RNA.
- The unmodified pegRNA must be simulated under identical conditions.
  Without that control, nothing about the linkers can be concluded.
- Hardware: this will not run on a laptop CPU. OpenMM on a single
  consumer GPU is the practical minimum; budget weeks of wall-clock
  time, and secure access before Stage 3 ends rather than after.

### 3.7 Define the analysis before running anything

Decide the observables first, so the simulations are designed to answer
the question rather than being mined afterwards.

- Per-region RMSD and RMSF, with the regions already defined in this
  repository (S, F1, F2, E, V).
- Base-pair occupancy over time, especially the SpCas9-recognised stem
  83-86 : 94-97 and the tevopreQ1 pseudoknot. This is the direct MD
  counterpart of the P(paired) numbers already computed.
- Local backbone geometry at the junction: P-to-P distance, torsion
  distributions, and whether they stay inside the ranges sampled by the
  isolated linker.
- Motif autonomy: contacts between the 3' motif and the rest of the
  molecule, which is the MD test of the hypothesis argued in
  `PegRNA3_pseudoknot_constrained_analysis.md`.
- For the complex: persistence of the SpCas9 contacts reported by Shuto
  et al. (Y64 and L99 to the 3' stem loop, R116 hydrogen bond to C96).

## Stage 4 — experiment, if at all possible (parallel)

One experimental result changes the character of the defence. In order
of increasing effort:

1. Denaturing gel showing that the click ligation produced a product of
   the expected length. Cheap and immediately convincing.
2. In vitro RT extension assay on the two-linker construct, reading out
   whether scaffold-derived incorporation is reduced. This directly tests
   prediction P2 and is the most valuable single experiment available.
3. Structure probing (SHAPE or DMS-MaPseq) of the modified and
   unmodified pegRNA, which validates the secondary-structure
   predictions and resolves the pseudoknot register empirically.
4. Cell-based editing efficiency, if a collaborator can run it.

A purely computational thesis is defensible, but a computational thesis
with one gel is much harder to argue with.

## Stage 5 — writing (last 2 months, started much earlier)

- Keep writing the weekly logs. They are already most of a methods
  section, and they carry dates, which protects the work's priority.
- Structure: definition of the construct, secondary-structure analysis
  with explicit model limitations, chemical and geometric validation of
  the linkers, MD, synthesis of the three predictions, limitations,
  outlook.
- The reproducible repository is itself a defensible contribution. Give
  it a release tag and a DOI through Zenodo and cite it in the thesis.
- Prepare answers to the three questions a committee will certainly ask:
  why should a secondary-structure prediction be trusted for a molecule
  that is protein-bound and pseudoknotted; where did the force-field
  parameters for the triazole come from and how were they validated; and
  what would falsify the conclusion.

## Sequencing summary

| Period | Focus | Output |
|---|---|---|
| Oct 2026 | Stage 0 and Stage 1 | defined edit, CI validation, pegLIT comparison, sgRNA control |
| Nov-Dec 2026 | Stage 2 | verified component inventory, linker products, span-versus-distance table, xTB reference profiles |
| Jan-Feb 2027 | Stage 3 | built and parameterised systems, validated parameters, equilibrated boxes |
| Feb-Apr 2027 | production MD and analysis | three replicas per system, occupancy and geometry analysis |
| Apr-May 2027 | writing | thesis, repository release, defence preparation |
| throughout | Stage 4 | any experimental data obtainable |

## Sources

- Shuto, Y. et al. Structural basis for pegRNA-guided reverse
  transcription by a prime editor. *Nature* **631**, 224-231 (2024).
  https://doi.org/10.1038/s41586-024-07497-8
  PDB 8WUS, 8WUT, 8WUU, 8WUV, 8YGJ; EMDB EMD-37858 to EMD-37861 and
  EMD-39253. Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC11222144/
- Nelson, J. W. et al. Engineered pegRNAs improve prime editing
  efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8930418/
- pegLIT. https://github.com/sshen8/peglit and https://peglit.liugroup.us/
- Structure and function of preQ1 riboswitches.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4177978/
- RDKit documentation. https://www.rdkit.org/docs/
- xTB documentation. https://xtb-docs.readthedocs.io/
- OpenMM documentation. https://openmm.org/documentation
