# Repeats specification

**Project:** 1, Repeats
**Repository:** `ga-repeats` (https://github.com/ptmcgrat/ga-repeats)
**Document version:** 0.11, 2026-11-06 (0.1: first draft from coordinator specification v0.4; 0.2: delivery text from v0.5; 0.3: relatedness context from v0.7; 0.4: the P1 genomes are downloaded, from v0.8; 0.5: the first conversation starts at the briefing; working with a student; 0.6: uploading wording corrected; 0.7: the student uploads the PDFs; 0.8: using the common library; 0.9: first conversation done; the owner's decisions on the question 1 pilot, the tools and Repbase; prior-art versions corrected; 0.10: masked genomes are not stored (coordinator answer); Repbase licence and location are not a concern; runtime and `ga_common` in the `ga-repeats` environment; first pilot run started; 0.11: the three pilot libraries built; comparisons, masking and manifests scripted (`repeats_009`); Dfam 4.0 classification reference; the student pushes from the laptop)
**Owner:** McGrath lab
**Status:** question 1 pilot: libraries built (`RI_v1`–`RI_v3`); comparisons and masking next, with `repeats_009`; **section 8 says where things stand**

This is the single specification document of project 1. It is always replaced whole.
It records what the coordinator specification already fixes, and leaves everything
else open. The methods are not designed here: the project designs them with the owner.

---

## 0. Where the first conversation starts, and who is working

**The person working on this project is a student, supervised by the owner (the professor).** The assistant
acts as a mentor and checker: it explains, asks questions, and checks the student's understanding against the
documentation and against evidence the student pastes. The student proposes; **the owner decides.** Where a choice
is the owner's, the assistant writes it into section 7 as an open question and does not settle it.

**Read at the start of every session:** this specification and `docs/repeats_briefing.md`.

**Before the first conversation, the student uploads the reference PDFs** listed in the briefing, Part D (items 1–5 are
required, 6–8 may follow). The assistant's first action is to list the documents it can see in the project, with the title
of each, and to ask for anything missing or unreadable. It does not go on until items 1–5 are there.

**The first conversation follows this order:**
1. The student reads Quah et al. 2025 (*Genome Research* 35:1094), in particular the Results section "Enrichment
   of TEs in SVs" and the Methods section "Transposable elements". The assistant goes through Part A of the
   briefing with the student, including the discussion points in A5.
2. **Question 1:** how many genomes should the library be built from (briefing, Part B)?
3. **Question 2:** how confident can we be that the library finds all transposon families, and how would we measure
   it (briefing, Part C)?
4. The student writes their proposals for questions 1 and 2 into sections 6 and 7 for the owner to review.

**No installation and no RepeatModeler run** until the gate 0 questions (section 4a) are answered and the owner has
agreed the plan for questions 1 and 2, unless the owner says otherwise.

*2026-10-28:* gate 0 is answered; the owner has agreed the question 1 pilot and chosen to evaluate question 2 after the
library is built (section 6). Installation and the pilot runs proceed gate by gate (section 4a).

Anything that needs another project goes to the coordinator (blocking items and integration requests only).

---

## 1. Purpose and scope

**Does:**
- build **one shared repeat library** for all assemblies of the panel;
- **mask every assembly** with it;
- produce, per assembly, the repeat annotation (contract `repeats`: RepeatMasker `.out`
  plus BED; interval, class, family, divergence).

**Does not:** align, annotate genes, or call structure; it designs no other project's work.
Projects A, 2a, 2b and 3 consume its output through the coordinator.

---

## 2. Inputs

| input | version | location | state |
|---|---|---|---|
| assemblies: one genome FASTA per panel assembly | panel P1 (25 assemblies, by accession version) | `<data root>/Genomes/<assembly_id>/<accession_version>/genome.fna.gz` (BGZF) | **downloaded 2026-10-26**: all 25 assemblies (with GFF3 for the 8 annotated), in Dropbox and on mzebra |
| the panel | P1, `config/genome_panel_P1.tsv` | coordinator repository; read with `ga_common.panel` | available |
| sequence names and lengths of each assembly | as fetched | `provenance.json` in the same folder | with the download |
| Repbase RepeatMasker Edition (final masking only) | 20181026 (`RepBaseRepeatMaskerEdition-20181026.tar.gz`) | held by the lab on utaka; to be copied to mzebra, outside the repository; where it sits in Dropbox is not a concern (owner) | licence not a concern (owner) |

- The registry holds no entries yet (coordinator spec, section 7). Until it does, the
  coordinator's `ga_common.paths` gives the location of each assembly; the project never
  builds a data path by hand.
- The common library is `ga_common`, in the `common/` folder of the coordinator repository (`ga-coordinator`). Install and use it
  as `docs/USING_COMMON.md` in that repository describes: on mzebra, `pip install -e ~/GenomeAlignments/ga-coordinator/common`, or
  pinned to the tag `common-v0.2.1`. It finds the coordinator's configuration through `GA_CONFIG_DIR` or the standard location
  `~/GenomeAlignments/ga-coordinator/config` (needs Python 3.9 or newer).
- Sequence names inside the genome files are exactly as NCBI delivered them (RefSeq `NC_...`
  for annotated assemblies, GenBank otherwise); outputs must use the same names.
- Nothing from the pilot is inherited. The family table `RepeatModeler_TE.xlsx` is recreated
  by this project; the old library is not needed.
- *Observed on mzebra on 2026-10-26, not confirmed as what produced the pilot runs:*
  RepeatMasker 4.1.0 and RepeatModeler 2.0.2a in the conda environment `GenomeSize`;
  TRF, rmblastn and cd-hit were not found on the base PATH. That environment is not used.
- *Chosen (section 6):* the tools run on mzebra from the Dfam TETools 2.01 image (one `.sif` file) under Apptainer 1.5.3
  with squashfuse 0.6.1, in the project's own conda environment `ga-repeats`, which also holds Python 3.11 and `ga_common`
  0.2.1. Apptainer 1.2.2 in conda `base` could not mount images (no squashfuse) and is not used; the `rnaseq` environment
  is not borrowed. The shared account's `.bashrc` sets `APPTAINER_CACHEDIR` for the RNA-seq pipeline; this project
  overrides it per session and does not edit `.bashrc`. CD-HIT 4.8.1 is in the image at `/opt/cd-hit/cd-hit-est` (not on
  its `PATH`).
  A separate copy of the 25 P1 assemblies on utaka (NCBI downloads of 2026-10-20) is not used.

---

## 2a. Relatedness context (coordinator specification, section 4b)

Context, not results and not a design. Sources: Malinsky et al. 2018, Blumer et al. 2025, Kumar et al. 2025;
the species table is `config/species_context.tsv` in the coordinator repository.

| code | species | ecogroup |
|---|---|---|
| MZ | *Maylandia zebra* (reference species) | mbuna |
| AC | *Astatotilapia calliptera* | generalist (stands in for the ancestor) |
| RC | *Rhamphochromis longiceps* | Rhamphochromis (pelagic) |
| DT | *Diplotaxodon limnothrissa* | Diplotaxodon (pelagic) |
| CV | *Copadichromis virginalis* | utaka |
| MC | *Mchenga conophoros* | shallow benthic |
| YH | *Aulonocara* sp. 'chitande type north' | deep benthic (true deep clade) |
| AS | *Aulonocara stuartgranti* | deep benthic |

- **The panel covers all seven ecogroups** of the Lake Malawi radiation, so a library build set can span them.
- **The sequences are very close:** about 0.1–0.25% average divergence between species pairs (SNP level), within a
  radiation that is roughly 0.4–0.8 million years old. That makes one shared library plausible; it is **untested**,
  because none of the three papers measures repeat content or repeat divergence between species.
- The group tree is not strictly bifurcating (extensive gene flow; the benthics may be of hybrid origin), so
  "closely related" is a group-level summary.
- Six large inversions (LG2, 9, 10, 11, 13, 20) distinguish some panel species from the reference. They matter
  for coordinates and alignment, not for masking.
- **Repeats drive the structural differences** (Quah et al. 2025, a pangenome of 8 long-read assemblies): 4.73–9.86% of
  base pairs differ structurally between species, against 0.1–0.25% at SNPs; the graph holds 33.1% more sequence than one
  genome; about 75% of structural-variant sequence is annotated as TE, in a genome that is about 36–38% TE. The paper does not
  test whether the repeat library is complete. See the briefing, Part A.

---

## 3. Outputs

Written under the data root, in a **version folder that is never overwritten**:

```
<data root>/RepeatInformation/versions/RI_v1/
  MANIFEST.yaml            panel version; assemblies used to BUILD the library;
                           assemblies MASKED; command(s); tool versions; date
  library/                 the shared repeat library
  masked/<assembly_id>/    RepeatMasker .out and BED (contract `repeats`); no genome file
```

- A different build subset, or a different library, is a new folder (`RI_v2`, ...). The build
  set and the masked set are recorded separately, so one library can be applied to a different
  subset of genomes.
- Every output names the panel version it was built from (in `MANIFEST.yaml`, and in the
  registry entry once it exists).
- Contract `repeats` (owner: this project): RepeatMasker `.out` + BED, per assembly. Exact
  fields are to be agreed (section 7).
- The **library itself** is an output of this project, but no contract lists it yet (section 7).
- **Question 1 pilot (decision 2026-10-28):** `RI_v1` (library built from MZ_MZ4f), `RI_v2` (MZ_NMetZ), `RI_v3` (YH_YH2f).
  Each `MANIFEST.yaml` also records `purpose: question 1 pilot`, the TETools image tag and its `sha256`, the Apptainer
  and `ga_common` versions and the random seed. Each library masks its own build genome and DT_004f, the held-out genome (six runs common to
  every design); whether each also masks the other pilot genomes is open (section 7). The final library is the next `RI_vN`.
- When libraries are merged, `library/` also holds the **cluster table**: every original family name, its source, its
  cluster, the representative, identity and coverage.
- **Masked genomes are not stored** (owner, 2026-10-30, answered through the coordinator). A soft- or hard-masked genome is
  regenerated on demand from `Genomes/<assembly_id>/<accession_version>/genome.fna.gz` and the BED; the coordinator adds a
  regeneration helper to the common library when a consuming project asks for a form.
- In `masked/<assembly_id>/`: `<assembly_id>.out` and `<assembly_id>.bed`. The BED is a pilot format (BED6, then class/family,
  divergence, RepeatMasker score and id) until the contract's fields are agreed.
- **Comparison results** (which families two libraries share, how similar they are, what the differences mask) are written
  to `~/scratch/ga-repeats/comparisons/<A>_vs_<B>/` until their permanent place is decided (section 7).

---

## 4. Methods

**Partly designed (2026-10-28):** the question 1 pilot below. Questions 1 and 2 are set out in
`docs/repeats_briefing.md` (Parts B and C).
Fixed by the coordinator specification:
- one shared library for all assemblies (decision 2026-10);
- whole genomes, not selected chromosomes (decision 2026-10-26).

**Prior art, recorded and not adopted:** Quah et al. 2025 annotated Lake Malawi genomes with RepeatModeler 2.0.3 and
RepeatMasker 4.1.2-p1 (Dfam TETools container 1.3), de novo libraries built with `-LTRStruct` and combined with
Repbase-derived libraries, `RepeatMasker -e rmblast -no_is`; which genomes the libraries came from is not stated.
*Correction (2026-10-28):* the TETools changelog shows that container 1.3 shipped RepeatModeler 2.0.1 and RepeatMasker
4.1.2; RepeatMasker 4.1.2-p1 came with 1.3.1 and RepeatModeler 2.0.3 with 1.5. The versions the paper names match
TETools 1.5, so its container version and its tool versions do not describe the same image.

Raised, not decided: which genomes build the final library (the pilot below informs it); whether existing family
names are kept now that the family table is recreated; the library policy when a new panel arrives (frozen library
and mask new genomes only, or rebuild and re-mask all).

**Question 1 pilot** (agreed 2026-10-28; libraries built by hand, 2026-10-30 to 2026-11-05; comparisons, masking and
manifests by the scripts of `repeats_009`, so no further step is manual; usage in `docs/pilot_scripts.md`):
1. RepeatModeler with `-LTRStruct` (TETools 2.01) on MZ_MZ4f, MZ_NMetZ and YH_YH2f, one library each (`RI_v1`–`RI_v3`).
   The known-family content RepeatModeler classifies against is recorded and identical for all three runs (Dfam 4.0:
   18,011 RepeatPeps proteins and 30,659 consensus sequences, built into `~/.RepeatModelerCache` during `RI_v1`, then reused).
2. Comparisons: between species (MZ_MZ4f against YH_YH2f) against within species (MZ_MZ4f against MZ_NMetZ).
3. The held-out genome DT_004f is masked with each library.
4. De-duplication and "new family": every family name is prefixed with its source (on a copy,
   `<assembly_id>-families.prefixed.fa`; RepeatModeler's own file is kept as produced); the pooled consensus sequences are
   clustered with cd-hit-est (an identity threshold, and coverage of the shorter sequence; options checked in the
   installed version); one representative per cluster (a classified label is preferred over Unknown; clusters whose
   members carry conflicting classes are flagged); the cluster table is kept. A family is new if its cluster holds no
   member from the other library. The thresholds are set after the pilot; the comparison is reported at several.
5. Reported (briefing B4), per class and in base pairs as well as families: the new families, how much of each genome
   they mask, and DT_004f's masked fraction under each library; plots of genome-wide class composition per genome and
   library (after Quah et al. Fig. 7A, genome-wide half only: the structural-variant half needs another project's calls).

**Final library** (decision 2026-10-28): the de novo library is merged with Repbase at the final masking stage only, with
the same de-duplication; which Repbase taxa are retained is open (section 7). The final build is scripted.

The plan (outputs, formats, checks) is agreed with the owner and recorded here **before any
code is written**.

---

## 4a. Working with a student: how the assistant behaves, and the checkpoint questions

**How the assistant behaves**
1. It treats the user as a student under the owner's supervision (section 0). It asks the question for a gate before
   giving the commands for that gate, and asks the student to paste evidence (the command, the log lines, a version
   or `ls` output) before moving on.
2. If an answer is wrong or vague, it points to the passage of the documentation or the briefing first. After a second
   miss it explains. It never lets a "yes" stand without the evidence.
3. It explains every command before the student runs it: what it does, what the output should look like, how to undo it.
   Scripts take throttling options (`--threads`, `--max-jobs`); long jobs run inside `screen` on the shared server;
   one person works on a repository at a time.
4. It never invents. If it is unsure whether an option exists or what a default is, it says so and has the student
   check the documentation of the **installed** version.
5. It records a two-line summary of the student's answers at the end of each gate, so the owner can review them.
6. Decisions that belong to the owner (the build set, the completeness criteria, tool versions, library sources)
   go into section 7, not into code.

**Checkpoint questions.** The hints are what a good answer contains; the assistant does not show them before the student
has tried.

*Gate 0: concepts, before any software*
1. What does RepeatModeler produce, and what does RepeatMasker do with it? *(A consensus library of families, seed
   alignments and a log; RepeatMasker uses a library to find and annotate or mask copies in a genome.)*
2. What is a TE family, and how does it differ from an insertion? *(A family is the copies descended from one ancestral
   source; an insertion is one copy at one locus. Genomes can share every family and still differ in many insertions.)*
3. Why must the 25 genomes not be concatenated into one database? *(Near-identical loci would appear many times and look
   repetitive.)*
4. In Quah et al., about three-quarters of structural-variant sequence is TE. What does that tell us, and what does it
   not tell us about whether a library is complete? *(Variation is largely TE copies; completeness was not tested, and
   it is not stated which genomes the libraries were built from.)*
5. Give two reasons why "0.1–0.25% divergence" understates how different these genomes are in repeats. *(It is a
   SNP-level figure; structural differences are 4.7–9.9% of base pairs; most extra sequence is private insertions.)*

*Gate 1: inputs*
6. Which file do you use, where is it, and how do you know it is the right accession version of panel P1? *(The Genomes
   folder, `provenance.json`, located through `ga_common.paths`, never a hand-built path.)*
7. The genome is stored BGZF-compressed. Can the tools read it directly, and what do you do if not? *(Check the
   documentation of the installed version; if you decompress, do it in scratch and keep sequence names as received.)*
8. What do the sequence and base counts printed by `BuildDatabase` say about your input, and do they match
   `provenance.json`? *(They must match the provenance record.)*

*Gate 2: running*
9. What does `-LTRStruct` add, and which extra programs does it need? *(Structural discovery of LTR elements; check
   that its dependencies are installed.)*
10. Where do you see the sampling rounds, how many did your run finish, and if it stopped early, why? *(The run log;
    rounds use larger samples; too little long sequence left stops it.)*
11. How many threads, and why? How do you limit them on the shared server, and how does the job survive a lost
    connection? What do you save from the run? *(`-threads`, the throttling flag, `screen`; keep the log and the
    families files.)*
12. Which output file is the library, and which holds the seed alignments?

*Gate 3: reading the results*
13. How many families came out, and what share is classified Unknown (families and base pairs)? What does the part
    after `#` in a family name mean? *(The class and superfamily assigned to the family.)*
14. How would you recognise an LTR element reported as two entries? *(An LTR entry and an internal entry.)*
15. Which TE classes do you expect to be under-detected here, and how will you check each? *(DNA/TIR, SINEs,
    LTR/Unknown, satellites; briefing C1 and C2.)*

*Gate 4: masking*
16. Which RepeatMasker command and options, and does it use only your library or also Dfam or Repbase? How do you
    know? *(The library and species options and the run header.)*
17. What do divergence and class/family in the `.out` file mean, and why do young insertions show low divergence?
    *(Divergence from the consensus; recent copies have had little time to accumulate mutations.)*
18. How many cores does `-pa 8` really use with your search engine? *(Check the documentation of the installed
    version; do not assume.)*

*Gate 5: before registering*
19. How will you decide the library is adequate, and what could you not measure? *(Briefing C2: held-out genome,
    residual repeats, class composition, the known insertions, the unknown fraction.)*
20. What goes in `MANIFEST.yaml`, and why is a version folder never overwritten? *(Panel version, build set, masked set,
    commands, tool and Dfam versions, date; a version is a record of what was done.)*
21. What would make you create `RI_v2` instead of changing `RI_v1`? *(Any change of library, genomes, tool versions
    or commands.)*

---

## 5. Checks (proposed; to be agreed)

- every output validates against the `repeats` contract version (once defined);
- every assembly in the masked set has output, and sequence names and lengths in the output
  match the assembly's `provenance.json`;
- every BED's sequence names are in the assembly's `provenance.json`, and its intervals lie within the recorded lengths, so a
  masked genome can always be regenerated (suggested by the coordinator, 2026-10-30);
- `MANIFEST.yaml` is present, complete, and names the panel version;
- the version folder is new: nothing is overwritten (the sync script refuses to);
- commands and tool versions are recorded so a run can be repeated;
- the completeness measurements of the briefing (Part C2) are reported with what could not be measured, and the unknown
  fraction is reported for the library (families) and for each genome (masked base pairs).

---

## 6. Decisions

| date | decision | reason |
|---|---|---|
| 2026-10 | One shared library for all assemblies | coordinator decision |
| 2026-10-26 | The folder is named `RepeatInformation`, with `versions/RI_vN` | repeat modelling may be run on different subsets of genomes, and every run must stay traceable |
| 2026-10-26 | `RepeatModeler_TE.xlsx` is recreated by this project | the pilot's library location is unknown; the old library is not needed |
| 2026-10-26 | Where the tools run is this project's choice | the owner wants flexibility |
| 2026-10-26 | Whole genomes | coordinator decision |
| 2026-10-26 | The first conversation starts at the briefing (the pangenome paper, then questions 1 and 2); no installation or run before the owner agrees the plan | the owner wants the biology and the two design questions settled before compute |
| 2026-10-26 | The assistant works with a student under the owner's supervision, with checkpoint questions (section 4a) | the owner's request |
| 2026-10-26 | The student uploads PDFs of the reference papers and web pages (briefing, Part D); the assistant confirms it can see items 1–5 before the first conversation goes on | the uploaded copies are always present and fixed to a version; web access is not assumed |
| 2026-10-28 | Question 1 pilot: separate libraries from MZ_MZ4f, MZ_NMetZ and YH_YH2f; between species (MZ_MZ4f vs YH_YH2f) against within species (MZ_MZ4f vs MZ_NMetZ); the held-out genome DT_004f masked with each library | measure whether a second lineage adds families before choosing the build set (briefing B4); DT_004f is pelagic, a subradiation neither library comes from, female and from the same submitter; owner |
| 2026-10-28 | RepeatModeler runs with `-LTRStruct` | LTR elements need structural discovery (briefing C1); owner |
| 2026-10-28 | The tools run on mzebra, from the Dfam TETools 2.01 image (RepeatModeler 2.0.9 and RepeatMasker 4.2.6 by its changelog, confirmed from the tools' own output at setup) under Apptainer 1.2.2 (conda `base`) | fixes for joining DNA-transposon fragments and for simple repeats lost in deep directories; the data root and sync are on mzebra; owner (Apptainer choice: student) |
| 2026-10-28 | Pilot libraries are `RI_v1`–`RI_v3`, each `MANIFEST.yaml` with `purpose: question 1 pilot` | keeps the `RI_vN` convention and every run traceable; owner |
| 2026-10-28 | The "new family" rule is set after the pilot, from the observed differences; the comparison is reported at several thresholds | owner |
| 2026-10-28 | Question 2: the library is built first; completeness is evaluated before the library is accepted | owner |
| 2026-10-28 | The final masking library merges the de novo library with Repbase, added at the final masking stage only | as Quah et al.; owner (placement proposed by the student) |
| 2026-10-28 | De-duplication between libraries: prefix names with their source, cluster with cd-hit-est, one representative per cluster, keep the cluster table (section 4) | family names are not comparable between runs; student |
| 2026-10-28 | The pilot is run by hand; the final build is scripted | quicker to start; scripts once question 1 is settled; student |
| 2026-10-29 | The Repbase licence is not a concern | owner |
| 2026-10-29 | The pilot's measurements are not reported to the coordinator; each version's `MANIFEST.yaml` records what was done | owner |
| 2026-10-29 | The runtime is Apptainer 1.5.3 with squashfuse 0.6.1 in the project's own conda environment `ga-repeats`, replacing 1.2.2 in `base` | 1.2.2 in `base` could not mount images (no squashfuse); 1.5.3 with squashfuse 0.6.1 was shown to run containers on this account; an own environment leaves `base` and the RNA-seq environment untouched; student |
| 2026-10-30 | `ga_common` 0.2.1 is installed in `ga-repeats`, with Python 3.11 | it was not installed in `ga-coord`; projects install their own copy (section 2); student |
| 2026-10-30 | Masked genomes are not stored: `masked/<assembly_id>/` holds the `.out` and the BED only; a masked genome is regenerated on demand from the genome and the BED | answer to project 1's blocking item, through the coordinator; owner |
| 2026-10-30 | Where the Repbase file sits in Dropbox is not a concern | answer to project 1's blocking item, through the coordinator; owner |
| 2026-11-05 | No manual steps: work runs in scratch (`~/scratch/ga-repeats`), and scripts move finished results into the version folder and write each `MANIFEST.yaml`; the comparisons, masking and manifests use the scripts of `repeats_009` (replacing "the pilot is run by hand" for these steps) | the owner wants no manual operations; building in scratch and moving into place only after checks is the coordinator's own pattern (coordinator specification, section 5) and keeps temporary files out of Dropbox; owner |
| 2026-11-05 | Long jobs run in `screen`, at 16 threads each, under `nice -n 10` | three runs use half of mzebra's 96 cores; a lost connection does not stop a run; student |
| 2026-11-06 | The student applies zips on the laptop and pushes to GitHub; mzebra only pulls (`git pull --ff-only`) | GitHub stays the backup and the server is never ahead of it; the shared server account does not commit; owner |

---

## 7. Open questions

| question | blocking |
|---|---|
| Contract for the shared library, and exact fields of the `repeats` contract (class and family naming, divergence, BED columns) | before first integration |
| Keep the existing family names, or take the names the recreated library gives? | no |
| **Question 1:** which genomes build the final library? The pilot is agreed (section 6); open are the "new family" thresholds (set after the pilot) and the build set that follows from the result | **yes, before the final library build** |
| **Question 2:** what completeness measurements, per TE class, and what counts as good enough (briefing, Part C)? Includes the unknown fraction and DNA, SINE and LTR/Unknown families | **yes, before the library is accepted** |
| Repbase merge (decided, section 6): which Repbase taxa are retained? The owner and the student decide | before the final masking |
| **Satellites and other tandem repeats:** does project 1 build and annotate them, or project 2b (centromeres, arrays)? (relayed to the coordinator) | before the library is accepted |
| Library policy for new panels: frozen or rebuilt? | before P2 |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) masked too? | no |
| Pilot masked set: the six runs common to every design (each library on its own genome and on DT_004f) may run. Are the cross-genome runs all six more (12 in total, symmetric), or only those briefing B4 strictly needs (7 or 8 in total)? | before the cross-genome masking |
| For the coordinator: AS_Sanger (`fAulStu2.1`) and AS_SAMEA123396442 (`fAulStu2.pri`) carry the same specimen name. Same individual? | before AS is used |
| Permanent place for comparison results: a `comparisons/` folder under `RepeatInformation` would need the coordinator to add it to `ga_common.paths` and `sync.py`; until then they stay in scratch | before the pilot is reported |

---

## 8. Status

First conversation done (2026-10-28); the owner's decisions are in section 6. The 25 P1 genomes are confirmed in the data
root on mzebra. Container set up: TETools 2.01 (`sha256` 8b805d16d279eb57e4948719d2acd8dcd75b31b5ce22fbadcd3a522307fe97ca) under
Apptainer 1.5.3 in `ga-repeats`; RepeatModeler 2.0.9 and RepeatMasker 4.2.6 confirmed from the tools' own output.

**Libraries built** (copied into their version folders 2026-11-05; manifests not yet written):

| version | built from | families | seed | run time | BuildDatabase |
|---|---|---|---|---|---|
| `RI_v1` | MZ_MZ4f | 2,396 (1,787 RepeatScout/RECON + 848 LTR pipeline − 239 redundant) | 1793333684 | 36:14:33 | 61 sequences, 972,958,384 bp |
| `RI_v2` | MZ_NMetZ | 2,398 (1,782 + 862 − 246) | 1793758093 | 33:53:39 | 75 sequences, 959,220,039 bp |
| `RI_v3` | YH_YH2f | 2,306 (1,773 + 693 − 160) | 1793758241 | 38:42:38 | 69 sequences, 976,509,173 bp |

Each ran five rounds with 16 threads and reached about 42% database coverage by the log's final coverage line; all three
were classified against the same Dfam 4.0 reference.
**Next:** install `repeats_009`; the comparisons (`RI_v1` with `RI_v2`, `RI_v1` with `RI_v3`); the six masking runs common to
every design; the owner's decision on the cross-genome runs; the summaries; the manifests; `sync.py push`. Nothing is ready
to integrate.

**Gate log** (section 4a, item 5):
- Gate 0, 2026-10-28: all five answered. Q2's definitions and Q3 were right first time; Q5 needed one pointer and the
  within-species part of Q2 a second prompt; Q1 (the run log; labels come from the library) and Q4 (why the SV figure
  cannot show completeness) were explained after two tries.
- Gates 1 and 2, 2026-11-05: question 6 explained (`ga_common` was not installed for this project; installed 2026-10-30);
  questions 7, 10 and 11 explained after going unanswered; question 8's counts recorded (61, 75 and 69 sequences), with the
  check against each genome run by `write_manifest.py`. Questions 9 and 12 were not asked separately: each log confirms
  `-LTRStruct`'s tools, and the families file, seed alignments and logs were kept.
- Gates 3 and 4, open: question 13 asked, not yet answered; question 16 asked with `repeats_009`; question 18 settled by
  RepeatMasker's own help (each `-pa` job uses 4 cores with RMBlast), which the scripts follow.

---

## 9. Environment and delivery (from coordinator specification, section 3a)

- Repository `ga-repeats` (https://github.com/ptmcgrat/ga-repeats): `~/GenomeAlignments/ga-repeats` on the laptop and on mzebra
  (outside Dropbox). One person works on it at a time (the server account is shared).
- **Code is delivered as numbered zips** (`repeats_001.zip`, ...) with `MANIFEST.txt`,
  `INSTRUCTIONS.md` and `files/`. On the laptop: `git pull`, unzip the exact file, `rsync`
  `files/` over the repository (excluding `.git` and `config/local.yaml`), `git rm` anything
  under DELETE, check `VERSION`, commit, push. On the server: `git pull --ff-only`. `repeats_001`–`repeats_008` were the coordinator's
  specification drafts; code starts at `repeats_009`. The student applies zips on the laptop and pushes (owner, 2026-11-06).
- **Data root:** Dropbox `Apps/GenomeAlignments`, mirrored on mzebra at
  `~/Temp/GenomeAlignments`. Move version folders with the coordinator's
  `scripts/sync.py` (dry run unless `--go`; never overwrites), for example
  `python ~/GenomeAlignments/ga-coordinator/scripts/sync.py push RepeatInformation/versions/RI_v1`.
- **Throttling:** every script takes `--threads` and `--max-jobs` (`GA_MAX_THREADS`,
  `GA_MAX_JOBS`); no limit unless asked.
- Install commands are printed in the conversation with every zip, and specifications are
  shared directly as files readable in the browser.
- Python 3.9 is the minimum for code that uses the common library.
- **Container and Repbase:** the TETools image is one `.sif` file (`~/containers/dfam-tetools-2.01.sif`), outside the data
  root and the repository; its tag and `sha256` go into every `MANIFEST.yaml`. Project 1's tools run in the conda
  environment `ga-repeats`. The Repbase files stay out of the repository; where they sit in Dropbox is not a concern (owner).
- **Documents in the Claude project:** this specification, `docs/repeats_briefing.md`, `genome_panel_P1.tsv` and `species_context.tsv`
  from the coordinator repository, the RepeatModeler and RepeatMasker documentation, the TETools changelog (release 2.01), and the papers and pages listed in the
  briefing, Part D, **uploaded by the student as PDFs** (items 1–5 required, 6–8 may follow).
