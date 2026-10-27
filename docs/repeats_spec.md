# Repeats specification

**Project:** 1, Repeats
**Repository:** `ga-repeats` *(name assumed)*
**Document version:** 0.5, 2026-10-26 (0.1: first draft from coordinator specification v0.4; 0.2: delivery text from v0.5; 0.3: relatedness context from v0.7; 0.4: the P1 genomes are downloaded, from v0.8; 0.5: the first conversation starts at the briefing; working with a student)
**Owner:** McGrath lab
**Status:** draft for discussion with the owner; no code written; **start at section 0 and `docs/repeats_briefing.md`**

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

- The registry holds no entries yet (coordinator spec, section 7). Until it does, the
  coordinator's `ga_common.paths` gives the location of each assembly; the project never
  builds a data path by hand.
- The common library is installed pinned:
  `pip install "git+https://github.com/ptmcgrat/ga-coordinator.git@common-v0.1.0#subdirectory=common"`
  (needs Python 3.9 or newer).
- Sequence names inside the genome files are exactly as NCBI delivered them (RefSeq `NC_...`
  for annotated assemblies, GenBank otherwise); outputs must use the same names.
- Nothing from the pilot is inherited. The family table `RepeatModeler_TE.xlsx` is recreated
  by this project; the old library is not needed.
- *Observed on mzebra on 2026-10-26, not confirmed as what produced the pilot runs:*
  RepeatMasker 4.1.0 and RepeatModeler 2.0.2a in the conda environment `GenomeSize`;
  TRF, rmblastn and cd-hit were not found on the base PATH. Where the tools run is the
  project's choice (the server is not assumed).

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
  masked/<assembly_id>/    RepeatMasker .out, BED, masked genome
```

- A different build subset, or a different library, is a new folder (`RI_v2`, ...). The build
  set and the masked set are recorded separately, so one library can be applied to a different
  subset of genomes.
- Every output names the panel version it was built from (in `MANIFEST.yaml`, and in the
  registry entry once it exists).
- Contract `repeats` (owner: this project): RepeatMasker `.out` + BED, per assembly. Exact
  fields are to be agreed (section 7).
- The **library itself** is an output of this project, but no contract lists it yet (section 7).

---

## 4. Methods

**Not designed yet.** The two questions to settle first are in `docs/repeats_briefing.md` (Parts B and C).
Fixed by the coordinator specification:
- one shared library for all assemblies (decision 2026-10);
- whole genomes, not selected chromosomes (decision 2026-10-26).

**Prior art, recorded and not adopted:** Quah et al. 2025 annotated Lake Malawi genomes with RepeatModeler 2.0.3 and
RepeatMasker 4.1.2-p1 (Dfam TETools container 1.3), de novo libraries built with `-LTRStruct` and combined with
Repbase-derived libraries, `RepeatMasker -e rmblast -no_is`; which genomes the libraries came from is not stated.

Raised, not decided: how the library is built (which genomes, for example one or more per ecogroup;
RepeatModeler with `-LTRStruct`?); whether existing family names are kept now that the family table is
recreated; the library policy when a new panel arrives (frozen library and mask new genomes
only, or rebuild and re-mask all).

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
- `MANIFEST.yaml` is present, complete, and names the panel version;
- the version folder is new: nothing is overwritten (the sync script refuses to);
- commands and tool versions are recorded so a run can be repeated;
- the completeness measurements of the briefing (Part C2) are reported with what could not be measured, and the unknown
  fraction is reported for the library and for the masked genome.

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

---

## 7. Open questions

| question | blocking |
|---|---|
| Contract for the shared library, and exact fields of the `repeats` contract (class and family naming, divergence, BED columns) | before first integration |
| Keep the existing family names, or take the names the recreated library gives? | no |
| **Question 1:** which genomes build the library (one, one per species, one per ecogroup, all)? The student proposes a comparison and a rule for "new family"; the owner agrees (briefing, Part B) | **yes, before any library build** |
| **Question 2:** what completeness measurements, per TE class, and what counts as good enough (briefing, Part C)? Includes the unknown fraction and DNA, SINE and LTR/Unknown families | **yes, before the library is accepted** |
| RepeatModeler with `-LTRStruct`, and which tool versions (mzebra has RepeatModeler 2.0.2a and RepeatMasker 4.1.0; newer releases and Dfam 4.0 exist)? | before the library is built |
| Combine the de novo library with Dfam or Repbase (Repbase has licence terms)? | no |
| **Satellites and other tandem repeats:** does project 1 build and annotate them, or project 2b (centromeres, arrays)? (relayed to the coordinator) | before the library is accepted |
| Library policy for new panels: frozen or rebuilt? | before P2 |
| Where the tools run, and how long jobs are kept alive on the shared server (no scheduler; `tmux`, `screen` and `nohup` not checked) | no |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) masked too? | no |
| Whether masked genomes are stored (size) or regenerated from the `.out` | no |

---

## 8. Status

Not started. Specification and briefing drafted; the P1 genomes are available. The first conversation starts at
section 0; the student's proposals for questions 1 and 2 are awaited.
Nothing is ready to integrate.

---

## 9. Environment and delivery (from coordinator specification, section 3a)

- Repository `ga-repeats`: `~/GenomeAlignments/ga-repeats` on the laptop and on mzebra
  (outside Dropbox). One person works on it at a time (the server account is shared).
- **Code is delivered as numbered zips** (`repeats_001.zip`, ...) with `MANIFEST.txt`,
  `INSTRUCTIONS.md` and `files/`. On the laptop: `git pull`, unzip the exact file, `rsync`
  `files/` over the repository (excluding `.git` and `config/local.yaml`), `git rm` anything
  under DELETE, check `VERSION`, commit, push. On the server: `git pull`.
- **Data root:** Dropbox `Apps/GenomeAlignments`, mirrored on mzebra at
  `~/Temp/GenomeAlignments`. Move version folders with the coordinator's
  `scripts/sync.py` (dry run unless `--go`; never overwrites), for example
  `python ~/GenomeAlignments/ga-coordinator/scripts/sync.py push RepeatInformation/versions/RI_v1`.
- **Throttling:** every script takes `--threads` and `--max-jobs` (`GA_MAX_THREADS`,
  `GA_MAX_JOBS`); no limit unless asked.
- Install commands are printed in the conversation with every zip, and specifications are
  shared directly as files readable in the browser.
- Python 3.9 is the minimum for code that uses the common library.
- **Documents in the Claude project:** this specification, `docs/repeats_briefing.md`, `genome_panel_P1.tsv` and `species_context.tsv`
  from the coordinator repository, the RepeatModeler and RepeatMasker documentation, and the papers listed in the
  briefing, Part D (web pages saved as PDF).
