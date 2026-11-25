# Repeats specification

**Project:** 1, Repeats
**Repository:** `ga-repeats` (https://github.com/ptmcgrat/ga-repeats)
**Document version:** 0.25, 2026-11-24 (0.1: first draft from coordinator specification v0.4; 0.2: delivery text from v0.5; 0.3: relatedness context from v0.7; 0.4: the P1 genomes are downloaded, from v0.8; 0.5: the first conversation starts at the briefing; working with a student; 0.6: uploading wording corrected; 0.7: the student uploads the PDFs; 0.8: using the common library; 0.9: first conversation done; the owner's decisions on the question 1 pilot, the tools and Repbase; prior-art versions corrected; 0.10: masked genomes are not stored (coordinator answer); Repbase licence and location are not a concern; runtime and `ga_common` in the `ga-repeats` environment; first pilot run started; 0.11: the three pilot libraries built; comparisons, masking and manifests scripted (`repeats_009`); Dfam 4.0 classification reference; the student pushes from the laptop; 0.12: the pilot's results and question 1's answer (one library per ecogroup); the merge rehearsals `RI_v4` and `RI_v5` and the merge settings; the owner's delegation; accuracy checks begun; 0.13: the coding-sequence check's result, the gene filter, Repbase in place with its metadata; 0.14: the gene filter's host proteins cleaned after its first report; the Repbase label check; 0.15: the gene filter's evidence refined after its second report; 0.16: gene-derived families tagged, not removed; panEDTA's copy-number filter and inconsistency measure adopted; the `RI_v7` rehearsal; 0.17: the copy-number filter set aside for now; the concatenation test; 0.18: the verification panel for question 2; 0.19: merge kept; tags by coverage only; 0.20: the first verification results; the known insertions anchored on published sequence; 0.21: label checks against the curated library; the findings document; 0.22: existing curation tools installed for comparison; 0.23: every evaluation against published data, like for like; the reference genome fAstCal1.2; 0.24: the curation pilots; 0.25: curation tools dropped; the two strategies with Repbase on fAstCal1.2)
**Owner:** McGrath lab
**Status:** question 1 answered; merge settings decided; `RI_v6` written (178 families tagged); the question 2 checks found class labels to be the main problem; the curation tools were tried and dropped for now; every evaluation is now against published data, like for like, on Quah et al.'s assembly fAstCal1.2, where our two library strategies, each combined with Repbase as Quah et al. did, are being compared; established findings in `docs/repeats_findings.md`; **section 8 says where things stand**

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
| Repbase RepeatMasker Edition (final masking; also the accuracy checks) | 20181026: `RepBaseRepeatMaskerEdition-20181026.tar.gz` (sha256 `213a15f9dafbd1302b1d544b82b9141eacab1c5be5b36233da7446f901375b2e`), holding the sequences only (`RMRBSeqs.embl`, 49,011 entries) | on mzebra in `~/Repbase/20181026`, outside the repository (copied from utaka 2026-11-13); where it sits in Dropbox is not a concern (owner) | licence not a concern (owner) |
| Repbase metadata (classes and species of those entries) | `RMRBMeta.embl`, release 20181026, 49,011 entries; shipped with RepeatMasker 4.1.0 (copied from the `GenomeSize` environment on mzebra) | `~/Repbase/20181026/RMRBMeta.embl` | CC BY 4.0; `repbase_match.py` joins it to the sequences by entry name and records both files' sha256 |
| Transposon proteins (gene filter) | Dfam 4.0 RepeatPeps, 18,011 proteins, as RepeatModeler cached it | `~/.RepeatModelerCache/PROT-Dfam_4.0/repeatpeps.lib` (the TETools image carries no libraries) | `gene_filter.py` records its sha256 |

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
  and `ga_common` versions and the random seed. Each library masks its own build genome and DT_004f, the held-out genome;
  two cross runs (`RI_v1` on YH_YH2f, `RI_v3` on MZ_MZ4f) separate lineage from run (section 6). The final library is the next `RI_vN`.
- A **merged version** (`merge_libraries.py`) has no single build genome. Its `library/` holds the merged library
  `<version>-families.prefixed.fa`, the **cluster table** `cluster_table.tsv` (every original family name, its source, its
  cluster, the representative, identity, how much of the family and of the representative the alignment covers, the
  family's own label, the cluster's label and the rule that set it) and `merge_record.json` (the sources and their library
  sha256, thresholds, command and counts).
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

**Merging** (decided 2026-11-13, after the rehearsals `RI_v4` and `RI_v5`; `merge_libraries.py`): the libraries' prefixed
families are pooled and clustered with cd-hit-est at 0.80 identity with local alignment (`-c 0.80 -G 0`), the alignment
covering at least 80% of both the shorter and the longer sequence (`-aS 0.80 -aL 0.80`); cd-hit-est's representative (the
longest member) is kept. Each cluster's label: (a) no classified member: Unknown; (b) classified members of one class: that
class, and Unknown members inherit it; (c) disagreeing members: the class given by the most source libraries; (d) a tie:
the representative's class if it is among the tied, else that of the longest classified member among them. The full label
(class/superfamily) is the representative's when its class is chosen, else the most frequent among the members of that
class. Clusters decided by (c) or (d) are flagged and inspected by hand. Requiring coverage of the longer sequence keeps a
family contained in a longer one (a nested insertion, or a chimeric consensus) as a family of its own: without it
(`RI_v4`), 30 classified families were relabelled by much shorter members they contained.

**Copy-number filter** (decided 2026-11-17, after panEDTA, Ou et al. 2024; `merge_libraries.py
--min-full-length-copies 3`): before pooling, a family is kept only if it has at least three full-length copies in its
own build genome. A copy is full-length by the 80-80-80 rule: its alignments cover at least 80% of the consensus, at
at least 80% identity, over at least 80 bp; RepeatMasker's fragments of one copy are joined by its element ID. Copies
are counted in each source's masking of its own build genome; every family's count is kept (`copy_number_table.tsv`). **Set aside for now (2026-11-17):** in the `RI_v7` rehearsal it dropped 1,290 of 7,100 families (18%), among them large LINE, LTR and DNA families that rarely occur full-length (40% of families of 8 kb or more were dropped); the dropped families accounted for 4.3-5.2% of DT_004f in their own libraries' masking (an upper bound). The option stays in `merge_libraries.py`, unused.

**Concatenation test** (2026-11-17; `merge_libraries.py --no-clustering`, `compare_concatenation.py`): the alternative to merging is to concatenate the libraries and let RepeatMasker's scores choose among near-duplicate families at each copy. `RI_v8` concatenates `RI_v1`-`RI_v3` (7,100 families, no clustering) and masks DT_004f; it is compared with `RI_v5` (the same sources merged) on coverage, class shares, label agreement, fragmentation (annotations and elements per Mb), how each `RI_v5` cluster's copies are split among its members, and run time. Concatenation is adopted if coverage rises clearly and consistency does not suffer; a middle ground is to collapse only near-identical families (95% identity or more).

**Gene filter** (decided 2026-11-13; `gene_filter.py`, after `gene_overlap.py` and `repbase_match.py`): de novo libraries
include families built from multi-copy host genes (NLRs, TRIMs, immune and olfactory receptors, zinc fingers, histones,
tubulins); `RI_v5` masked 13–14% of the coding sequence of three annotated genomes. Evidence A: on every annotated genome
checked, a family is the main masker of at least 5 genes of known multi-copy gene families (half or more of their coding
sequence masked), or masks coding sequence in at least 20 named genes ("uncharacterized" gene models, mostly built from
transposons, are not counted). Evidence B: blastx of the family's consensus matches a
host protein at evalue 1e-10 or below, with at least twice the bitscore of its best transposon-protein match
(RepeatPeps). The host proteins are the longest protein per gene of MZ_NMK2024a, translated from its genome and
annotation, named genes only: products called "uncharacterized" (mostly gene models built from transposons), products
naming a transposon protein, and proteins that blastp matches to RepeatPeps (domesticated transposon genes such as
ZBED or NYNRIN) are left out. How much of the family the host match covers is recorded. (The first report, with
uncharacterized and transposon-derived host proteins kept, flagged 936 families by evidence B, most of them real
transposons; nothing was written from it.) Evidence B counts only when it is strong: the host match covers at least 30%
of the family (decided 2026-11-17: coverage only; until then it also counted when the host gene belonged to a known
multi-copy gene family, whatever the coverage, which tagged 126 more families on matches covering a median 8% of them;
that rule stays available as `--b-family-exception`). The known gene families (NLRs, TRIMs, immune, olfactory, taste and
trace amine receptors, protocadherins, zinc fingers, histones, tubulins) still name the tag categories. A family with only weak evidence B (a small
host-like stretch in what is otherwise a transposon) is kept and listed; cutting the stretch out of its consensus is a
possible later step. Before any version is written, the report estimates what the removals cost: the masked coding
sequence and the genome masking the removed families account for. Families flagged by both are tagged; families flagged by one are tagged unless a curated Repbase transposon covers at least half of them (decided 2026-11-17: tagged, not removed). A tagged family keeps its sequence and class and is renamed `<family>_hostgene-<CATEGORY>#<class>`, for example `MZ_MZ4f_rnd-1_family-163_hostgene-NLR#LTR/Gypsy`; CATEGORY is NLR, TRIM, IMMUNE, OR, TAAR, TASTE, PCDH, ZNF, HIST, TUB, OTHER (another named host gene) or MULTI (no protein match, but coding sequence in many named genes), from the best host protein, else from the gene family it mainly masks. Masking and class shares are unchanged by tagging, so they stay comparable with Quah et al., who kept such families in their library; any analysis can include or exclude tagged families. The thresholds are reviewed on the first report before any version is written.

**Final library**: built from one assembly per ecogroup, seven libraries in total, merged as above (decision 2026-11-11);
merged with Repbase at the final masking stage only, with the same de-duplication (decision 2026-10-28). Which assemblies,
and which Repbase taxa, are open (section 7). The final build is scripted.

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

**The principle for every evaluation** (decided 2026-11-23, section 6): every evaluation uses published, independent data
as its reference, compared like for like (the same assembly, the same definitions); our own libraries are only ever what
is evaluated, never the reference. Published references differ in strength: those resting on independent evidence (the
insertions of Carleton et al. 2020, confirmed by sequencing; Repbase's curated cichlid families) weigh more than the
output of another annotation pipeline (the percentages of Quah et al. 2025), which shows consistency with a published
method rather than correctness. Exception: controls whose truth is known by construction (a shuffled chromosome, where
any masking is spurious).

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
| 2026-11-09 | Cross-genome masking: `RI_v1` on YH_YH2f and `RI_v3` on MZ_MZ4f only, not all six cross runs | enough to tell whether YH_YH2f's lower LTR total is lineage or run; owner and student |
| 2026-11-11 | **Question 1:** the final library is built from one assembly per ecogroup, seven in total, RepeatModeler run on each | the pilot (section 8): lineage does not limit what a library masks, but each further library adds masked sequence, with diminishing returns (+2.41, then +1.16 points of DT_004f); owner and student |
| 2026-11-12 | Merge rehearsal `RI_v4`: `RI_v1`–`RI_v3` merged (0.80 identity, 80% of the shorter sequence, local alignment, label rule a–d), DT_004f masked | test the merge before seven libraries depend on it; student |
| 2026-11-13 | The owner leaves the project's decisions to the student | stated by the student; whether this covers the decisions this document reserves to the owner by name is to be confirmed (section 7) |
| 2026-11-13 | Merge rehearsal `RI_v5`: as `RI_v4`, plus 80% coverage of the longer sequence (`-aL 0.80`) | in `RI_v4`, short families contained in long ones (nested insertions, chimeric consensus sequences) joined their clusters and relabelled them; student |
| 2026-11-13 | Merge settings for the final build: those of `RI_v5` (section 4, Merging) | `RI_v5`: conflicting clusters 132 → 16, no classified family relabelled, labels agree with each source on 86–87% of shared sequence, 99.4% of the sources' masking kept (43.91% of DT_004f, against 44.15%); 5,288 families against 3,849; student |
| 2026-11-13 | Repbase is kept on mzebra in `~/Repbase/20181026`, with the RepeatMasker metadata file (`RMRBMeta.embl`) beside it | the package holds sequences only; the matching release's classes come from RepeatMasker 4.1.0's metadata, already on mzebra; student |
| 2026-11-13 | Gene filter (section 4), rehearsed on `RI_v5` as `RI_v6`, before the known insertions and the Repbase label check | `RI_v5` masks 13–14% of coding sequence, 3.3–3.7% in known multi-copy gene families; the later checks also test that the filter keeps real transposons; student |
| 2026-11-14 | Gene filter refined (section 4): evidence A counts named genes only; evidence B counts only when strong (host match covering 30% of the family, or a known multi-copy gene family); weak evidence B alone is kept and listed; the report estimates the cost | the second report: the host match covered a median 11% of the families removed on protein evidence alone (105 of them LTR elements), and 98 of evidence A's 132 families came only through a gene count that included uncharacterized models; student |
| 2026-11-17 | Gene-derived families are tagged, not removed: `<family>_hostgene-<CATEGORY>#<class>`, categories in section 4 | Quah et al. kept them in their library (they filtered their gene set instead); tagging keeps masking and class shares comparable while letting analyses exclude them; the class after `#` stays a transposon class because RepeatMasker reads it as the classification; student |
| 2026-11-17 | panEDTA's copy-number filter before merging (section 4): at least 3 full-length (80-80-80) copies in the family's own build genome | removes low-copy and incomplete families that would accumulate across seven libraries, and keeps the library and masking time down (Ou et al. 2024); student |
| 2026-11-17 | panEDTA's annotation-inconsistency measure for the labels (question 2) | a reference-free check of label consistency (Ou et al. 2024); to be built after reading EDTA's own script; student |
| 2026-11-17 | Rehearsal `RI_v7`: `RI_v1`–`RI_v3` merged with the copy-number filter and the merge settings, then the gene check, the Repbase match and tagging | test the copy-number filter and the tagging before the final build; student |
| 2026-11-17 | The copy-number filter is set aside for now | in `RI_v7` it dropped 18% of families, among them large LINE, LTR and DNA families with few full-length copies; the dropped families accounted for 4.3-5.2% of DT_004f in their own libraries' masking; student |
| 2026-11-17 | Concatenation test `RI_v8`: `RI_v1`-`RI_v3` concatenated without clustering, DT_004f masked, compared with `RI_v5` | concatenation loses no family, but near-duplicates may split one family's copies, add label noise and fragment annotations, and masking takes longer; test before choosing; student |
| 2026-11-17 | Question 2 verification panel: the known insertions (*rx1* and *mitfa* first; *tbx2a* and *fhl2b* once their genes are identified), lineage-specific loss on the build genomes, masking per chromosome with LG3, an approximate divergence landscape, curated-family recovery against Repbase's cichlid entries (Flynn et al. 2020's categories), the inconsistency measure, residual repeats on the final library, the same-genome comparison with Quah et al. | covers briefing C2 and adds the curated-family and lineage checks; the first five run on `RI_v5` now (`repeats_023`); student |
| 2026-11-17 | How Repbase joins the final library (taxa, and precedence in clusters) is decided after the verification | the checks show where the de novo library is weak; student |
| 2026-11-17 | The libraries are merged (section 4), not concatenated | the `RI_v8` test: concatenation gained 0.09 points of DT_004f but split one family's copies across several names in 980 of 1,150 multi-member clusters (about 164 Mb), with labels and fragmentation unchanged; student |
| 2026-11-17 | Evidence B counts by coverage only: strong when the host match covers at least 30% of the family | with the gene-family exception, 126 more families were tagged on matches covering a median 8% of them, mostly LTR and LINE families, likely transposons carrying a gene-like stretch; a tag tells other projects "not a transposon"; student |
| 2026-11-18 | Established findings are kept in `docs/repeats_findings.md`, apart from this specification and from hypotheses, and updated whenever a result becomes solid | to frame a paper; student |
| 2026-11-18 | `-LTRStruct` stays (decided 2026-10-28); mislabelled or composite consensus sequences are dealt with after RepeatModeler | the mislabelling is not specific to the LTR pipeline (findings B11: 42 of 54 disagreements come from RepeatScout/RECON families), and RepeatModeler 2.0.9 offers no option capping the LTR internal length (its LTR options are `-LTRStruct` and `-ltr_retriever_dir`); student |
| 2026-11-19 | Existing tools come before our own rules: the curation tools MCHelper (Orozco-Arias et al. 2024) and TEtrimmer (Qian et al. 2025) are run on a source library and compared, with the checks of question 2, against the raw library and later against rule-based label correction | the label problem (findings B13, B14) concerns every family, so a fix must be systematic; published, benchmarked tools first; student |
| 2026-11-20 to 2026-11-23 | The curation tools run in their own conda environments, called from `ga-repeats` through `conda run`; `ga-repeats` is unchanged: `repeat_family_merging_mchelper` (MCHelper 1.7.0, commit cd80885, its authors' exact environment `MCHelper.yml`: Python 3.8.16, BLAST 2.10.1, HMMER 3.3.2) and `repeat_family_merging_tetrimmer` (TEtrimmer 1.7.2 from bioconda: Python 3.10.14, RepeatModeler 2.0.8, RepeatMasker 4.2.3, RMBlast 2.14.1, Perl 5.32.1). Tool sources, databases and test runs are under `<scratch>/tools` and `<scratch>/tool_tests` | one environment cannot be solved: TEtrimmer's `pfam_scan` needs zlib below 1.3, RMBlast 2.17 (needed by RepeatMasker 4.2.4 and RepeatModeler 2.0.9) needs zlib 1.3.2 or later, and an unconstrained solve fell back to RepeatModeler 1.0.8 with a Perl that does not run on mzebra; how they were built: `docs/tool_environments.md`; student |
| 2026-11-23 | Every evaluation uses published, independent data as its reference, compared like for like (the same assembly, the same definitions); our libraries are only ever what is evaluated. Exception: controls whose truth is known by construction | the masked shares of our libraries (40.3-43.9% of DT_004f) are above every published figure for Lake Malawi cichlids (35-37% from long-read assemblies, Conte et al. 2019; 36.32-38.29%, Quah et al. 2025), and only comparisons with published data can say whether that is an excess or an error; student |
| 2026-11-23 | Comparisons with Quah et al. 2025 use their assembly, *A. calliptera* fAstCal1.2 (GCF_900246225.1), not AC_fAstCal68 (GCA_964374335.1, a different, more contiguous assembly): it is a reference genome, downloaded into `<scratch>/references/AC_fAstCal1_2/` only, listed in `references/reference_genomes.tsv`, used only for comparisons with published TE annotations of it, never to build or curate a library nor for alignment, and not added to the shared genome store (no request to the coordinator). Comparisons of library variants are masked on it from now on; results on DT_004f stay as recorded | Quah et al.'s TE figures (36.55%: DNA 12.09%, LINE 8.32%, LTR 7.38%, Unknown 8.54%, SINE, helitrons and retroposons <1% together) are for fAstCal1.2, so a comparison on another assembly mixes assembly and library; student |
| 2026-11-24 | The curation tools MCHelper and TEtrimmer are dropped for now; their environments and scripts stay | in the pilots (section 8) neither fixed the composite family at *tbx2a* (769); MCHelper gives a class to every family it cannot classify, and called the satellite family 9801 a Helitron; TEtrimmer would need about 10 days per library on 16 threads; student |
| 2026-11-24 | Two library strategies are compared on fAstCal1.2, each combined with Repbase as Quah et al. 2025 combined their de novo libraries with Repbase's RepeatMasker libraries: the three libraries concatenated without clustering (`RI_v8`, Quah et al.'s strategy) plus Repbase = `RI_v9`, and the merged, tagged library (`RI_v6`, our strategy) plus Repbase = `RI_v10`. Repbase is added unclustered (`add_to_library.py`), so RepeatMasker keeps the best-scoring family at each locus | to follow Quah et al.'s procedure, so that our results can be set beside theirs, and because a Repbase consensus that scores better may replace a mislabelled family of ours at a locus; student |
| 2026-11-24 | The Repbase entries added are those of Repbase 20181026 assigned to *A. calliptera* or to one of its ancestral clades (RepeatMasker's `-species` rule; sister taxa such as *O. niloticus* and *Danio rerio* left out), the lineage taken from RepeatMasker 4.1.0's own taxonomy (`queryTaxonomyDatabase.pl`, run read-only in the `GenomeSize` environment; `taxonomy.dat` sha256 6840a4813a78684cb215e95c90e21cd343d82f420b285d1bdd38dc8f1d7a89f0) | Quah et al. do not say which Repbase entries they used; this is the usual way of deriving a species' RepeatMasker library from Repbase; all 49,011 entries span 1,067 species and clades, mostly plants and invertebrates; proposed by the assistant, used by the student |
| 2026-11-24 | Ensembl's soft-masked fAstCal1.2 (release 103) is a published reference for which bases are repeats (it carries no classes), fetched into the reference folder and checked against Ensembl's checksums | NCBI provides no RepeatMasker run for fAstCal1.2; Ensembl's masking is an independent published annotation of the same assembly; student |

---

## 7. Open questions

| question | blocking |
|---|---|
| Contract for the shared library, and exact fields of the `repeats` contract (class and family naming, divergence, BED columns) | before first integration |
| Keep the existing family names, or take the names the recreated library gives? | no |
| **Question 1:** answered (section 6, 2026-11-11): one assembly per ecogroup. Open: which assembly in each ecogroup; whether `RI_v1` (mbuna) and `RI_v3` (deep benthic) serve as two of the seven; and which genome is held out to check the final library if the Diplotaxodon library is built from DT_004f. The pilot's "new family" thresholds changed none of its conclusions (0.80–0.95), so no separate rule is set | **yes, before the final library build** |
| **Question 2:** what completeness measurements, per TE class, and what counts as good enough (briefing, Part C)? Includes the unknown fraction and DNA, SINE and LTR/Unknown families. Begun on `RI_v5` (section 8): coding sequence masked on DT_005m; the known insertions (*rx1*, *tbx2a*, *mitfa*, *fhl2b*); class shares on AC_fAstCal68 against Quah et al.; labels against Repbase | **yes, before the library is accepted** |
| Does the owner's delegation (2026-11-13) cover the decisions this document reserves to him by name: the build set, the completeness criteria, tool versions and library sources? | before the final library build |
| `write_manifest.py` for merged versions (`RI_v4`, `RI_v5`, the final library): there is no RepeatModeler run; the manifest records the sources, the merge settings and the merge record | before merged versions are synced |
| Push `RI_v1`–`RI_v3` (and the rehearsal versions) to Dropbox as records? | no |
| Where is *tbx2a* in the annotated genomes? NCBI names only *tbx2b*; a gene described as TBX3 on chromosome 10 may be it (needed for the known insertions) | no |
| Does the final masking add `-a` (alignments), for an exact (Kimura) divergence landscape? | before the final masking |
| How to choose between `RI_v9` and `RI_v10`: proposed are agreement with Ensembl's masked bases (both directions), the class at the published insertions (*tbx2a*: hAT-8; *mitfa*: Rex1; *fhl2b*: the AFC-SINE, once that check is built), and, for information only, the distance from Quah et al.'s composition; a tie to the merge (`RI_v10`), which keeps one family per element. With Repbase in both, the curated cichlid Repbase entries cannot grade them (they are members); the same comparison without Repbase (`RI_v8`, `RI_v6`) can still use them | before the results |
| How to choose a curation tool: proposed is to judge by agreement with the curated library, composites, Unknown share and run time, and, instead of a limit on lost coverage, to require that the bases the curated library masks (known TE copies) stay masked within 0.5 points, with any other lost coverage broken down (in coding sequence, in the shuffled control, elsewhere); the student asked to reconsider coverage limits because our masked shares exceed published ones | after the comparison with Quah et al. on fAstCal1.2 |
| The shuffled-chromosome control (allowed by the principle as an exception): which chromosome, and how it is shuffled (preserving dinucleotide composition) | next package |
| The settings of the tool comparison: proposed in `repeats_027` are MCHelper `-r A -a F` with its defaults and TE-Aid off (its plots serve manual inspection only), BUSCO set `actinopterygii_odb10.2020-08-05` (the release of its authors' own test set; `odb12` sets exist), TEtrimmer `--classify_all` with its defaults | after the pilot |
| Whether to give TEtrimmer the curated Repbase entries through its `--curatedlib` (it drops input families more than 95% identical to them) | with the Repbase decision |
| How to correct mislabelled and composite consensus sequences (relabel from curated matches, trim the other element, or correct labels after masking), once the composite check and the comparison against the curated library have measured the problem | before the final build |
| The *mitfa* insertion in MC_MC1m (2,087 bp intron 1 against 678 bp without it): our library masks 1,388 bp of it as LTR/Gypsy and LTR/Pao, the published annotation (Carleton et al. 2020, CENSOR) is mainly Rex1 (LINE); which is right? | no |
| What is MZ_MZ4f_rnd-5_family-9801 (Unknown, 13.4 kb consensus, 12.9 Mb of DT_004f at 3-4% divergence in joined, tandem-like annotations)? Possibly of interest to project 2b, through the coordinator once characterized | no |
| The host-gene tag in the shared `repeats` contract, and whether masked genomes used by other projects should include tagged families (through the coordinator) | before first integration |
| How much coding sequence may the final library mask? (no published threshold; per gene, about 4% of human genes carry transposon sequence in coding regions) | **yes, before the library is accepted** |
| The earlier RepeatMasker annotation of LG2, LG10 and LG13 in seven annotated assemblies (an old RepeatMasker, a library from GT3a alone): another de novo library, so not an accuracy reference; worth a consistency comparison? | no |
| Repbase merge (decided, section 6): which Repbase taxa are retained? The owner and the student decide | before the final masking |
| **Satellites and other tandem repeats:** does project 1 build and annotate them, or project 2b (centromeres, arrays)? (relayed to the coordinator) | before the library is accepted |
| Library policy for new panels: frozen or rebuilt? | before P2 |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) masked too? | no |
| For the coordinator: AS_Sanger (`fAulStu2.1`) and AS_SAMEA123396442 (`fAulStu2.pri`) carry the same specimen name. Same individual? | before AS is used |
| Permanent place for comparison results: a `comparisons/` folder under `RepeatInformation` would need the coordinator to add it to `ga_common.paths` and `sync.py`; until then they stay in scratch | before the pilot is reported |

---

## 8. Status

First conversation done (2026-10-28); decisions in section 6. The 25 P1 genomes are in the data root on mzebra. Container:
TETools 2.01 (`sha256` 8b805d16d279eb57e4948719d2acd8dcd75b31b5ce22fbadcd3a522307fe97ca) under Apptainer 1.5.3 in `ga-repeats`;
RepeatModeler 2.0.9, RepeatMasker 4.2.6 and CD-HIT 4.8.1, confirmed from the tools' own output.

**Versions** (manifests not yet written; nothing synced to Dropbox):

| version | what | families | DT_004f masked |
|---|---|---|---|
| `RI_v1` | RepeatModeler on MZ_MZ4f (seed 1793333684, 36:14:33; BuildDatabase 61 sequences, 972,958,384 bp) | 2,396 | 40.81% |
| `RI_v2` | RepeatModeler on MZ_NMetZ (seed 1793758093, 33:53:39; 75 sequences, 959,220,039 bp) | 2,398 | 40.33% |
| `RI_v3` | RepeatModeler on YH_YH2f (seed 1793758241, 38:42:38; 69 sequences, 976,509,173 bp) | 2,306 | 40.60% |
| `RI_v4` | merge of `RI_v1`–`RI_v3`, `-aS 0.80` (rehearsal) | 3,849 | 43.72% |
| `RI_v5` | merge of `RI_v1`–`RI_v3`, `-aS 0.80 -aL 0.80` (rehearsal; the final build's settings) | 5,288 | 43.91% |

Each RepeatModeler run: five rounds, 16 threads, about 42% database coverage, classified against Dfam 4.0. Masked
percentages include RepeatMasker's own simple repeats and low-complexity sequence.

**Question 1 pilot, results** (`repeats_009`–`repeats_014`; the comparison results stay in scratch, section 7):
- Families: 36% of family groups hold families from both libraries within species, 32% between lineages (cd-hit-est,
  0.80); the 4-point gap holds from 0.80 to 0.95. 95–96% of families have a match in the other library, 57–62% over at
  least 80% of their length.
- Masking: every library masks every genome at 39.9–41.6%; YH_YH2f 41.56% with MZ_MZ4f's library and 41.60% with its own.
  A second library masks 2.1–2.8% of a genome that the first misses, within species and between lineages alike; on
  DT_004f a second library adds 2.41 points, a third 1.16.
- Labels (one label per base): about four-fifths of the sequence two libraries both mask gets the same class (81.8%
  within species, 79.2–80.3% between lineages); LTR differs most between lineages. YH_YH2f's lower LTR total is mostly
  the same sequence labelled differently; the sequence part comes from `RI_v1`, which masks 0.27–0.38 points of LTR the
  others leave unmasked.
- The answer is in section 6 (2026-11-11).

**Merge rehearsals** (`repeats_015`–`repeats_016`): `RI_v5` keeps 99.4% of what its sources mask together (44.15% of
DT_004f); it has 16 conflicting clusters (132 in `RI_v4`), no classified family relabelled (30 in `RI_v4`), and labels that
agree with each source on 86.2–86.9% of shared sequence (82.6–85.1% for `RI_v4`). 1,499 of its 5,288 families are Unknown.

**Accuracy checks on `RI_v5`** (`repeats_017`): it masks 43.49% of AC_fAstCal68, transposons 41.26% (Quah et al.:
36.55% of *A. calliptera*; DNA 15.39% against 12.09%, not yet explained). It masks 13.0–14.0% of the coding sequence of
DT_005m, MZ_NMK2024a and AC_fAstCal68: about 5–6% in "uncharacterized" gene models (mostly by transposon-classified
families), 3.3–3.7% in known multi-copy gene families (false positives), about 4% elsewhere. Masking took about 2.8 h
per genome with one library, 3.9 h with `RI_v4` and 4.5 h with `RI_v5`.

**Repbase label check** (`repeats_018`): 990 of `RI_v5`'s 5,288 families (19%) have a Repbase match covering at least
80% of them; their classes agree for DNA 99.0%, LINE 98.2%, LTR 81.4%, SINE 50% (14 families). 154 Unknown families have
such a match and could take its class. 1,539 families have Repbase transposon support (a transposon covering half or more).

**Gene filter reports on `RI_v5`** (nothing written yet): the first flagged 936 families by evidence B (host proteins not
cleaned); the second, with cleaned host proteins, flagged 487 by evidence B and 132 by evidence A, 98 by both, and would
have removed 491 families. The 98 flagged by both match NLRC3, NLRs, TRIMs, butyrophilin, Fc and olfactory receptors.

**`RI_v7`** (copy-number filter): merged, not masked; the filter is set aside (section 4).

**`RI_v8`** (concatenation test, DT_004f): 44.00% masked against `RI_v5`'s 43.91%; labels 97.9% the same; fragmentation
unchanged (1,371 against 1,363 annotations per Mb); in 980 of `RI_v5`'s 1,150 multi-member clusters two or more members each
hold 10% of the cluster's bases; run time 5:10 against 4:30. **Tagging report on `RI_v5`:** 304 families tagged by the
rule with the gene-family exception, 178 by coverage only; tagged families mask 6.1-6.4% of coding sequence and 2.05% of DT_004f.

**Divergence landscape, DT_004f, `RI_v5`** (approximate: mismatches, not Kimura; nested annotations each counted, so
classes sum to 45.8% against 43.9% masked): DNA 16.0% of the genome (median divergence 5%), LINE 10.6% (4%), LTR 8.1% (4%),
SINE 0.7% (8%), Unknown 8.8% (10%). About half of DNA, LINE and LTR bases lie within 5% of their consensus and 14-20% within
1%; an older wave at 12-26% holds 34% of DNA and 25% of LINE and LTR bases. Unknown is skewed old (42% at 12-26%) and has a
spike at 3-4% (1.54% of the genome, three to five times its neighbouring bins), to be traced to its families. SINEs have few
young copies (8% within 1%).

**Verification on `RI_v5`, first results** (2026-11-17 and -18):
- *Tagging by coverage only:* 178 families tagged (Unknown 100, LTR 28, LINE 27, DNA 16, Satellite 6, Retroposon 1; categories
  OTHER 67, NLR 31, MULTI 28, IMMUNE 15, TRIM 12, PCDH 7, ZNF 6, OR 5, HIST 2, TUB 2, TAAR 2, TASTE 1); tagged families mask
  5.4-5.7% of coding sequence and 1.26% of DT_004f (the earlier rule: 304 families, 2.05%).
- *Lineage-specific loss:* on its own build genome, `RI_v5` misses 0.107% of what `RI_v1` masks on MZ_MZ4f and 0.122% of what
  `RI_v3` masks on YH_YH2f, and masks 3.2-3.4 points more; 86-87% of shared bases get the same class.
- *LG3:* the most-masked chromosome in all eleven masked assemblies (63-70%, against 43.5-45% genome-wide). Across ten of them its
  length varies 64-87 Mb while its unmasked part stays 23-27 Mb: each extra Mb is 0.86 Mb of masked sequence. DT_004f's LG3 is
  27% gap (28.4 Mb of N) and 67.5% masked over non-gap sequence; DT_004f is 2.98% gap overall. Three assemblies list the
  mitochondrion (MT) as a chromosome.
- *Unknown spike at 3-4%:* one family, MZ_MZ4f_rnd-5_family-9801 (13.4 kb consensus), masks 12.9 Mb of DT_004f in 102
  annotations averaging about 127 kb: joined copies of a tandemly repeated consensus.
- *Known insertions:* *mitfa* intron 1 is 678 bp without the insertion (Carleton et al. 2020, Fig. S4C). MZ_NMK2024a,
  AC_fAstCal68, RC_006pf, YH_YH1m and AS_Sanger have a 676-679 bp intron with no repeats (no insertion, as expected for MZ and
  AC); MC_MC1m's is 2,087 bp (the insertion), 1,388 bp of it masked as LTR/Gypsy and LTR/Pao (section 7); DT_005m and CV_CV1m
  lack that transcript. At *rx1*, no Rex-Babar in any of the eight 5 kb windows; the state needs the sequence-anchored check.

**Verification, second results** (2026-11-18; details in `docs/repeats_findings.md`): `RI_v6` written (5,288 families,
178 tagged). The insertion sites, anchored on Carleton et al.'s Table S4, find the *mitfa* insertion in MC_MC1m only (exactly
1,408 bp) and the *tbx2a* insertion in MZ_NMK2024a only; both are masked, but labelled LTR where the literature reports
Rex1 (LINE) and hAT-8 (DNA). Curated-family recovery: 147 of 165 cichlid entries at least present (89%), but our
best-matching family has the curated class for only 93 of them (63%; DNA 58 of 106). The Unknown spike is a satellite
(237 bp period; best Repbase match TZSAT). With gaps excluded, LG3 is the most-masked chromosome in all eleven assemblies.
The RepeatMasker edition of Repbase records no provenance for its entries ("Culled from Repbase").

**Curation tools** (2026-11-19 to 2026-11-23; decisions in section 6): both installed and run on their own test data
in `<scratch>/tool_tests`. TEtrimmer 1.7.2 kept the same 5 of its 9 test families as its stored standard output, with no
class in conflict; it classified 2 the standard leaves Unknown and kept one more family (the standard was made by an older
version, without the evaluation levels 1.7.2 writes, and our run used `--classify_all`). Its first run spent 9 of 13
minutes downloading Pfam (7.4 GB, kept). MCHelper 1.7.0 turned its 21-family toy library into 11 families in 21 minutes
on 8 threads; its repository holds expected results only for its full test inputs, so the toy output cannot be compared
with them (its clustering and subfamily splitting depend on the whole library); a run of the full fruit-fly input is an
optional check.

**Comparison with published data** (2026-11-23; decisions in section 6): Quah et al. 2025 report TEs for *A. calliptera*
fAstCal1.2 (GCF_900246225.1; PacBio CLR with Illumina, contig N50 4.4 Mb, 0.88 Gb): 36.55%, of which DNA 12.09%, LINE
8.32%, LTR 7.38%, Unknown 8.54%, SINEs, helitrons and retroposons together under 1%; their libraries were RepeatModeler
2.0.3 (`-LTRStruct`) libraries combined with Repbase's RepeatMasker libraries, masked with RepeatMasker 4.1.2-p1 (`-e
rmblast -no_is`). Our AC_fAstCal68 (GCA_964374335.1, Wellcome Sanger 2024, contig N50 28.2 Mb) is a different assembly,
so the 4.7-point gap of findings C4 mixes assembly and library. Our libraries do not yet contain Repbase, which theirs did.

**Curation pilots** (2026-11-23 to 2026-11-24; 100 families of `RI_v1` on MZ_MZ4f, 16 threads; a class-stratified
sample, seed 1, always including families 276, 769 and 9801). TEtrimmer 1.7.2: 77 families out in 9.70 h (LTR 38, DNA 18,
LINE 12, Unknown 8, Retroposon 1). MCHelper 1.7.0: 97 out in 2.64 h, every one classified (Wicker class I 61, class II 36).
At this rate the full `RI_v1` (2,396 families) would take TEtrimmer about 232 h on 16 threads and MCHelper about 63 h.
Of the 23 Unknown families going in, TEtrimmer left 8 Unknown and MCHelper none: MCHelper does not remove families it
cannot classify (as its toy run had suggested) but assigns them a class, and it called the satellite family 9801 a
Helitron (`CLASSII/HELITRON`), against the evidence of findings B7 (237 bp tandem period; best Repbase matches TZSAT and
ONSAT, both satellites); TEtrimmer kept 9801 as Unknown. Neither tool relabelled the *tbx2a* composite 769: TEtrimmer split
it into three subfamilies, all LTR/Gypsy; MCHelper kept it as LTR/Gypsy, marked incomplete. Family 276 stayed LTR/Gypsy
in both. MCHelper relabelled some LTR-pipeline families as DNA transposons (family 168 became four `CLASSII/TIR`
subfamilies); whether that is right needs evidence. These are pilot observations, not findings.

**The reference genome** (2026-11-24): fAstCal1.2 fetched from NCBI, 249 sequences, 880,445,564 bp, every file matching
NCBI's md5 sums; NCBI provides no RepeatMasker run for it. Ensembl release 103 holds its soft-masked genome (toplevel and
per chromosome). RepeatMasker 4.1.0's taxonomy gives *A. calliptera*'s lineage in 35 levels, from the species through
Cichlidae, Teleostei, Actinopterygii and Vertebrata to Eukaryota.

**Next** (`repeats_029`): `repbase_lineage.py` (the Repbase subset of the lineage), `add_to_library.py` (`RI_v9` = `RI_v8`
+ subset, `RI_v10` = `RI_v6` + subset), masking of fAstCal1.2 with both, `fetch_ensembl.py` and `ensembl_masking.py`
(masked bases against Ensembl's), `te_composition.py` (classes against Quah et al.'s), `site_annotation.py --with-version`
(which family wins at the *tbx2a* and *mitfa* insertions once Repbase is in the library), and `compare_variants.py`.
Then *fhl2b* on the reference genome. Nothing is ready to integrate.

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
