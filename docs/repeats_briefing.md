# Repeats project: briefing for the first conversation

**For:** the student and the assistant of project 1 (Repeats). **Owner:** McGrath lab.
**Version:** 1.2, 2026-10-26. Written by the coordinator; the project decides its own methods.

This briefing is where the first conversation starts. It has two questions that have to be settled
with the owner before any RepeatModeler run, and one paper that has to be gone over first because it
is the closest prior work on this exact system.

1. **How many genomes should the repeat library be built from?** All 25, one per species (8), or fewer?
2. **How confident can we be that the library finds all the transposon families?** What is known to be
   resolved poorly, and how would we measure it?

Everything below separates three kinds of statement: **what a source says** (with where), **my reading**
(the coordinator's inference, marked as such) and **what is not known**. Nothing has been run yet.

---

## Part A. The pangenome paper, to be gone over first

**Quah, F. X., Almeida, M. V., Blumer, M., Yuan, C. U., et al. (2025).** Lake Malawi cichlid pangenome
graph reveals extensive structural variation driven by transposable elements. *Genome Research* 35:1094–1107.
Open access: https://genome.cshlp.org/content/35/5/1094 (doi 10.1101/gr.279674.124).

The student reads the paper's Results section "Enrichment of TEs in SVs" and the Methods section
"Transposable elements" in full before the conversation continues. The summary below is a map, not a
substitute.

### A1. What they did

- Eight long-read assemblies of seven Lake Malawi species (six new, two published: *A. calliptera* fAstCal1.2
  and *M. zebra* UMD2a), joined into a pangenome graph with minigraph, with *A. calliptera* as the backbone.
- The new assemblies are contig-level (contig N50 0.6–8.4 Mb; the ONT assemblies are the most fragmented,
  BUSCO 92.5–93.9%). **Ours are far better: contig N50 22–41 Mb across the 25 panel assemblies.**
- **TE annotation method** (Methods): RepeatModeler v2.0.3 and RepeatMasker v4.1.2-p1 from the Dfam TETools
  container v1.3. De novo families came from `BuildDatabase` and `RepeatModeler` with `-LTRStruct`. The de novo
  libraries were **combined with Repbase-derived RepeatMasker libraries**, and all genomes were annotated with
  RepeatMasker (`-e rmblast -no_is`).

### A2. What they found (numbers as reported)

| finding | value |
|---|---|
| extra sequence in the 7-genome graph, relative to the backbone | 291.47 Mb, **33.1%** |
| share of that extra sequence private to a single assembly | 62.5% (182.16 of 291.47 Mb) |
| structural-variant alleles that are singletons | 60% |
| base pairs that differ structurally between species (pairwise) | **4.73–9.86%** (mean 7.11%) |
| the same measure between two *Rhamphochromis* individuals (same species) | 2.96–3.68% |
| *A. calliptera* genome annotated as TE | **36.55%** |
| TE share in the other species (genome-wide) | 36.32–38.29% |
| TE classes in *A. calliptera* | DNA transposons 12.09%, LINEs 8.32%, LTR 7.38%; SINEs, helitrons and retroposons together under 1%; **unknown 8.54%** |
| structural-variant sequence annotated as TE | **74.65%** (72.26–75.74% across species), 2.04× the genome-wide share |
| presence/absence TE insertion events | 254,015 events in 90,451 variant bubbles |
| most frequent insertion subclasses | LTR/Gypsy 38,216; **LTR/Unknown 28,080**; DNA/TcMar 21,518; DNA/hAT 18,947; LINE/Rex-Babar 18,292; LINE/L2 13,077 |
| genes with a polymorphic TE bubble in the body or 2 kb upstream | 55.9% (29.7% if near the start); no functional enrichment |

Known cichlid TE insertions were recovered: a haplochromine-specific SINE upstream of *fhl2b* (egg-spot
pattern), a nested insertion upstream of *rx1* (opsin expression) and TE insertions in the first intron of *mitfa*.

### A3. What the authors say they cannot conclude

- They **decline to compare TE composition between species** because assembly quality varies.
- PCR validation was hard in repetitive regions (primer design, Sanger reads).

### A4. What the paper does not tell us (the gaps our project has to fill)

- **Which genomes the de novo libraries were built from.** The Methods say "libraries" (plural) and do not
  say whether each genome got its own, or how the libraries were merged.
- **Whether any TE family is missing.** Completeness of the library was not tested. The TE shares are
  reported for the combined library only.
- **Family-level differences between species.** The paper counts insertions (copies at loci), not which
  families are present or absent in each lineage. A new family private to one lineage would be invisible here.
- **Satellites and other tandem repeats.** Not treated separately.

### A5. Discussion points for the student (the coordinator's reading, to be challenged)

1. **Why this matters to us.** Structural variants between these genomes are about three-quarters TE-derived.
   Project 3 (Variant finder) has to call transposon insertions with a class, and project 2b has to tell
   inversion breakpoints and repeat arrays from assembly noise. Both depend on this library.
2. **SNP divergence understates the difference.** Species differ by 0.1–0.25% at SNPs but 4.7–9.9% of
   base pairs structurally. A "very close" panel is close in SNPs, not necessarily in repeats.
3. **Copies differ more than families, probably.** Most extra sequence is private insertions, which a
   shared library can mask if the family is in the library. Whether every family is in the library is
   exactly question 1 and question 2 below. The paper does not answer it.
4. **Unknowns are a visible gap.** 8.54% of the genome is annotated as unknown TE, and LTR/Unknown is the
   second most frequent insertion subclass (28,080 events). That is a known resolution gap, in this system.
5. **A small class can carry large biology.** SINEs are under 1% of the genome, yet one SINE insertion
   underlies the egg-spot pattern. Shares of the genome are not the right yardstick for whether a class matters.
6. **Comparability.** Their tool versions (RepeatModeler 2.0.3, RepeatMasker 4.1.2-p1) are close to those seen on
   mzebra (RepeatModeler 2.0.2a, RepeatMasker 4.1.0). Newer releases exist (section C4), so a different
   result from ours could come from the version, the genomes or the library, and we must record all three.
7. **Their genomes were fragmented; ours are not.** More repeat copies resolve into chromosomes in ours.
   Their caution about species-level TE comparisons may weigh less for us, but that is a hypothesis.

---

## Part B. Question 1: how many genomes for the library?

### B1. Facts about the software

- RepeatModeler **samples** the genome in rounds instead of reading all of it. Round 1 samples about 40 Mb
  (RepeatScout); later rounds (RECON) grow by roughly a factor of 3. A logged 2.0.4 run showed samples of
  about 40, 10, 30, 90 and 270 Mb. If too little long sequence remains it stops early. The sample sizes and
  rounds completed are in the run log (`rmod.log`), which must be saved.
- So **25 genomes would not give 25 times the sampling**. A single genome of about 1 Gb is already only partly sampled.
- **Never concatenate the genomes into one `BuildDatabase`.** This is the coordinator's reasoning, not a
  documented warning: with 0.1–0.25% SNP divergence, every ordinary locus would appear about 25 times and could
  be modeled as a repeat.
- Output per run: `<name>-families.fa` (consensus library), `<name>-families.stk` (seed alignments), `rmod.log`.

### B2. Facts about the panel

- 25 assemblies, 8 species, all seven ecogroups covered (species table in the coordinator repository).
- Contig N50 ranges 22.14–40.94 Mb. Within a species it differs by only a few Mb, so **assembly quality does
  not separate the candidates**. The annotated assembly of a species is usually not the one with the highest N50.
- Annotated assemblies (RefSeq): MZ_NMK2024a, AC_fAstCal68, CV_CV1m, DT_005m, MC_MC1m, RC_006pf, YH_YH1m,
  AS_Sanger. These are one per species.

### B3. The options

| option | runs | notes |
|---|---|---|
| one genome (for example the reference) | 1 | cheapest; may miss families whose copies sit mostly in other lineages |
| one genome per species | 8 | libraries merged and de-duplicated afterwards; covers every species and, with 2 deep-benthic species, all 7 ecogroups |
| one genome per ecogroup | 7 | same idea, fewer runs |
| every genome separately | 25 | not a larger sample per run; most cost, least new information expected |

### B4. How to decide (to be designed with the owner, measured not assumed)

The measurement the student can run: build a library from one genome, run RepeatModeler on a second genome from
a different ecogroup, and compare the two libraries. Report at least: how many families in the second run are
new (not covered by the first library, by a stated similarity rule), how much of each genome those new families
mask, and how the masked fraction of a held-out genome changes when it is masked with each library. If new
families are few and mask little sequence, one genome is enough for now; if not, add genomes until the gain stops.

The decision, and the similarity rule for "new", are the student's to propose and the owner's to agree. Record
them in the specification (section 6) with the numbers behind them.

---

## Part C. Question 2: how confident can we be that all families are found?

**Short answer: no confidence yet, and "all" is not something a run can show.** Completeness has to be measured,
and the measurement has limits. Nothing has been run.

### C1. Where the literature and this paper point to weak spots

| area | what is known | source |
|---|---|---|
| unclassified families | 8.54% of the *A. calliptera* genome is "unknown" TE; LTR/Unknown is the second most frequent insertion subclass | Quah et al. 2025 |
| DNA (TIR) transposons | largest class in *A. calliptera* (12.09%); in a published animal benchmark RepeatModeler2 reached about 39% sensitivity for TIR elements against about 81% for LTR and 76% for non-LTR | Ou et al. 2024 (Genome Biology, response to a commentary) |
| SINEs | under 1% of the genome but biologically important here (egg spots); SINE annotation is described as marginal for RepeatModeler2 and EDTA | Quah et al. 2025; Ou et al. 2024 |
| LTR elements | need `-LTRStruct`; the LTR and the internal part are often separate entries, the LTR part often with an unknown superfamily | RepeatModeler2 paper; Goubert et al. 2022 |
| tandem repeats and satellites | RepeatModeler targets interspersed repeats; satellites need separate tools; centromeric satellites matter to project 2b | RepeatModeler documentation; Blumer et al. 2025 |
| consensus quality | automatic consensus sequences are drafts; a reliable library needs manual curation | Goubert et al. 2022 |

Reading (the coordinator's): DNA transposons are the largest class in these genomes and also the class the
benchmark found hardest, so a sensitivity gap would matter most there. That is a hypothesis for the student to test.

### C2. Ways to build confidence (each one is a measurement)

1. **Held-out genome.** Mask a genome that was not used to build the library; how much sequence is masked
   compared with masking it with a library built from itself?
2. **Residual repeats.** After masking, does RepeatModeler (or a high-copy k-mer count) still find repeats in the
   masked genome?
3. **Class composition against the literature.** Compare class shares with the published *A. calliptera* values
   (DNA 12.09%, LINE 8.32%, LTR 7.38%, unknown 8.54%); large differences need an explanation (version, library, genome).
4. **The known insertions.** Are the *fhl2b* SINE, *rx1* and *mitfa* TE insertions recovered and classified?
5. **Unknown fraction.** Report it for the library (families) and for the masked genome (base pairs), and decide what to try for the
   unknowns (curation, a second classifier, similarity to Dfam).
6. **Divergence landscape.** Do young (low-divergence) and old families both look sensible per class?

### C3. What "good enough" means

That is the owner's decision, per class, and it is recorded in the specification (section 5, checks). The project
must report, for every run, what was measured and what could not be measured.

### C4. Tool versions

- mzebra's conda environment `GenomeSize` holds RepeatMasker 4.1.0 and RepeatModeler 2.0.2a (seen 2026-10-26).
- The RepeatMasker site (news of 11 June 2026) announces a RepeatModeler maintenance release that works with
  RepeatMasker 4.2.4 and the new Dfam 4.0 FamDB format; the Dfam TETools container follows it. mzebra has no
  docker or singularity, so a conda or source install would be needed.
- The student records exact versions (and the Dfam version, if Dfam is used) for every run, and the owner
  decides whether to move to current releases **before** the library is built.

### C5. Whether to combine with Dfam or Repbase

Quah et al. combined their de novo library with Repbase-derived libraries. Repbase has licence terms; Dfam ships
with RepeatMasker. Whether to combine, and with what, changes family names and class assignments, so it is tied to
the open question on keeping existing family names (specification section 7).

---

## Part D. Sources, in reading order

1. **Quah et al. 2025**, Genome Research 35:1094. https://genome.cshlp.org/content/35/5/1094 . *Read first, in full.*
2. **RepeatModeler README**: https://github.com/Dfam-consortium/RepeatModeler and https://www.repeatmasker.org/RepeatModeler .
   Commands, `-LTRStruct`, threads, output files.
3. **Goubert's workshop tutorial** (RepeatModeler2 then RepeatMasker, step by step, with further reading):
   https://www.repeatmasker.org/~cgoubert/UofA_Cyverse_Workshops/2024_TE_Workshop/TE_workshop_V1.html
4. **Flynn et al. 2020**, RepeatModeler2, PNAS. https://pmc.ncbi.nlm.nih.gov/articles/PMC7196820
5. **Goubert et al. 2022**, a beginner's guide to manual curation of transposable elements, Mobile DNA.
   https://doi.org/10.1186/s13100-021-00259-7
6. **Ou et al. 2024**, response to a commentary on TE annotation benchmarks, Genome Biology.
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10762968/ (the sensitivity table by TE class).
7. **Dfam classification** (how families are named): https://www.dfam.org/classification
8. **Carleton et al. 2020**, Movement of transposable elements contributes to cichlid diversity, Molecular Ecology
   (doi 10.1111/mec.15685). Cited by Quah et al.; not read by the coordinator.

### Uploading the sources (the student does this, before the first conversation)

The student saves each source below **as a PDF** and uploads it to the Repeats project, so the assistant has the same
fixed text at the start of every session. Items 1–5 are needed before the first conversation continues; items 6–8 may follow.

| # | source | how to get the PDF | suggested file name |
|---|---|---|---|
| 1 | Quah et al. 2025 | on the journal page, "Download PDF" | `Quah2025_GenomeRes_pangenome.pdf` |
| 2 | RepeatModeler README | browser: Print, then Save as PDF | `RepeatModeler_README.pdf` |
| 3 | Goubert workshop tutorial | browser: Print, then Save as PDF | `Goubert_RM2_workshop.pdf` |
| 4 | Flynn et al. 2020 | the PDF link on the PubMed Central page | `Flynn2020_RepeatModeler2.pdf` |
| 5 | Goubert et al. 2022 | on the journal page, "Download PDF" | `Goubert2022_manual_curation.pdf` |
| 6 | Ou et al. 2024 | the PDF link on the PubMed Central page | `Ou2024_GenomeBiol_TE_benchmark.pdf` |
| 7 | Dfam classification page | browser: Print, then Save as PDF | `Dfam_classification.pdf` |
| 8 | Carleton et al. 2020 | the journal page; use the Georgia Tech library if it is behind a paywall | `Carleton2020_MolEcol_TE_cichlids.pdf` |

How:
- **Papers:** use the publisher's or PubMed Central's own PDF, not a print of the web page.
- **Web pages:** in the browser, Print, then choose Save as PDF. Open the result and check that the text can be selected
  (it should not be a picture of the page) and that the whole page is there.
- **Upload:** in the Repeats project, add the PDFs to the project's files.
- **Check:** ask the assistant, "List the documents in this project and give the title of each." Re-upload anything that
  is missing or unreadable.
- **Version note:** run `RepeatModeler -version` and `RepeatMasker -v` once the tools are available, and tell the assistant
  whether the uploaded README matches the installed version. If it does not, say so.

---

## Part E. Questions this raises for other projects (relayed through the coordinator)

- **Satellites and tandem repeats:** does project 1 build and annotate them, or does project 2b, which already
  owns centromeres and arrays? (Open in the coordinator specification.)
- **Repeat annotation as an input to alignment** (project 2a) and to variant classification (project 3): which
  class and family levels do they need?
