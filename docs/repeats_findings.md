# Repeats (project 1): established findings

**Version 4, 2026-11-23** (C4: Quah et al.'s figures are for another assembly). Version 3, 2026-11-19: B5 checked independently; B12–B14 added; C2 resolved as B12; C5 added. A record of what our data show, kept for framing a paper. It is updated whenever a result
becomes solid, and it is not the specification: decisions and open questions stay in `repeats_spec.md`.

**What counts as established:** measured on our data, with a stated method, from recorded inputs that the repository's
scripts can reproduce. Numbers are given as measured, with their caveats. Anything suspected but not yet tested goes in
the last section, not here.

**Common setting** unless stated otherwise:
- **Libraries:** RepeatModeler 2.0.9 with `-LTRStruct`, one per build assembly: `RI_v1` (MZ_MZ4f), `RI_v2` (MZ_NMetZ),
  `RI_v3` (YH_YH2f).
- **Merged library:** `RI_v5`, the three merged. `RI_v6` is `RI_v5` with gene-derived families tagged.
- **Masking:** RepeatMasker 4.2.6, `-engine rmblast -no_is`, from the Dfam TETools 2.01 image. "One label per base"
  means the highest-scoring annotation covering each base.
- **Genomes:** panel P1. DT_004f was never used to build a library and serves as the held-out genome.

---

## A. Building the library (question 1)

**A1. Libraries built from different assemblies mask about the same share of a genome, but not the same sequence.**
- **Share:** each single library masks every genome tested at 39.9–41.6%. On held-out DT_004f, the three give
  40.3–40.8%.
- **Different sequence:** each masks 2.1–2.8% of a genome that another library leaves unmasked.
- **Diminishing returns:** on DT_004f, one library masks 40.6% on average; adding a second gains 2.4 points, a third
  1.2 more, to 44.15% for all three together.

**A2. RepeatModeler models the same repeats differently from one run to the next, even within a species.**
- **Shared families:** two libraries of *M. zebra* (MZ_MZ4f, MZ_NMetZ) share only 36.2% of family groups at 80% identity;
  libraries from different lineages share 32.3%.
- **Robust to the threshold:** the within- versus between-lineage gap stays about 4 points from 80% to 95% identity.
- **Yet nearly everything has a counterpart:** 95–96% of families have a match in another library, 57–62% over at least
  80% of their length, at a median divergence of 2.2–2.9%.

**A3. Class labels disagree between libraries, most for LTR elements.**
- **Agreement:** 81.8% of bases masked by two libraries get the same class within a species, 79.2–80.3% between
  lineages.
- **Kinds of disagreement:** 54–59% are two different transposon classes, 33–37% Unknown against a class.
- **Concentration:** the top 50 family pairs carry 30–39% of the disagreement.

**A4. Merging at 80% identity over 80% of both families keeps families and labels intact.** Setting: cd-hit-est
`-c 0.80 -aS 0.80 -aL 0.80`, local alignment.
- **Conflicts:** 16 label conflicts among 5,288 clusters (0.3%), and no classified family overruled.
- **Coverage:** `RI_v5` masks 43.91% of DT_004f, against 44.15% for its three sources together, a loss of 0.25 points.
- **Faithful representatives:** members cover their cluster's representative over a median 99.9% (at least 79.5%).
- **Against one-sided coverage** (`-aS 0.80` only, `RI_v4`): 132 conflicts, 30 classified families relabelled, and a
  loss of 0.44 points.

**A5. Merging beats concatenating the libraries.** `RI_v8`, the three libraries concatenated, against `RI_v5`, on DT_004f:
- **Coverage:** concatenation gains only 0.09 points (44.00% against 43.91%).
- **Unchanged:** labels (97.9% the same on shared bases) and fragmentation (1,371 against 1,363 interspersed annotations
  per Mb).
- **The cost:** in 980 of the 1,150 multi-member clusters, two or more members each hold at least 10% of the cluster's
  bases. About 164 Mb of annotation would divide one family's copies among several names.
- **Time:** masking takes 5:10 h against 4:30 h.

**A6. Merging costs a build genome almost nothing of its own masking.** On each build genome, the merged library against
that genome's own library:

| build genome | own library | merged | own library masks, merged doesn't | same class on shared bases |
|---|---|---|---|---|
| MZ_MZ4f | 40.43% | 43.54% | 0.107% | 86.1% |
| YH_YH2f | 41.60% | 44.92% | 0.122% | 87.1% |

**A7. A full-length copy-number filter (panEDTA's rule) does not suit these genomes.** The rule keeps a family only with
at least 3 full-length copies (80-80-80) in its own genome.
- **Families dropped:** 1,290 of 7,100 (18%). That includes 40% of families of 8 kb or more, and 74% of the dropped
  families have one or two full-length copies.
- **What was lost:** large LINE, LTR and DNA families that rarely occur full-length, such as L1 at 3.0 Mb and
  hAT-Charlie at 3.5 Mb on DT_004f.
- **Cost:** in their own libraries' masking, the dropped families account for 4.3–5.2% of DT_004f (an upper bound).

**A8. Masking time grows with library size.** About 2.75–3.4 h per genome with one library (about 2,400 families);
4.5 h with 5,288 families; 5.2 h with 7,100. All at 16 threads per genome.

## B. Accuracy and completeness (question 2)

**B1. The de novo library masks host genes.**
- **Coding sequence masked by `RI_v5`:** 12.9–14.0% in three annotated genomes (DT_005m, MZ_NMK2024a, AC_fAstCal68).
- **In known multi-copy gene families:** 3.3–3.7% of coding sequence (NLRs, TRIMs, immune, olfactory and other
  receptors, zinc fingers, histones, tubulins).

**B2. 178 families are built from host genes, and they're a small share of the genome.** Of `RI_v5`'s 5,288 families,
178 are judged gene-derived and tagged in `RI_v6`. The evidence: gene overlap on three genomes, and a host-protein match
covering at least 30% of the family that beats any transposon protein. None has Repbase transposon support.
- **Their footprint:** they mask 5.4–5.7% of coding sequence, but only 1.26% of DT_004f.
- **Counting small gene-family matches too** (under 30% of the family) would add 126 families, mostly LTR and LINE.
  Those mask 0.79% of the genome but little coding sequence, so they behave like transposons carrying a gene-like
  stretch.

**B3. Labels agree with curated Repbase except for LTR and SINE.** Of the 990 families with a Repbase match covering at
least 80% of them, 79.2% share Repbase's class:

| class | same class as Repbase |
|---|---|
| DNA | 99.0% |
| LINE | 98.2% |
| LTR | 81.4% |
| SINE | 50% (7 of 14) |

Another 154 Unknown families have such a Repbase match.

**B4. The library recovers 89% of curated cichlid families.** Measured against Repbase's 165 cichlid entries, in Flynn et
al. 2020's categories:

| class | perfect | good | present | not found | total |
|---|---|---|---|---|---|
| DNA | 65 | 15 | 26 | 9 | 115 |
| LINE | 8 | 2 | 8 | 6 | 24 |
| LTR | 4 | 0 | 15 | 1 | 20 |
| SINE | 1 | 0 | 1 | 1 | 3 |
| Satellite | 0 | 0 | 2 | 1 | 3 |
| **all** | **78** | **17** | **52** | **18** | **165** |

- **147 of 165 (89%) are recovered at least as present,** 95 (58%) as perfect or good.
- **DNA transposons are recovered 92% of the time.**
- **The missing include 6 LINE/L2 entries.**

**B5. Published insertions are detected, and an independent annotation confirms two misclassifications.** The probes
are Carleton et al. 2020's Table S4 sequences; every one aligns end to end at 96–100% identity in all 11 genomes tested.
Each inserted sequence was also annotated with all of Repbase, independently of our library.

| locus | which genomes carry the insertion | Repbase annotation of the inserted sequence | our annotation |
|---|---|---|---|
| *tbx2a* upstream | MZ_NMK2024a only; MZ_MZ4f, also *M. zebra*, lacks it | **hAT-8_AFC (DNA/hAT-Ac), 97% identity over 1,073 bp** (positions 114–1186 of the 1,332 bp probe), as published | the same stretch (114–1201) as **LTR/Gypsy**, by `MZ_MZ4f_ltr-1_family-769`, a composite LTR family (B14) |
| *mitfa* intron 1 | MC_MC1m only, exactly 1,408 bp between the flanks; none in the other 10, including the negative controls (two *M. zebra*, *A. calliptera*) | **Rex1 (LINE/Rex-Babar) over 629–1177**: Rex1-1_PM (lamprey) at 81–84%, Rex1-6_OL (medaka) at 97% over 59 bp, as published; nothing elsewhere in the insertion | 629–1177 as **LTR/Pao** (`MZ_NMetZ_ltr-1_family-392`); 6–617 and 1179–1408 as LTR/Gypsy (`MZ_MZ4f_ltr-1_family-276`), where Repbase finds nothing |

- **Detection is good:** 1,388 of the 1,408 *mitfa* bases and the whole hAT-8 stretch at *tbx2a* are masked.
- **The library already has the right family for hAT-8_AFC:** a perfect match, `YH_YH2f_rnd-5_family-4640`, though
  labelled Unknown. At *tbx2a*, RepeatMasker preferred the composite LTR family.
- **At *rx1*,** 10 of the 11 genomes carry the deletion, like *M. zebra*; AC_fAstCal68 is B12.

**B6. LG3 is the most-masked chromosome, and its length differences are repeat.** In all 11 assemblies, LG3 is masked
63.5–70.5% over non-gap sequence, against 43.5–45.3% genome-wide. Its non-gap length varies from 64.0 to 86.2 Mb, but its
unmasked part stays at 22.9–26.6 Mb.

**B7. A satellite family masks 13–20 Mb per genome, but is labelled Unknown.** The family is `MZ_MZ4f_rnd-5_family-9801`.
- **Its consensus:** 13.4 kb, internally repeated with a 237 bp period.
- **Repbase:** its best match is TZSAT, a satellite, covering 75% of the family.
- **Its arrays:** 99% of its bases lie in arrays of 10 kb or more, the largest 0.4–1.4 Mb. It's present on most
  chromosomes, often near one end.
- **Per genome:** 12.6–20.2 Mb in the four genomes profiled.
- **The landscape spike:** it produces 12.9 of the 15.4 Mb of Unknown at 3–4% divergence in DT_004f.
- **Genome-wide:** of the 40.6 Mb of DT_004f that curated satellite entries mask, `RI_v5` labels 52.3% Unknown, 20.1 Mb of
  it through this family (B13).

**B8. The divergence landscape shows two waves of activity.** Approximate (mismatch-based) landscape of DT_004f:
- **Recent wave:** about half of DNA, LINE and LTR bases lie within 5% of their consensus, and 14–20% within 1%.
- **Older wave:** at 12–26% divergence, holding 34% of DNA bases and about 25% of LINE and LTR bases.
- **Unknown is older** (median divergence 10%), and **SINEs have few young copies** (8% within 1%).

**B9. `RI_v5` masks 43.5–45.3% of each genome's non-gap sequence** across 11 assemblies.

**B11. The family that best matches a curated cichlid family often carries the wrong class.** Of the 147 curated
cichlid families recovered (B4), the family of ours covering most of each one has the same class in 93 (63%).

| curated class | same class / recovered |
|---|---|
| DNA | 58 / 106 |
| LINE | 17 / 18 |
| LTR | 17 / 19 |
| SINE | 1 / 2 |
| Satellite | 0 / 2 |

- **Where the 54 disagreements come from:** 42 are RepeatScout/RECON families, 12 LTR-pipeline families.
- **Not all are errors:** 19 are labelled Unknown, which means unclassified rather than wrong.
- **One example:** 11 LTR-pipeline families labelled LTR best-match curated DNA transposons.
- **How this fits B3:** B3 says that families *labelled* DNA are almost always DNA. B11 says that many curated DNA
  families are best represented by families labelled otherwise. The first is precision, the second recall.
- **The satellite again:** the *O. niloticus* satellite ONSAT is matched best by `MZ_MZ4f_rnd-5_family-9801`, the
  satellite of B7.

**B12. AC_fAstCal68 carries an insertion at *rx1*, containing Rex1-5_AFC.**
- **Size:** 840 bp lie between the flanks: the 413 bp ancestral stretch plus about 427 bp.
- **Repbase annotation:** Rex1-5_AFC (LINE/Rex-Babar) at 265–398 (90% identity) and 428–659 (98%). That is the element
  Carleton et al. 2020 report in the long-palette species' *rx1* insertion.
- **Our annotation:** LINE, with no subclass (`MZ_NMetZ_rnd-5_family-4263`), over 243–659, with an (AC)n between.
- **The other 10 genomes carry the deletion.**
- **Caveat:** Carleton et al. report *A. calliptera* as ancestral, so the state differs between individuals or
  populations of *A. calliptera*. This rests on one assembly.

**B13. Against the curated library, genome-wide, detection is near complete but classes often differ.** DT_004f was
masked with the 165 curated cichlid entries alone, and compared base by base with `RI_v5`.
- **The curated entries alone mask 23.17% of DT_004f:** DNA 10.99%, LINE 6.55%, satellites 4.08%, LTR 1.33%, SINE 0.22%.
- **`RI_v5` masks 98.9–99.5% of those bases in every class.**

| curated class | given the same class by `RI_v5` | given instead |
|---|---|---|
| DNA | 71.5% | Unknown 12.8%, LTR 9.5%, LINE 4.7% |
| LINE | 79.2% | LTR 9.0%, Other 4.3%, Unknown 3.7% |
| LTR | 81.5% | Unknown 7.3%, LINE 7.0%, DNA 4.1% |
| SINE | 50.6% | LINE 15.1%, LTR 12.6%, DNA 11.5% |
| satellites | 39.3% | Unknown 52.3% |

- **The largest single source of disagreement is the satellite family of B7:** Unknown on 20.1 Mb of curated satellite.
- **LTR-pipeline families give 22% of the curated bases that receive another class.**

**B14. 65 families are composite: their consensus contains an element of another class.** A read-only check of every
family against the curated cichlid entries, at 90% identity or more over at least 300 bp, found 243 families carrying
another class's entry:

| status | families | footprint on DT_004f | what it means |
|---|---|---|---|
| **composite** | **65** | 14.67 Mb (1.47%) | a near-complete element of another class inside |
| relabel | 13 | 2.06 Mb | the family *is* that element; one of the 13 (874 kb) is only a naming difference, PLE against LINE/Penelope |
| Unknown with class evidence | 78 | 11.08 Mb | 67 of them DNA |
| partial | 87 | 16.83 Mb | listed only |

- **The composites:** 35 come from the LTR pipeline and 30 from RepeatScout/RECON. Most are LTR or LINE families
  containing DNA transposons.
- **Their consensus sequences are implausibly long:** among the largest, a "SINE" of 11.4 kb and a "LINE/I" of 22.7 kb.
- **Several are among the largest sources of disagreement in B13.**
- **The positive controls:**
  - The *tbx2a* family (769) is flagged for the hAT element its insertion carries, a genuine pass.
  - The *mitfa* family 392 is flagged too, but for another hAT-Ac element. The Rex1 at that site matches Repbase only
    distantly (81–84%), below the check's threshold.
  - Family 276 carries no curated element, and none was found in its part of the insertion either.

**B10. Gaps in assemblies matter for per-chromosome comparisons.** DT_004f is 2.98% gap, including 28.4 Mb (27%) of its
LG3; the other assemblies are 0–1.8% gap. Three assemblies list the mitochondrion as a chromosome.

## C. Not yet established

Each of these has a test planned or open, and none should be cited as a finding.

- **C1. Composite consensus sequences may explain the misclassifications in B5 and B11.** Some consensus sequences may
  contain other elements: an LTR element carrying a Rex1 or hAT copy, or a LINE carrying a DNA transposon. They would
  outscore the dedicated families on those elements' copies and pass on their own label. B11 shows this would not be
  limited to the LTR pipeline. It fits the RepeatModeler2 paper's note that nested (mosaic) structures are not always
  resolved; RepeatModeler 2.0.9 offers no option capping the LTR internal length. The planned test is a composite-family
  check over all families. B14 shows composites exist (65 families, 1.47% of DT_004f), and several are among the largest
  sources of disagreement. How much of the disagreement in B13 they explain is still to be measured.
- **C2.** Resolved: now B12.
- **C3. Family 9801 may be SATA,** the centromeric satellite described in tilapias, whose monomer is 237 bp as far as I
  know. Its arrays' positions fit a centromeric satellite, but position alone doesn't prove it.
- **C5. What the rest of MC's *mitfa* insertion is remains open.** About 840 bp lie outside its Rex1 stretch. Our
  library labels them LTR/Gypsy (family 276), and Repbase finds nothing there.
- **C4. The 4.7-point gap to Quah et al. 2025 is unexplained, and partly an assembly difference.** On AC_fAstCal68 we
  annotate 41.26% as transposons; Quah et al. report 36.55% for *A. calliptera*, but for fAstCal1.2 (GCF_900246225.1;
  PacBio CLR, contig N50 4.4 Mb, 0.88 Gb), not for AC_fAstCal68 (GCA_964374335.1, contig N50 28.2 Mb). Their libraries
  also included Repbase; ours do not yet. The planned test: our libraries on fAstCal1.2 itself, class by class against
  their figures (DNA 12.09%, LINE 8.32%, LTR 7.38%, Unknown 8.54%).
