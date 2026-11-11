# Question 1 pilot scripts (repeats_009; cluster reading fixed in repeats_010; overlap check and figures in repeats_011; gain curve in repeats_012; label cross-tabulation in repeats_013)

Run on mzebra, from the repository root, in the conda environment `ga-repeats`
(`conda activate ga-repeats`). Every script prints what it would do; it acts only with
`--go`, never overwrites anything, takes `--threads` and `--max-jobs` (environment
`GA_MAX_THREADS`, `GA_MAX_JOBS`), and finds data through `ga_common`. Defaults: image
`~/containers/dfam-tetools-2.01.sif`, scratch `~/scratch/ga-repeats`, 16 threads.

| script | answers | writes |
|---|---|---|
| `compare_libraries.py --a RI_v1 --b RI_v2` | which families are shared or found in only one library (cd-hit-est at 0.80, 0.85, 0.90, 0.95 identity, coverage of the shorter sequence 0.80); how similar the best matches are (RepeatMasker, each library's consensus sequences against the other) | `<scratch>/comparisons/RI_v1_vs_RI_v2/` |
| `mask_genomes.py --version RI_v1 --assemblies MZ_MZ4f DT_004f` | gate 4: masks genomes with one version's library only; checks the BED against the genome | `.out` and `.bed` in the version's `masked/<assembly_id>/`; work and run records in `<scratch>/<version>/masking/` |
| `summarize_masking.py --comparison <scratch>/comparisons/RI_v1_vs_RI_v2 --threshold 0.80` | what the differences mean in the genome: base pairs per class, and base pairs masked by families found in only one library | `masking_summary_c0.80.tsv` in the comparison folder |
| `compare_masking.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3` | the second check: on one genome masked by several libraries, the bases masked by both and by only one, whole genome and per class | `<scratch>/comparisons/masking_overlap/<assembly>.tsv` |
| `masking_gain.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3 --plot` | does adding libraries keep adding masked sequence? For every order of adding them: what each adds beyond the earlier ones (per class) and the running total; briefing B4's "until the gain stops" | `<scratch>/comparisons/masking_gain/<assembly>.tsv` and, with `--plot`, `<assembly>_gain.png` and `.pdf` |
| `label_agreement.py --assembly YH_YH2f --versions RI_v1 RI_v3 --plot` | do two libraries give the same bases the same class? Class matrix in base pairs (of the bases one library labels with a class, what the other calls them, or unmasked), agreement on bases masked by both, and the family pairs carrying most disagreement | `<scratch>/comparisons/label_agreement/<assembly>_<A>_vs_<B>.tsv`, `..._families.tsv` and, with `--plot`, `..._matrix.png` and `.pdf` |
| `plot_pilot.py --comparisons <scratch>/comparisons/RI_v1_vs_RI_v2 <scratch>/comparisons/RI_v1_vs_RI_v3` | figures: shared clusters by threshold and by class, best matches, class composition per library and genome (after Quah et al. Fig. 7A), masking overlap; needs matplotlib | `<scratch>/comparisons/plots/` (PNG and PDF) |
| `write_manifest.py --version RI_v1` | the version's record: build, checks, tools, library checksums, masked set, commands | the version's `MANIFEST.yaml`, once, after its masking |

Facts the scripts rely on, read from the installed tools (2026-11): RepeatModeler 2.0.9,
RepeatMasker 4.2.6, CD-HIT 4.8.1 at `/opt/cd-hit/cd-hit-est` in the image; RepeatMasker's
`-pa` runs jobs of 4 cores each with RMBlast; library names are `name#class/subclass`;
the `.out` columns are those of RepeatMasker's help, section 3.1. cd-hit-est word lengths
per threshold follow the CD-HIT user's guide.

Known limits: the comparison results live in scratch until their permanent place is
decided (spec section 7); the BED format is a pilot format (BED6 plus class, divergence,
score and RepeatMasker id) until the `repeats` contract fields are agreed; `library/` and
`masked/` are located inside the version folder that `ga_common` gives for
`repeat_manifest`, following spec section 3.
