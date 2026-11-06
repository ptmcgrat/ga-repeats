# Question 1 pilot scripts (repeats_009; cluster reading fixed in repeats_010)

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
