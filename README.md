# ga-alignment: project 2a (Alignment)

Pairwise alignment of the 24 P1 assemblies to the reference `MZ_NMK2024a` with minimap2 (PAF with the `cs`
tag), and an unclipped pangenome graph of all 25 with Minigraph-Cactus. The authority is the specification,
`alignment_spec.md`, kept in the Claude project; the conventions are the coordinator's `docs/DELIVERY.md` and
`docs/USING_COMMON.md`.

## Where things live on mzebra

| what | where |
|---|---|
| this repository | `~/GenomeAlignments/ga-alignment` (mzebra only pulls; nobody commits here) |
| conda environment | `ga-alignment`: Apptainer 1.5.3 with squashfuse 0.6.1, Python 3.11, minimap2 2.31, `ga_common` at `common-v0.2.1` |
| Cactus image | `~/containers/cactus_v3.3.0.sif` (pull cache and tmp in `~/containers/cache` and `~/containers/tmp`) |
| working space | `~/scratch/ga-alignment` |
| results | stay in scratch until the coordinator adds a data-root folder for them (`ga_common.paths`, `sync.py`) |

## Scripts

Every script is a dry run unless given `--go`, and never overwrites a result.

- `scripts/check_env.py`: records the evidence that the environment works (tool versions from their own
  output, the image's sha256, AVX2, disk, and which planned options and configuration attributes the installed
  versions have), and exports the conda environment. Writes `~/scratch/ga-alignment/env/<UTC timestamp>/`.

Tests: `python -B -m unittest discover -s tests` (stand-in data only; run on the laptop and on mzebra).

## Status

alignment_001: environment only. No genome has been read and no alignment has been run.
