# Repeats specification

**Project:** 1, Repeats
**Repository:** `ga-repeats` *(name assumed)*
**Document version:** 0.1, 2026-10-26 (first draft, from coordinator specification v0.4)
**Owner:** McGrath lab
**Status:** draft for discussion with the owner; no code written

This is the single specification document of project 1. It is always replaced whole.
It records what the coordinator specification already fixes, and leaves everything
else open. The methods are not designed here: the project designs them with the owner.

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
| assemblies: one genome FASTA per panel assembly | panel P1 (25 assemblies, by accession version) | `<data root>/Genomes/<assembly_id>/<accession_version>/genome.fna.gz` (BGZF) | **not yet downloaded** (coordinator task) |
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

**Not designed yet.** Fixed by the coordinator specification:
- one shared library for all assemblies (decision 2026-10);
- whole genomes, not selected chromosomes (decision 2026-10-26).

Raised, not decided: how the library is built (which genomes; RepeatModeler with
`-LTRStruct`?); whether existing family names are kept now that the family table is
recreated; the library policy when a new panel arrives (frozen library and mask new genomes
only, or rebuild and re-mask all).

The plan (outputs, formats, checks) is agreed with the owner and recorded here **before any
code is written**.

---

## 5. Checks (proposed; to be agreed)

- every output validates against the `repeats` contract version (once defined);
- every assembly in the masked set has output, and sequence names and lengths in the output
  match the assembly's `provenance.json`;
- `MANIFEST.yaml` is present, complete, and names the panel version;
- the version folder is new: nothing is overwritten (the sync script refuses to);
- commands and tool versions are recorded so a run can be repeated.

---

## 6. Decisions

| date | decision | reason |
|---|---|---|
| 2026-10 | One shared library for all assemblies | coordinator decision |
| 2026-10-26 | The folder is named `RepeatInformation`, with `versions/RI_vN` | repeat modelling may be run on different subsets of genomes, and every run must stay traceable |
| 2026-10-26 | `RepeatModeler_TE.xlsx` is recreated by this project | the pilot's library location is unknown; the old library is not needed |
| 2026-10-26 | Where the tools run is this project's choice | the owner wants flexibility |
| 2026-10-26 | Whole genomes | coordinator decision |

---

## 7. Open questions

| question | blocking |
|---|---|
| **The P1 genomes are not downloaded yet** (coordinator: first `fetch.py --go` on mzebra pending) | **yes, for any run** |
| Contract for the shared library, and exact fields of the `repeats` contract (class and family naming, divergence, BED columns) | before first integration |
| Keep the existing family names, or take the names the recreated library gives? | no |
| Which genomes build the library; RepeatModeler with `-LTRStruct`? | no |
| Library policy for new panels: frozen or rebuilt? | before P2 |
| Where the tools run, and how long jobs are kept alive on the shared server (no scheduler; `tmux`, `screen` and `nohup` not checked) | no |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) masked too? | no |
| Whether masked genomes are stored (size) or regenerated from the `.out` | no |

---

## 8. Status

Not started. Specification drafted; waiting for the P1 download and for the owner's plan.
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
- Python 3.9 is the minimum for code that uses the common library.
