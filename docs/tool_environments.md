# Curation tool environments (MCHelper, TEtrimmer)

Built on mzebra 2026-11-20 to 2026-11-23 (specification section 6). The project's own scripts run in `ga-repeats`,
which is unchanged; they call each tool in its own environment through `conda run -n <environment>`.
Tool sources, databases and test runs are under `~/scratch/ga-repeats/tools` and `~/scratch/ga-repeats/tool_tests`.

| environment | tool | versions inside |
|---|---|---|
| `repeat_family_merging_mchelper` | MCHelper 1.7.0, repository commit `cd808850183db5718222549c84e84378dab8e610` | its authors' exact `MCHelper.yml`: Python 3.8.16, BLAST 2.10.1, HMMER 3.3.2, CD-HIT 4.8.1, MAFFT 7.505, CIAlign 1.0.18 |
| `repeat_family_merging_tetrimmer` | TEtrimmer 1.7.2 (bioconda) | Python 3.10.14, RepeatModeler 2.0.8, RepeatMasker 4.2.3, RMBlast 2.14.1 (`blastn` 2.14.1), Perl 5.32.1, `pfam_scan.pl` |

The exported specifications (`conda env export`, 2026-11-23) are in this repository:
`envs/repeat_family_merging_mchelper.yml` and `envs/repeat_family_merging_tetrimmer.yml`. To rebuild an environment
exactly: `conda env create -f envs/<file>`. The TEtrimmer environment also holds TEsorter 2.0.0 (a dependency of
TEtrimmer), a classifier by transposon protein domains, should a classification step need it.

## Why two environments

1. **Not in `ga-repeats`**: TEtrimmer 1.7.2 requires Python below 3.11; `ga-repeats` has Python 3.11.
2. **Not in one shared environment**: TEtrimmer needs `pfam_scan`, which through BioPerl needs an old samtools binding
   built against zlib below 1.3, while RepeatMasker 4.2.4 and RepeatModeler 2.0.9 need RMBlast 2.17 or later, which
   needs zlib 1.3.2 or later; MCHelper's dependencies (OpenCV, poppler, R) also want the newer zlib. Without version
   constraints the solver fell back to RepeatModeler 1.0.8 and RepeatMasker 4.0.6 with Perl 5.22, whose Perl does not
   run on mzebra (`libnsl.so.1` missing).
3. **MCHelper**: its `MCHelper.yml` pins every package of its authors' machine on Python 3.8, so it gets that
   environment exactly.

## How they were built

```bash
mkdir -p ~/scratch/ga-repeats/tools && cd ~/scratch/ga-repeats/tools
git clone https://github.com/gonzalezlab/MCHelper.git          # needs git-lfs for db/allDatabases.clustered_rename.fa.zip
git clone https://github.com/qjiangzhao/TEtrimmer.git          # used for its tests/ only

# TEtrimmer: the packaged release, with RepeatModeler2-era versions required
conda create -y -n repeat_family_merging_tetrimmer -c conda-forge -c bioconda \
  python=3.10 tetrimmer=1.7.2 "repeatmodeler>=2.0.4" "repeatmasker>=4.1.4"

# MCHelper: its authors' environment
conda env create -n repeat_family_merging_mchelper -f MCHelper/MCHelper.yml

# MCHelper's repository and databases (its README); git-lfs from conda-forge if the clone's checkout fails
cd MCHelper
git lfs install --local && git lfs pull && git reset --hard HEAD
cd db
unzip -o '*.zip'
makeblastdb -in allDatabases.clustered_rename.fa -dbtype nucl          # 51,127 sequences
wget https://urgi.versailles.inrae.fr/download/repet/profiles/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz
tar xvf ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm.tar.gz
mv ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm Pfam35.0.hmm             # 19,903 profiles (grep -c '^NAME')
```

TEtrimmer downloads Pfam itself into the folder given with `-pd`; the first run put it in
`~/scratch/ga-repeats/tools/tetrimmer_pfam` (7.4 GB, 9 minutes).

MCHelper's BUSCO genes (proposed: `actinopterygii_odb10.2020-08-05`, the release of its authors' own test set):

```bash
mkdir -p ~/scratch/ga-repeats/tools/busco && cd ~/scratch/ga-repeats/tools/busco
wget https://busco-data.ezlab.org/v4/data/lineages/actinopterygii_odb10.2020-08-05.tar.gz
tar xzf actinopterygii_odb10.2020-08-05.tar.gz
cat actinopterygii_odb10/hmms/*.hmm > actinopterygii_odb10.hmm
```

## Smoke tests (`~/scratch/ga-repeats/tool_tests`, 2026-11-23)

- **TEtrimmer**, its repository's `tests/test_input.fa` (9 families) on `tests/test_genome.fa`, `-t 8 --classify_all`:
  13 min 24 s, of which 9 min were the Pfam download. Output: 6 families. Its stored standard output
  (`Standard_TEtrimmer_output_for_test_file.zip`) has 5; all 5 are in ours with no class in conflict; ours classifies
  two the standard leaves Unknown (LTR/Gypsy, LTR/Copia) and keeps one more (DNA/TcMar-Fot1). The standard was made by
  an older version (no evaluation level in its headers).
- **MCHelper**, its README's toy library (`Test_dir/fasta_input/TE_lib_toy.fa`, 21 families) on the fruit-fly genome
  with BUSCO `diptera_odb10.2020-08-05`, `-r A -t 8 -a F`: 20 min 56 s; 11 families out. The repository's expected
  results belong to its full test inputs (326 families from the full fruit-fly library, 451 from REPET), not to the toy
  library; MCHelper's clustering and subfamily splitting depend on the whole library, so the toy output cannot be compared
  with them. A run of the full fruit-fly library against the 326 is an optional check.

## Output formats the import handles

- MCHelper: `curated_sequences_NR.fa`, labels in Wicker's form (`#CLASSI/LTR/GYPSY`, `#CLASSII/TIR/P`), `_inc` for
  families judged incomplete, `_s_N` for subfamilies split off.
- TEtrimmer 1.7.2: `TEtrimmer_consensus_merged.fasta`, headers `name#Class/Subclass<TAB>TEtrimmer<TAB><evaluation>`;
  it writes hyphens in names as underscores (`rnd-1_family-1` becomes `rnd_1_family_1`).
