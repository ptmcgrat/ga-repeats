#!/usr/bin/env python3
"""Build one ecogroup's library with RepeatModeler, exactly as the pilot libraries RI_v1-RI_v3 were built (spec section 4).

The pilot's method (as write_manifest.py records it), in one script, run from the work folder <scratch>/<version>/:
  1. the genome decompressed to <build>.fna
  2. BuildDatabase -name <build> <build>.fna                                           -> builddatabase.log
  3. nice -n 10 RepeatModeler -database <build> -threads <N> -LTRStruct (TETools 2.01)  -> run.out
  4. the five outputs copied, never overwriting, into <version folder>/library/: <build>-families.fa,
     <build>-families.stk, <build>-rmod.log, run.out, builddatabase.log
  5. <build>-families.prefixed.fa beside them: every family name prefixed with "<build>_" (RepeatModeler's own file is
     kept as produced)
  6. build_record.json in the work folder: assembly, accession, commands, start and finish, hours, threads, seed, families
RepeatModeler classifies against the Dfam 4.0 content cached in ~/.RepeatModelerCache during RI_v1, the same for every
build. A reference genome is refused (spec section 6, 2026-11-23). The decompressed genome is removed afterwards; the
RepeatModeler work folder (RM_*) is kept in scratch.

A run takes about 36 hours: start it inside screen. After it: mask_genomes.py, then write_manifest.py once the
version's masking is complete.

Usage (from the repository root, in ga-repeats):
  python scripts/build_library.py --version RI_v11 --assembly MZ_NMK2024a          # dry run
  python scripts/build_library.py --version RI_v11 --assembly MZ_NMK2024a --go     # inside screen
"""
import argparse
import json
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

OUTPUTS = ["{b}-families.fa", "{b}-families.stk", "{b}-rmod.log", "run.out", "builddatabase.log"]


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the new version, e.g. RI_v11")
    parser.add_argument("--assembly", required=True, help="the build genome, an id of the panel, e.g. MZ_NMK2024a")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def prefix_names(lines, build: str):
    """RepeatModeler's families file -> the same with every name prefixed '<build>_'."""
    for line in lines:
        yield (">" + build + "_" + line[1:]) if line.startswith(">") else line


def seed_of(text: str) -> Optional[str]:
    m = re.search(r"Random Number Seed: (\d+)", text)
    return m.group(1) if m else None


def main(argv=None) -> None:
    args = parse_args(argv)
    build = args.assembly
    if c.is_reference(build):
        c.stop(f"{build} is a reference genome: kept for comparisons only, never used to build a library")
    accession = c.panel_accession(build)
    genome = c.genome_path(build, accession)
    if not genome.exists():
        c.stop(f"{genome} not found")
    target = c.version_folder(args.version)
    if target.exists():
        c.stop(f"{target} already exists; a version is never changed")
    work = (args.scratch / args.version).resolve()
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; look at it, then clear it with: rm -r {work}")
    fasta = f"{build}.fna"
    bd = c.container(args.sif, "BuildDatabase", "-name", build, fasta, nice=False)
    rm = c.container(args.sif, "RepeatModeler", "-database", build, "-threads", str(args.threads), "-LTRStruct")
    print(f"{args.version}: RepeatModeler on {build} ({accession}), {args.threads} threads, in {work}\n"
          f"then into {target / 'library'}/")
    print(f"  cd {work}\n  (the genome decompressed to {fasta})")
    c.run(bd, False, work / "builddatabase.log")
    c.run(rm, False, work / "run.out")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to build (inside screen: about 36 hours).")
        return

    work.mkdir(parents=True)
    os.chdir(work)
    started = datetime.now()
    t0 = time.time()
    print(f"decompressing {genome} ...", flush=True)
    c.decompress(genome, work / fasta)
    c.run(bd, True, work / "builddatabase.log")
    c.run(rm, True, work / "run.out")
    hours = (time.time() - t0) / 3600
    names = [p.format(b=build) for p in OUTPUTS]
    missing = [n for n in names if not (work / n).exists() or (work / n).stat().st_size == 0]
    if missing:
        c.stop(f"RepeatModeler finished but these are missing or empty: {', '.join(missing)}; see {work / 'run.out'}")
    lib = target / "library"
    for n in names:
        c.copy_new(work / n, lib / n)
    original = lib / f"{build}-families.fa"
    with open(original) as src, open(lib / f"{build}-families.prefixed.fa", "x") as dst:
        dst.writelines(prefix_names(src, build))
    (work / fasta).unlink()
    families = sum(1 for line in open(original) if line.startswith(">"))
    seed = seed_of((work / f"{build}-rmod.log").read_text())
    record = {"version": args.version, "assembly": build, "accession": accession, "genome": str(genome),
              "commands": [f"cd {work}", f"zcat {genome} > {fasta}", c.shell_text(bd) + " > builddatabase.log 2>&1",
                           c.shell_text(rm) + " > run.out 2>&1"],
              "threads": args.threads, "random_seed": seed, "families": families,
              "started": f"{started:%Y-%m-%d %H:%M:%S}", "finished": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "hours": round(hours, 2),
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(work / "build_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"\n{args.version}: {families:,} families from {build} in {hours:.1f} h (seed {seed})")
    print(f"library: {lib}/\nnext: python scripts/mask_genomes.py --version {args.version} --assemblies AC_fAstCal1_2 --go")


if __name__ == "__main__":
    main()
