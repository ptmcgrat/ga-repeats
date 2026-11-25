#!/usr/bin/env python3
"""Curate a library with an existing tool, MCHelper or TEtrimmer (question 2: existing tools before our own rules).

Runs one tool on one version's library against a genome, by default the library's build genome, in the tool's own
conda environment (spec section 6, 2026-11-23): repeat_family_merging_mchelper or repeat_family_merging_tetrimmer,
called through `conda run`; this script itself runs in ga-repeats. With --pilot N it runs on N families only, to time
the tool and look at its output before a full run: the --include families present in the library (by default the
composite families at the tbx2a and mitfa insertions and the satellite family) plus a random sample of the rest,
stratified by class (fixed --seed).

The commands (both tools' own --help, 2026-11-23):
  MCHelper 1.7.0   python3 MCHelper.py -r A -a F --te_aid N -t <threads> --input_type fasta -l <library>
                   -g <genome> -o <out> -b <BUSCO HMMs>
                   Everything else is MCHelper's default: 16 extension rounds of 500 bp, and a family is processed
                   only with at least one full-length fragment (a copy covering 94% of the consensus). TE-Aid draws
                   plots for manual inspection only (MCHelper's README), so it is off.
  TEtrimmer 1.7.2  TEtrimmer -i <library> -g <genome> -o <out> -t <threads> --classify_all -pd <Pfam folder> -l <log>
                   Everything else is TEtrimmer's default.

Usage (from the repository root, in ga-repeats):
  python scripts/curate_library.py --tool tetrimmer --version RI_v1 --pilot 100          # dry run
  python scripts/curate_library.py --tool tetrimmer --version RI_v1 --pilot 100 --go

Work, never overwritten: <scratch>/curation/<tool>/<version>_<assembly>[_pilot<N>]/ with input.fa, out/ (the tool's
output), tool.log, and record.json (tool and version, environment, command, families in and out, times). The genome
copy and anything the tool indexes beside it are removed after the run.
"""
import argparse
import json
import os
import random
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

TOOLS = ("mchelper", "tetrimmer")
INCLUDE = ["MZ_MZ4f_ltr-1_family-769", "MZ_MZ4f_ltr-1_family-276", "MZ_NMetZ_ltr-1_family-392",
           "MZ_MZ4f_rnd-5_family-9801"]
OUTPUT = {"mchelper": Path("out") / "curated_sequences_NR.fa", "tetrimmer": Path("out") / "TEtrimmer_consensus_merged.fasta"}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tool", required=True, choices=TOOLS)
    parser.add_argument("--version", required=True, help="the library to curate, e.g. RI_v1")
    parser.add_argument("--assembly", help="genome to curate against (default: the version's build genome; "
                                           "required for a merged version)")
    parser.add_argument("--pilot", type=int, help="curate only this many families (stratified sample)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--include", nargs="*", default=INCLUDE, help="families always in a pilot, where present")
    parser.add_argument("--env", help="the tool's conda environment (default repeat_family_merging_<tool>)")
    parser.add_argument("--conda", default=os.environ.get("CONDA_EXE", "conda"), help="conda executable")
    parser.add_argument("--mchelper-dir", type=Path, help="MCHelper's repository (default <scratch>/tools/MCHelper)")
    parser.add_argument("--busco", type=Path,
                        help="BUSCO HMMs for MCHelper (default <scratch>/tools/busco/actinopterygii_odb10.hmm)")
    parser.add_argument("--pfam-dir", type=Path, help="TEtrimmer's Pfam folder (default <scratch>/tools/tetrimmer_pfam)")
    c.add_common_args(parser)
    args = parser.parse_args(argv)
    args.env = args.env or f"repeat_family_merging_{args.tool}"
    args.mchelper_dir = args.mchelper_dir or args.scratch / "tools" / "MCHelper"
    args.busco = args.busco or args.scratch / "tools" / "busco" / "actinopterygii_odb10.hmm"
    args.pfam_dir = args.pfam_dir or args.scratch / "tools" / "tetrimmer_pfam"
    return args


def stratified_sample(families: Sequence[c.Family], n: int, include: Sequence[str], seed: int) -> List[c.Family]:
    """The include families present (by name before '#'), plus a class-stratified random sample of the rest, n in all."""
    wanted = set(include)
    chosen = [f for f in families if f.family in wanted]
    rest = [f for f in families if f.family not in wanted]
    k = max(0, min(n - len(chosen), len(rest)))
    groups: Dict[str, List[c.Family]] = defaultdict(list)
    for f in sorted(rest, key=lambda f: f.name):
        groups[la.group(f.label)].append(f)
    quota = {g: k * len(v) / len(rest) for g, v in groups.items()} if rest else {}
    take = {g: int(q) for g, q in quota.items()}
    for g in sorted(quota, key=lambda g: (-(quota[g] - take[g]), g))[: k - sum(take.values())]:
        take[g] += 1
    rng = random.Random(seed)
    for g in sorted(groups):
        chosen += rng.sample(groups[g], take.get(g, 0))
    return chosen


SUFFIXES = (r"_inc$", r"_s_\d+$", r"_\d{2}$")     # MCHelper: incomplete, split subfamily; TEtrimmer: subfamily _01


def normalize(name: str) -> str:
    """Family name before '#', with hyphens as underscores (TEtrimmer writes rnd-1_family-1 as rnd_1_family_1)."""
    return name.split("#", 1)[0].replace("-", "_")


def source_family(name: str, inputs: Dict[str, str]) -> Optional[str]:
    """The input family an output family derives from: an exact match first, then with the tools' suffixes removed one
    at a time (so a family numbered _20 is never mistaken for TEtrimmer's subfamily suffix)."""
    current = normalize(name)
    while True:
        if current in inputs:
            return inputs[current]
        for pattern in SUFFIXES:
            stripped = re.sub(pattern, "", current)
            if stripped != current:
                current = stripped
                break
        else:
            return None


def rm_class(tool: str, label: str) -> str:
    """An output label as a RepeatMasker class (MCHelper's Wicker labels mapped)."""
    return c.top_class(la.wicker_label(label) if tool == "mchelper" else label)


def tool_command(args, input_fa: Path, genome: Path, out: Path, work: Path) -> List[str]:
    base = ["nice", "-n", "10", args.conda, "run", "-n", args.env, "--no-capture-output"]
    if args.tool == "mchelper":
        return base + ["python3", str(args.mchelper_dir / "MCHelper.py"), "-r", "A", "-a", "F", "--te_aid", "N",
                       "-t", str(args.threads), "--input_type", "fasta", "-l", str(input_fa), "-g", str(genome),
                       "-o", str(out), "-b", str(args.busco)]
    return base + ["TEtrimmer", "-i", str(input_fa), "-g", str(genome), "-o", str(out), "-t", str(args.threads),
                   "--classify_all", "-pd", str(args.pfam_dir), "-l", str(work / "tetrimmer.log")]


def version_command(args) -> List[str]:
    base = [args.conda, "run", "-n", args.env, "--no-capture-output"]
    if args.tool == "mchelper":
        return base + ["python3", str(args.mchelper_dir / "MCHelper.py"), "--version"]
    return base + ["TEtrimmer", "--version"]


def tool_version(args) -> str:
    """The tool's own report of its version (last non-empty line), or why it could not be read."""
    try:
        result = subprocess.run(version_command(args), capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired) as error:
        c.stop(f"cannot run {args.tool} in environment {args.env}: {error}")
    lines = [l.strip() for l in (result.stdout + result.stderr).splitlines() if l.strip()]
    if result.returncode != 0 or not lines:
        c.stop(f"{args.tool} in environment {args.env} did not report its version:\n" + "\n".join(lines[-5:]))
    return lines[-1]


def main(argv=None) -> None:
    args = parse_args(argv)
    build, _, library = c.library_files(args.version)
    assembly = args.assembly or (None if build == args.version else build)
    if not assembly:
        c.stop(f"{args.version} is a merged or derived version: give the genome to curate against with --assembly")
    if c.is_reference(assembly):
        c.stop(f"{assembly} is a reference genome: kept for comparisons only, never used to build or curate a library "
               f"(spec section 6, 2026-11-23)")
    genome = c.genome_path(assembly, c.panel_accession(assembly))
    needs = [args.mchelper_dir / "MCHelper.py", args.busco] if args.tool == "mchelper" else [args.pfam_dir]
    for f in needs + [genome]:
        if not f.exists():
            c.stop(f"{f} not found")
    families = c.read_families(library)
    chosen = stratified_sample(families, args.pilot, args.include, args.seed) if args.pilot else list(families)
    stem = f"{args.version}_{assembly}" + (f"_pilot{args.pilot}" if args.pilot else "")
    work = args.scratch / "curation" / args.tool / stem
    if work.exists():
        c.stop(f"{work} exists; results are never overwritten. To run again: rm -r {work}")
    input_fa, fasta, out = work / "input.fa", work / f"{assembly}.fna", work / "out"
    command = tool_command(args, input_fa, fasta, out, work)
    by_class = Counter(la.group(f.label) for f in chosen)
    print(f"{args.tool} on {args.version} ({len(chosen):,} of {len(families):,} families: "
          f"{', '.join(f'{k} {v}' for k, v in by_class.most_common())}) against {assembly}, {args.threads} threads, "
          f"in environment {args.env}")
    if args.pilot:
        present = [f.family for f in chosen if f.family in set(args.include)]
        print(f"pilot: always included {len(present)} ({', '.join(present) or 'none present'}); the rest sampled "
              f"by class with seed {args.seed}")
    print(f"work: {work}")
    c.run(command, False, work / "tool.log")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    version_text = tool_version(args)
    print(f"{args.tool} reports: {version_text}", flush=True)
    work.mkdir(parents=True)
    names = {f.name for f in chosen}
    with open(library) as src, open(input_fa, "x") as dst:
        keep = False
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split()[0] in names
            if keep:
                dst.write(line)
    c.decompress(genome, fasta)
    started = datetime.now()
    t0 = time.time()
    c.run(command, True, work / "tool.log")
    hours = (time.time() - t0) / 3600
    for f in work.glob(fasta.name + "*"):
        f.unlink()
    result = work / OUTPUT[args.tool]
    if not result.exists():
        c.stop(f"{args.tool} finished but wrote no {result}; see {work / 'tool.log'}")
    out_families = c.read_families(result)
    labels = Counter(rm_class(args.tool, f.label) for f in out_families)
    inputs = {normalize(f.family): f.family for f in chosen}
    sources = Counter(source_family(f.family, inputs) for f in out_families)
    represented = {k for k in sources if k}
    record = {"tool": args.tool, "tool_version": version_text, "environment": args.env, "version": args.version,
              "assembly": assembly, "pilot": args.pilot, "seed": args.seed if args.pilot else None,
              "included": sorted(f.family for f in chosen if f.family in set(args.include)) if args.pilot else [],
              "families_in": len(chosen), "families_out": len(out_families), "labels_out": dict(labels.most_common()),
              "inputs_represented": len(represented), "inputs_split": sum(1 for k, v in sources.items() if k and v > 1),
              "inputs_absent": len(chosen) - len(represented), "outputs_without_source": sources.get(None, 0),
              "command": c.shell_text(command), "threads": args.threads, "output": str(OUTPUT[args.tool]),
              "started": f"{started:%Y-%m-%d %H:%M:%S}", "finished": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "hours": round(hours, 2),
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(work / "record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"\n{args.tool}: {len(chosen):,} families in, {len(out_families):,} out in {hours:.2f} h "
          f"({', '.join(f'{k} {v}' for k, v in labels.most_common(8))}); inputs kept {len(represented):,} "
          f"(split {record['inputs_split']:,}), absent {record['inputs_absent']:,}")
    print(f"per family: python scripts/curation_report.py --from {work} --go")
    print(f"library: {result}\nrecord:  {work / 'record.json'}")
    if not args.pilot:
        print(f"next: python scripts/import_curated.py --from {work} --version <new version>")


if __name__ == "__main__":
    main()
