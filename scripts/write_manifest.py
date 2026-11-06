#!/usr/bin/env python3
"""Write a version's MANIFEST.yaml (spec section 3; gate 5, question 20).

Collects, for one version folder:
  - the build: assembly, accession version, start and finish, random seed, threads,
    sampling (from the RepeatModeler log), family counts;
  - a check that the genome's sequence and base counts equal BuildDatabase's and the
    RepeatModeler log's (gate 1, question 8);
  - tools: the image tag and checksum, Apptainer, RepeatModeler, its dependencies,
    the LTR pipeline, the Dfam version RepeatModeler classified against, ga_common;
  - the library files with their checksums;
  - the masked set: each masked/<assembly_id>/ with its RepeatMasker command (from
    mask_genomes.py's run record) and a fresh BED check;
  - the commands, as run.

A version folder's files are never replaced, so the manifest is written once, when the
folder is complete (after its masking). Without --go it is only printed.

Example:
  python scripts/write_manifest.py --version RI_v1          # print it
  python scripts/write_manifest.py --version RI_v1 --go     # write it (once)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import socket
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

IMAGE_TAG = "docker://dfam/tetools:2.01"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="e.g. RI_v1")
    parser.add_argument("--purpose", default="question 1 pilot")
    c.add_common_args(parser)
    return parser.parse_args(argv)


# ---------------------------------------------------------------- log parsing

def parse_builddatabase(text: str) -> Tuple[int, int]:
    match = re.search(r"Number of sequences \(bp\) added to database:\s*(\d+)\s*\(\s*(\d+)\s*bp\s*\)", text)
    if not match:
        c.stop("no 'Number of sequences (bp) added to database' line in builddatabase.log")
    return int(match.group(1)), int(match.group(2))


def _one(pattern: str, text: str, flags: int = 0) -> Optional[str]:
    match = re.search(pattern, text, flags)
    return match.group(1).strip() if match else None


def parse_rmod_log(text: str) -> Dict[str, object]:
    info: Dict[str, object] = {
        "repeatmodeler": _one(r"RepeatModeler Version (\S+)", text),
        "search_engine": _one(r"Search Engine = (.+)", text),
        "threads": _one(r"Threads = (\d+)", text),
        "dependencies": _one(r"Dependencies: (.+)", text),
        "random_seed": _one(r"Random Number Seed: (\d+)", text),
        "database_sequences": _one(r"- Sequences = (\d+)", text),
        "database_bases": _one(r"- Bases = (\d+)", text),
        "round1_sample_bp": _one(r"Final Sample Size = (\d+) bp", text),
        "repeatscout_recon_families": _one(r"RepeatScout/RECON discovery complete: (\d+) families", text),
        "ltr_pipeline_families": _one(r"- (\d+) LTRPipeline families", text),
        "ltr_redundant_removed": _one(r"Removed (\d+) redundant LTR families", text),
        "program_time": _one(r"Program Time: (\S+)", text),
        "dfam": _one(r"(?:PROT|CONS)-Dfam_([\d.]+)", text),
        "rm_directory": _one(r"Using output directory = (\S+)", text),
    }
    ltr = re.search(r"LTR Structural Analysis: (\w+)\s*(\(.*?\))?", text, re.S)
    info["ltr_structural_analysis"] = (" ".join(ltr.group(0).split()[3:]) if ltr else None)
    info["rounds"] = len(re.findall(r"^RepeatModeler Round # \d+", text, re.M))
    coverage = re.findall(r"Input Database Coverage: (\d+) bp out of (\d+) bp \( ([\d.]+) % \)", text)
    info["final_coverage"] = (f"{coverage[-1][0]} of {coverage[-1][1]} bp ({coverage[-1][2]}%)"
                              if coverage else None)
    dfam_sizes = re.findall(r"added (\d+) sequences", text)
    info["dfam_reference_sizes"] = dfam_sizes or None
    return info


def run_start(rm_directory: Optional[str]) -> Optional[datetime]:
    """RepeatModeler names its working folder RM_<pid>.<Day><Mon><dd><HH><MM><SS><YYYY>."""
    if not rm_directory:
        return None
    match = re.search(r"RM_\d+\.(\w{3})(\w{3})(\d{2})(\d{2})(\d{2})(\d{2})(\d{4})", rm_directory)
    if not match:
        return None
    _, month, day, hh, mm, ss, year = match.groups()
    return datetime.strptime(f"{year} {month} {day} {hh}:{mm}:{ss}", "%Y %b %d %H:%M:%S")


def add_duration(start: Optional[datetime], hms: Optional[str]) -> Optional[datetime]:
    if not start or not hms:
        return None
    h, m, s = (int(x) for x in hms.split(":"))
    return start + timedelta(hours=h, minutes=m, seconds=s)


# ---------------------------------------------------------------- YAML

def _scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def to_yaml(value: object, indent: int = 0) -> List[str]:
    pad = "  " * indent
    lines: List[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, (dict, list)) and item:
                lines.append(f"{pad}{key}:")
                lines += to_yaml(item, indent + 1)
            else:
                lines.append(f"{pad}{key}: {_scalar(item) if not isinstance(item, (dict, list)) else '[]'}")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                sub = to_yaml(item, indent + 1)
                lines.append(f"{pad}- {sub[0].strip()}")
                lines += sub[1:]
            else:
                lines.append(f"{pad}- {_scalar(item)}")
    else:
        lines.append(f"{pad}{_scalar(value)}")
    return lines


# ---------------------------------------------------------------- tools

def apptainer_version() -> str:
    try:
        out = subprocess.run(["apptainer", "--version"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        out = "unknown"
    env = os.environ.get("CONDA_DEFAULT_ENV", "unknown")
    return f"{out} (conda environment {env})"


def ga_common_version() -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:
        return "unknown"
    for name in ("ga_common", "ga-common"):
        try:
            return version(name)
        except PackageNotFoundError:
            continue
    return "unknown"


# ---------------------------------------------------------------- main

def main(argv=None) -> None:
    args = parse_args(argv)
    vdir = c.version_folder(args.version)
    build_id, original, prefixed = c.library_files(args.version)
    lib = vdir / "library"
    accession = c.panel_accession(build_id)
    genome = c.genome_path(build_id, accession)
    rmod_log = lib / f"{build_id}-rmod.log"
    info = parse_rmod_log(rmod_log.read_text())
    bd_sequences, bd_bases = parse_builddatabase((lib / "builddatabase.log").read_text())

    print(f"{args.version}: reading {genome} to count sequences and bases (a minute or two)")
    genome_lengths = c.fasta_lengths(genome)
    genome_sequences, genome_bases = len(genome_lengths), sum(genome_lengths.values())
    counts_agree = (genome_sequences == bd_sequences == int(info["database_sequences"] or -1)
                    and genome_bases == bd_bases == int(info["database_bases"] or -1))

    start = run_start(info["rm_directory"])  # type: ignore[arg-type]
    finish = add_duration(start, info["program_time"])  # type: ignore[arg-type]
    families_final = len(c.fasta_lengths(original))
    sizes = info["dfam_reference_sizes"]
    threads = info["threads"]

    commands = [
        f"(by hand, session functions prep / start_run / finish_run; in {args.scratch / args.version})",
        f"zcat {genome} > {build_id}.fna",
        f"apptainer exec --cleanenv {args.sif} BuildDatabase -name {build_id} {build_id}.fna",
        f"nice -n 10 apptainer exec --cleanenv {args.sif} RepeatModeler -database {build_id} "
        f"-threads {threads} -LTRStruct > run.out 2>&1   (inside screen)",
        f"cp -n {build_id}-families.fa {build_id}-families.stk {build_id}-rmod.log run.out builddatabase.log "
        f"{lib}/",
        f"awk -v p=\"{build_id}_\" '/^>/ {{sub(/^>/, \">\" p)}} {{print}}' {original} > {prefixed}",
    ]

    masked = []
    masked_root = vdir / "masked"
    for folder in sorted(p for p in masked_root.glob("*") if p.is_dir()) if masked_root.exists() else []:
        assembly = folder.name
        record_file = args.scratch / args.version / "masking" / assembly / "run_record.json"
        record = json.loads(record_file.read_text()) if record_file.exists() else {}
        bed = folder / f"{assembly}.bed"
        lengths = c.fasta_lengths(c.genome_path(assembly, c.panel_accession(assembly)))
        with open(bed) as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]
        problems = c.check_bed(rows, lengths)
        masked.append({
            "assembly_id": assembly,
            "accession_version": c.panel_accession(assembly),
            "files": [p.name for p in sorted(folder.iterdir())],
            "command": record.get("command", "MISSING: no run record from mask_genomes.py"),
            "started": record.get("started"), "finished": record.get("finished"),
            "annotations": record.get("annotations"), "masked_bp": record.get("masked_bp"),
            "genome_bp": sum(lengths.values()),
            "bed_check": "passed" if not problems else "FAILED: " + "; ".join(problems[:3]),
        })
        commands.append(record.get("command", f"MISSING: masking command for {assembly}"))

    manifest = {
        "version": args.version,
        "purpose": args.purpose,
        "panel_version": c.PANEL_VERSION,
        "host": socket.gethostname(),
        "written": f"{datetime.now():%Y-%m-%d %H:%M}",
        "build_set": [{"assembly_id": build_id, "accession_version": accession, "genome": str(genome)}],
        "masked_set": [{"assembly_id": m["assembly_id"], "accession_version": m["accession_version"]}
                       for m in masked],
        "build": {
            "started": f"{start:%Y-%m-%d %H:%M:%S}" if start else None,
            "finished": f"{finish:%Y-%m-%d %H:%M:%S}" if finish else None,
            "program_time": info["program_time"],
            "random_seed": info["random_seed"],
            "threads": threads,
            "rounds": info["rounds"],
            "round1_sample_bp": info["round1_sample_bp"],
            "final_database_coverage": info["final_coverage"],
            "families_repeatscout_recon": info["repeatscout_recon_families"],
            "families_ltr_pipeline": info["ltr_pipeline_families"],
            "ltr_families_removed_as_redundant": info["ltr_redundant_removed"],
            "families_final": families_final,
        },
        "checks": {
            "genome": f"{genome_sequences} sequences, {genome_bases} bp",
            "builddatabase": f"{bd_sequences} sequences, {bd_bases} bp",
            "repeatmodeler_log": f"{info['database_sequences']} sequences, {info['database_bases']} bp",
            "counts_agree": counts_agree,
            "masked_bed_checks": [f"{m['assembly_id']}: {m['bed_check']}" for m in masked],
        },
        "tools": {
            "image": IMAGE_TAG,
            "image_file": str(args.sif),
            "image_sha256": c.sha256(args.sif) if args.sif.exists() else "MISSING",
            "apptainer": apptainer_version(),
            "repeatmodeler": info["repeatmodeler"],
            "search_engine": info["search_engine"],
            "dependencies": info["dependencies"],
            "ltr_structural_analysis": info["ltr_structural_analysis"],
            "dfam": (f"{info['dfam']} (FamDB in the image; RepeatModeler classification cache "
                     f"~/.RepeatModelerCache/PROT-Dfam_{info['dfam']} and CONS-Dfam_{info['dfam']})"),
            "dfam_classification_reference": (
                f"{sizes[0]} RepeatPeps proteins, {sizes[1]} consensus sequences (built during this run)"
                if sizes and len(sizes) >= 2 else "built during an earlier run; reused from the cache"),
            "ga_common": ga_common_version(),
        },
        "library": [{"file": p.name, "sha256": c.sha256(p)} for p in sorted(lib.iterdir()) if p.is_file()],
        "masked": masked,
        "commands": commands,
    }
    text = "\n".join(to_yaml(manifest)) + "\n"
    print(text)

    problems = []
    if not counts_agree:
        problems.append("the genome, BuildDatabase and RepeatModeler counts disagree")
    if not masked:
        problems.append("nothing is masked yet; the manifest is written once the version folder is complete")
    if any(not m["bed_check"].startswith("passed") for m in masked):
        problems.append("a masked BED failed its check")
    if any(str(m["command"]).startswith("MISSING") for m in masked):
        problems.append("a masking run record is missing")
    if problems:
        if args.go:
            c.stop("; ".join(problems))
        print("NOT READY TO WRITE: " + "; ".join(problems))
        return
    if not args.go:
        print("dry run: nothing was written. Add --go to write it (once).")
        return
    if (vdir / "MANIFEST.yaml").exists():
        c.stop(f"{vdir / 'MANIFEST.yaml'} already exists; a version folder's files are never replaced")
    with open(vdir / "MANIFEST.yaml", "x") as handle:
        handle.write(text)
    print(f"written: {vdir / 'MANIFEST.yaml'}")


if __name__ == "__main__":
    main()
