#!/usr/bin/env python3
"""Mask genomes with one version's library (gate 4; spec sections 3 and 5).

For each assembly, in <scratch>/<version>/masking/<assembly_id>/:
  1. decompress the genome (names unchanged);
  2. run RepeatMasker with the version's prefixed library only (-lib; no -species,
     so no Dfam or Repbase families are added), RMBlast, -no_is;
  3. convert the .out to BED and check that every BED sequence name is in the genome
     and every interval lies within that sequence (a masked genome can then always
     be regenerated, spec section 5);
  4. only then copy the .out and the BED into <version folder>/masked/<assembly_id>/.
Masked genomes are not stored (decision 2026-10-30). Nothing is ever overwritten.

RepeatMasker's -pa: with RMBlast each job uses 4 cores (RepeatMasker -help, 1.3), so
--threads 16 gives -pa 4 per assembly; with --max-jobs 2, each assembly gets half.

Example (runs for hours: start it inside screen):
  python scripts/mask_genomes.py --version RI_v1 --assemblies MZ_MZ4f DT_004f        # dry run
  python scripts/mask_genomes.py --version RI_v1 --assemblies MZ_MZ4f DT_004f --go
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="library version, e.g. RI_v1")
    parser.add_argument("--assemblies", nargs="+", required=True, help="panel assembly ids to mask")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def mask_one(args: argparse.Namespace, assembly: str, library: Path, pa: int) -> str:
    accession = c.panel_accession(assembly)
    genome = c.genome_path(assembly, accession)
    work = args.scratch / args.version / "masking" / assembly
    dest = c.version_folder(args.version) / "masked" / assembly
    fasta = work / f"{assembly}.fna"
    command = c.container(args.sif, "RepeatMasker", "-lib", str(library), "-engine", "rmblast",
                          "-pa", str(pa), "-no_is", "-dir", str(work), str(fasta))
    print(f"[{assembly}] genome {genome} ({accession})")
    print(f"[{assembly}] work in {work}; results to {dest}")
    if not genome.exists():
        c.stop(f"[{assembly}] genome not found: {genome}")
    if work.exists():
        c.stop(f"[{assembly}] {work} already exists")
    if dest.exists():
        c.stop(f"[{assembly}] {dest} already exists; masked results are never overwritten")
    if not args.go:
        print(f"[{assembly}] would decompress the genome to {fasta}")
        c.run(command, go=False, log=work / "repeatmasker.log")
        return f"{assembly}: dry run"

    work.mkdir(parents=True)
    started = datetime.now()
    c.decompress(genome, fasta)
    lengths = c.fasta_lengths(fasta)
    c.run(command, go=True, log=work / "repeatmasker.log")
    out_file = work / f"{assembly}.fna.out"
    if not out_file.exists():
        c.stop(f"[{assembly}] RepeatMasker wrote no {out_file.name}; see {work / 'repeatmasker.log'}")
    hits = list(c.parse_out(out_file))
    if not hits:
        c.stop(f"[{assembly}] {out_file.name} holds no annotations")
    bed_rows = [c.hit_to_bed(h) for h in hits]
    problems = c.check_bed(bed_rows, lengths)
    if problems:
        c.stop(f"[{assembly}] BED check failed:\n  " + "\n  ".join(problems))
    bed_file = work / f"{assembly}.bed"
    c.write_tsv(bed_file, ["#chrom", "start", "end", "family", "score", "strand", "class_family",
                           "divergence", "sw_score", "rm_id"], bed_rows)
    masked_bp = c.union_length((r[0], int(r[1]), int(r[2])) for r in bed_rows)
    genome_bp = sum(lengths.values())

    c.copy_new(out_file, dest / f"{assembly}.out")
    c.copy_new(bed_file, dest / f"{assembly}.bed")
    header = [line.strip() for line in (work / "repeatmasker.log").read_text().splitlines()[:40]
              if line.strip()]
    record = {
        "version": args.version, "assembly_id": assembly, "accession_version": accession,
        "genome": str(genome), "library": str(library), "library_sha256": c.sha256(library),
        "command": c.shell_text(command), "started": f"{started:%Y-%m-%d %H:%M}",
        "finished": f"{datetime.now():%Y-%m-%d %H:%M}",
        "repeatmasker_log_head": header,
        "annotations": len(hits), "masked_bp": masked_bp, "genome_bp": genome_bp,
        "bed_check": "passed: every BED name is a genome sequence and every interval lies within it",
        "outputs": [str(dest / f"{assembly}.out"), str(dest / f"{assembly}.bed")],
    }
    (work / "run_record.json").write_text(json.dumps(record, indent=2) + "\n")
    return (f"{assembly}: {len(hits)} annotations, {masked_bp:,} of {genome_bp:,} bp masked "
            f"({100 * masked_bp / genome_bp:.2f}%); BED check passed; copied to {dest}")


def main(argv=None) -> None:
    args = parse_args(argv)
    _, _, library = c.library_files(args.version)
    jobs = max(1, min(args.max_jobs, len(args.assemblies)))
    pa = max(1, (args.threads // jobs) // c.RMBLAST_CORES_PER_JOB)
    print(f"{args.version}: library {library}")
    print(f"{len(args.assemblies)} assemblies, {jobs} at a time, RepeatMasker -pa {pa} each "
          f"(about {jobs * pa * c.RMBLAST_CORES_PER_JOB} cores in total)")
    if len(set(args.assemblies)) != len(args.assemblies):
        c.stop("an assembly is listed twice")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results = list(pool.map(lambda a: mask_one(args, a, library, pa), args.assemblies))
    print("\n".join(results))
    if not args.go:
        print("dry run: nothing was written. Add --go to run.")


if __name__ == "__main__":
    main()
