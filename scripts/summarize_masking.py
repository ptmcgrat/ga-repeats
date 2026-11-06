#!/usr/bin/env python3
"""What the differences between two libraries mean in the genome (question 1 pilot).

Reads a comparison made by compare_libraries.py (<A>_vs_<B>/) and the masking results
in each version's masked/<assembly_id>/ (made by mask_genomes.py). For every library
(A, B) and every genome it has masked, it reports per class:
  - base pairs masked;
  - base pairs masked by families found in only that library at the chosen threshold
    (their cluster holds no family of the other library);
  - genome size and the fractions.
Rows for DT_004f give the held-out genome's masked fraction under each library.

Example:
  python scripts/summarize_masking.py --comparison ~/scratch/ga-repeats/comparisons/RI_v1_vs_RI_v2 --threshold 0.80
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--comparison", type=Path, required=True,
                        help="a folder made by compare_libraries.py, e.g. .../comparisons/RI_v1_vs_RI_v2")
    parser.add_argument("--threshold", type=float, default=0.80, help="which clustering to use (default 0.80)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def only_in_one(cluster_tsv: Path) -> Dict[str, Set[str]]:
    """source -> family names (part before '#') whose cluster holds no other source."""
    clusters: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    with open(cluster_tsv, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            clusters[row["cluster"]].append((row["source"], c.split_name(row["family"])[0]))
    result: Dict[str, Set[str]] = defaultdict(set)
    for members in clusters.values():
        sources = {s for s, _ in members}
        if len(sources) == 1:
            result[next(iter(sources))].update(f for _, f in members)
    return result


def read_bed(path: Path) -> List[Tuple[str, int, int, str, str]]:
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3], c.top_class(f[6])))
    return rows


def summarize(bed: List[Tuple[str, int, int, str, str]], unique: Set[str]) -> Dict[str, Tuple[int, int]]:
    by_class: Dict[str, list] = defaultdict(list)
    by_class_unique: Dict[str, list] = defaultdict(list)
    for name, start, end, family, cls in bed:
        for key in (cls, "ALL"):
            by_class[key].append((name, start, end))
            if family in unique:
                by_class_unique[key].append((name, start, end))
    return {key: (c.union_length(by_class[key]),
                  c.union_length(by_class_unique[key]) if by_class_unique[key] else 0)
            for key in by_class}


def main(argv=None) -> None:
    args = parse_args(argv)
    meta = json.loads((args.comparison / "meta.json").read_text())
    cluster_tsv = args.comparison / f"clusters_c{args.threshold:.2f}.tsv"
    if not cluster_tsv.exists():
        c.stop(f"no clustering at {args.threshold:.2f} in {args.comparison}")
    target = args.comparison / f"masking_summary_c{args.threshold:.2f}.tsv"
    if target.exists():
        c.stop(f"{target} already exists")
    unique = only_in_one(cluster_tsv)

    plan = []
    for version, build_id in ((meta["a"], meta["id_a"]), (meta["b"], meta["id_b"])):
        masked = c.version_folder(version) / "masked"
        for folder in sorted(p for p in masked.glob("*") if p.is_dir()):
            bed = folder / f"{folder.name}.bed"
            if bed.exists():
                plan.append((version, build_id, folder.name, bed))
    if not plan:
        c.stop("no masking results found yet; run mask_genomes.py first")
    for version, build_id, assembly, bed in plan:
        print(f"{version} ({build_id}) on {assembly}: {bed}")
    if not args.go:
        print("dry run: nothing was written. Add --go to run.")
        return

    genome_bp: Dict[str, int] = {}
    rows = []
    for version, build_id, assembly, bed in plan:
        if assembly not in genome_bp:
            genome = c.genome_path(assembly, c.panel_accession(assembly))
            genome_bp[assembly] = sum(c.fasta_lengths(genome).values())
        size = genome_bp[assembly]
        stats = summarize(read_bed(bed), unique.get(build_id, set()))
        for cls in sorted(k for k in stats if k != "ALL") + ["ALL"]:
            masked_bp, unique_bp = stats[cls]
            rows.append([version, build_id, assembly, cls, masked_bp, unique_bp, size,
                         f"{100 * masked_bp / size:.3f}", f"{100 * unique_bp / size:.3f}"])
    c.write_tsv(target, ["library_version", "library_built_from", "assembly", "class", "masked_bp",
                         "masked_bp_by_families_only_in_this_library", "genome_bp", "percent_masked",
                         "percent_masked_by_families_only_in_this_library"], rows)
    print(f"written: {target}")
    for row in rows:
        if row[3] == "ALL":
            print(f"  {row[0]} on {row[2]}: {row[7]}% masked, {row[8]}% by families only in {row[1]}")


if __name__ == "__main__":
    main()
