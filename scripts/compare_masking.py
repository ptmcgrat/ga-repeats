#!/usr/bin/env python3
"""Second check: on one genome, which bases does one library mask and another not?

"Masked by families found in only one library" (summarize_masking.py) is not the same
as "sequence the other library misses": the other library may mask the same bases
under other family names. This script compares the masked bases themselves. For every
pair of versions that have masked the assembly, it reports:
  - bases masked by both, by only the first, by only the second (whole genome);
  - per class of each library: its masked bases, and how many of them the other
    library does not mask at all.

Inputs: <version folder>/masked/<assembly>/<assembly>.bed (made by mask_genomes.py).
Result: <scratch>/comparisons/masking_overlap/<assembly>.tsv. Never overwritten.

Example:
  python scripts/compare_masking.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3        # dry run
  python scripts/compare_masking.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3 --go
"""
from __future__ import annotations

import argparse
import itertools
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

MAIN_CLASSES = ("DNA", "LINE", "LTR", "SINE", "Unknown")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly", required=True, help="a genome masked by every version listed, e.g. DT_004f")
    parser.add_argument("--versions", nargs="+", required=True, help="two or more versions, e.g. RI_v1 RI_v2 RI_v3")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def masked_sets(rows) -> Tuple[Dict[str, list], Dict[str, Dict[str, list]]]:
    """Merged intervals of all annotations, and per class."""
    by_class: Dict[str, list] = defaultdict(list)
    for name, start, end, _family, cls in rows:
        by_class[cls].append((name, start, end))
    everything = c.merge_by_sequence((r[0], r[1], r[2]) for r in rows)
    return everything, {cls: c.merge_by_sequence(spans) for cls, spans in by_class.items()}


def pair_rows(assembly: str, a: str, b: str, sets: Dict[str, tuple], genome_bp: int) -> List[list]:
    all_a, classes_a = sets[a]
    all_b, classes_b = sets[b]
    rows = []
    for side, mine, other_all, my_classes in ((a, all_a, all_b, classes_a), (b, all_b, all_a, classes_b)):
        for cls in sorted(my_classes) + ["ALL"]:
            spans = mine if cls == "ALL" else my_classes[cls]
            masked = c.merged_length(spans)
            alone = c.subtract_length(spans, other_all)
            rows.append([assembly, f"{a} vs {b}", side, cls, masked, alone, genome_bp,
                         f"{100 * masked / genome_bp:.3f}", f"{100 * alone / genome_bp:.3f}"])
    return rows


def main(argv=None) -> None:
    args = parse_args(argv)
    if len(set(args.versions)) < 2:
        c.stop("give at least two different versions")
    beds = {}
    for version in args.versions:
        bed = c.version_folder(version) / "masked" / args.assembly / f"{args.assembly}.bed"
        if not bed.exists():
            c.stop(f"{version} has not masked {args.assembly}: {bed} not found")
        beds[version] = bed
        print(f"{version}: {bed}")
    out_dir = args.scratch / "comparisons" / "masking_overlap"
    target = out_dir / f"{args.assembly}.tsv"
    pairs = list(itertools.combinations(args.versions, 2))
    print(f"pairs: {', '.join(f'{a} vs {b}' for a, b in pairs)}")
    print(f"result: {target}")
    if target.exists():
        c.stop(f"{target} already exists; results are never overwritten")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    genome_bp = sum(c.fasta_lengths(genome).values())
    sets = {}
    for version, bed in beds.items():
        print(f"reading {version} ...", flush=True)
        sets[version] = masked_sets(c.read_bed(bed))

    rows = []
    for a, b in pairs:
        rows += pair_rows(args.assembly, a, b, sets, genome_bp)
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["assembly", "pair", "library", "class", "masked_bp", "masked_bp_not_masked_by_other",
                         "genome_bp", "percent_masked", "percent_not_masked_by_other"], rows)
    print(f"written: {target}\n")
    for a, b in pairs:
        all_rows = {r[2]: r for r in rows if r[1] == f"{a} vs {b}" and r[3] == "ALL"}
        only_a, only_b = all_rows[a][5], all_rows[b][5]
        both = all_rows[a][4] - only_a
        print(f"{args.assembly}, {a} vs {b}: masked by both {100 * both / genome_bp:.2f}%, "
              f"only by {a} {100 * only_a / genome_bp:.2f}%, only by {b} {100 * only_b / genome_bp:.2f}%")
        for side in (a, b):
            parts = [f"{r[3]} {r[8]}%" for r in rows
                     if r[1] == f"{a} vs {b}" and r[2] == side and r[3] in MAIN_CLASSES]
            print(f"    not masked by the other, by {side}'s class: " + ", ".join(parts))


if __name__ == "__main__":
    main()
