#!/usr/bin/env python3
"""Does adding libraries keep adding masked sequence? (briefing B4: add genomes until the gain stops)

On one genome masked by several libraries, for every order of adding them:
  step 1   the bases the first library masks;
  step k   the bases the k-th library masks that none of the earlier ones did, split by
           the k-th library's classes; and the running total of bases masked by any of them.
Pairwise comparisons (compare_masking.py) show what a second library adds; this shows
whether a third adds less. Every order is reported, because the gain depends on which
library comes first.

Inputs: <version folder>/masked/<assembly>/<assembly>.bed (made by mask_genomes.py).
Results: <scratch>/comparisons/masking_gain/<assembly>.tsv, and with --plot a figure next to
it (needs matplotlib). Never overwritten.

Example:
  python scripts/masking_gain.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3          # dry run
  python scripts/masking_gain.py --assembly DT_004f --versions RI_v1 RI_v2 RI_v3 --plot --go
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path
from typing import Dict, FrozenSet, List

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import compare_masking as cm  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly", required=True, help="a genome masked by every version listed, e.g. DT_004f")
    parser.add_argument("--versions", nargs="+", required=True, help="two to six versions, e.g. RI_v1 RI_v2 RI_v3")
    parser.add_argument("--plot", action="store_true", help="also draw the gain curve (needs matplotlib)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def union(sets: List[Dict[str, list]]) -> Dict[str, list]:
    return c.merge_by_sequence((name, s, e) for spans in sets for name, items in spans.items() for s, e in items)


def gain_rows(assembly: str, sets: Dict[str, tuple], genome_bp: int) -> List[list]:
    """sets: version -> (all masked, merged; per class, merged). One row per order, step and class."""
    versions = list(sets)
    unions: Dict[FrozenSet[str], Dict[str, list]] = {frozenset(): {}}

    def union_of(members: FrozenSet[str]) -> Dict[str, list]:
        if members not in unions:
            unions[members] = union([sets[v][0] for v in sorted(members)])
        return unions[members]

    rows = []
    for order in itertools.permutations(versions):
        label = " > ".join(order)
        for step, version in enumerate(order, start=1):
            before = union_of(frozenset(order[:step - 1]))
            after = union_of(frozenset(order[:step]))
            cumulative = c.merged_length(after)
            all_spans, by_class = sets[version]
            for cls in sorted(by_class):
                added = c.subtract_length(by_class[cls], before)
                rows.append([assembly, label, step, version, cls, added, "", genome_bp,
                             f"{100 * added / genome_bp:.3f}", ""])
            added = c.subtract_length(all_spans, before)
            rows.append([assembly, label, step, version, "ALL", added, cumulative, genome_bp,
                         f"{100 * added / genome_bp:.3f}", f"{100 * cumulative / genome_bp:.3f}"])
    return rows


def step_means(rows: List[list]) -> Dict[int, float]:
    """Mean percent added at each step, over all orders (ALL rows)."""
    by_step: Dict[int, List[float]] = {}
    for r in rows:
        if r[4] == "ALL":
            by_step.setdefault(r[2], []).append(float(r[8]))
    return {step: sum(v) / len(v) for step, v in sorted(by_step.items())}


def main(argv=None) -> None:
    args = parse_args(argv)
    if not 2 <= len(set(args.versions)) <= 6 or len(set(args.versions)) != len(args.versions):
        c.stop("give two to six different versions")
    beds = {}
    for version in args.versions:
        bed = c.version_folder(version) / "masked" / args.assembly / f"{args.assembly}.bed"
        if not bed.exists():
            c.stop(f"{version} has not masked {args.assembly}: {bed} not found")
        beds[version] = bed
        print(f"{version}: {bed}")
    out_dir = args.scratch / "comparisons" / "masking_gain"
    target = out_dir / f"{args.assembly}.tsv"
    figures = [out_dir / f"{args.assembly}_gain.{ext}" for ext in ("png", "pdf")] if args.plot else []
    n_orders = len(list(itertools.permutations(args.versions)))
    print(f"{n_orders} orders of adding {len(args.versions)} libraries; result: {target}"
          + (f"; figure: {figures[0].with_suffix('')}.png and .pdf" if figures else ""))
    existing = [p for p in [target] + figures if p.exists()]
    if existing:
        c.stop("already exists, never overwritten: " + ", ".join(map(str, existing)))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    genome_bp = sum(c.fasta_lengths(genome).values())
    sets = {}
    for version, bed in beds.items():
        print(f"reading {version} ...", flush=True)
        sets[version] = cm.masked_sets(c.read_bed(bed))
    rows = gain_rows(args.assembly, sets, genome_bp)
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["assembly", "order", "step", "library_added", "class", "added_bp", "cumulative_bp",
                         "genome_bp", "percent_added", "percent_cumulative"], rows)
    print(f"written: {target}\n")

    for order in dict.fromkeys(r[1] for r in rows):
        steps = [r for r in rows if r[1] == order and r[4] == "ALL"]
        print(f"{order}: " + " -> ".join(f"{r[3]} +{r[8]}% (total {r[9]}%)" for r in steps))
    means = step_means(rows)
    print("\nmean added, over all orders: " + ", ".join(f"step {s}: {m:.2f}%" for s, m in means.items()))
    last = max(means)
    print("classes of the bases added at the last step (mean over orders): " + ", ".join(
        f"{cls} {sum(float(r[8]) for r in rows if r[2] == last and r[4] == cls) / n_orders:.3f}%"
        for cls in ("DNA", "LINE", "LTR", "SINE", "Unknown")))

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            c.stop("matplotlib is not installed in this environment; the table was written, the figure was not")
        fig, ax = plt.subplots(figsize=(6, 4))
        for order in dict.fromkeys(r[1] for r in rows):
            steps = [r for r in rows if r[1] == order and r[4] == "ALL"]
            ax.plot([r[2] for r in steps], [float(r[9]) for r in steps], color="#b4b2a9", linewidth=1, marker="o",
                    markersize=3)
        mean_total, total = [], 0.0
        for s, m in means.items():
            total += m
            mean_total.append(total)
        ax.plot(list(means), mean_total, color="#2a78d6", linewidth=2.5, marker="o", label="mean over orders")
        for s, (m, t) in enumerate(zip(means.values(), mean_total), start=1):
            if s > 1:
                ax.annotate(f"+{m:.2f}%", (s, t), textcoords="offset points", xytext=(6, -12), fontsize=8)
        ax.set_xticks(list(means))
        ax.set_xlabel("Libraries combined")
        ax.set_ylabel(f"{args.assembly} masked by any of them (%)")
        ax.legend(frameon=False, loc="lower right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        for path in figures:
            fig.savefig(path, dpi=200)
        print(f"figure: {figures[0]} and .pdf")


if __name__ == "__main__":
    main()
