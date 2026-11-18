#!/usr/bin/env python3
"""Divergence landscape (question 2, briefing C2.6): masked bases per class by divergence from the consensus.

Do young (low-divergence) and old families both look sensible per class? Each annotation's divergence is
RepeatMasker's "perc div." from the .out file, carried in the BED's divergence column: the share of mismatches to
the consensus, NOT the Kimura divergence with CpG adjustment that RepeatMasker's own landscape tools compute from
alignments (our maskings ran without -a, so no alignments exist). Treat this as an approximate landscape; for an
exact one the final masking would add -a.

Simple repeats and low-complexity annotations are left out. Overlapping annotations (nested elements) are each
counted, as landscapes built from alignments also do.

Usage (from the repository root):
  python scripts/divergence_landscape.py --version RI_v5 --assembly DT_004f          # dry run
  python scripts/divergence_landscape.py --version RI_v5 --assembly DT_004f --go

Writes, never overwriting, in <scratch>/comparisons/divergence/: <version>_<assembly>_landscape.tsv (class, bin,
bases, share of the genome) and, if matplotlib is available, <version>_<assembly>_landscape.png.
"""
import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

SIMPLE = {"Simple_repeat", "Low_complexity"}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a version that has masked the assembly, e.g. RI_v5")
    parser.add_argument("--assembly", required=True, help="e.g. DT_004f")
    parser.add_argument("--bin", type=float, default=1.0, help="bin width in percent divergence (default 1)")
    parser.add_argument("--max", type=float, default=50.0, help="last bin collects everything from here (default 50)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def landscape(lines: Iterable[str], width: float, top: float) -> Dict[str, Counter]:
    """class group -> {bin start: bases}, from BED lines (divergence in column 8)."""
    out: Dict[str, Counter] = defaultdict(Counter)
    last = int(top // width)
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if c.top_class(f[6]) in SIMPLE:
            continue
        b = min(int(float(f[7]) // width), last)
        out[la.group(f[6])][round(b * width, 6)] += int(f[2]) - int(f[1])
    return out


def weighted_median(bins: Counter) -> float:
    total, running = sum(bins.values()), 0
    for b in sorted(bins):
        running += bins[b]
        if running >= total / 2:
            return b
    return 0.0


def main(argv=None) -> None:
    args = parse_args(argv)
    bed = c.version_folder(args.version) / "masked" / args.assembly / f"{args.assembly}.bed"
    if not bed.exists():
        c.stop(f"{args.version} has not masked {args.assembly}: {bed} not found")
    out_dir = args.scratch / "comparisons" / "divergence"
    stem = f"{args.version}_{args.assembly}_landscape"
    targets = [out_dir / f"{stem}.tsv", out_dir / f"{stem}.png"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    print(f"{args.assembly} masked by {args.version}: bins of {args.bin:g}% divergence up to {args.max:g}%\nresults: {out_dir}/{stem}.*")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome_bp = sum(c.fasta_lengths(c.genome_path(args.assembly, c.panel_accession(args.assembly))).values())
    with open(bed) as handle:
        land = landscape(handle, args.bin, args.max)
    rows = []
    for g in la.GROUPS:
        for b in sorted(land.get(g, {})):
            rows.append([g, b, land[g][b], f"{100 * land[g][b] / genome_bp:.4f}"])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["class", "divergence_bin_start", "bp", "percent_of_genome"], rows)
    print(f"\n{args.assembly} ({genome_bp:,} bp); divergence = mismatches to the consensus (.out), not Kimura-corrected")
    for g in la.GROUPS:
        bins = land.get(g)
        if not bins:
            continue
        total = sum(bins.values())
        young = sum(n for b, n in bins.items() if b < 5)
        print(f"  {g:8s} {100 * total / genome_bp:5.2f}% of the genome; median divergence {weighted_median(bins):.0f}%; "
              f"under 5%: {100 * young / total:.0f}% of its bases")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"table: {targets[0]} (matplotlib not available: no plot)")
        return
    fig, ax = plt.subplots(figsize=(9, 4.5))
    starts = sorted({b for g in land for b in land[g]})
    bottom = [0.0] * len(starts)
    for g in la.GROUPS:
        heights = [100 * land.get(g, {}).get(b, 0) / genome_bp for b in starts]
        ax.bar(starts, heights, width=args.bin * 0.9, bottom=bottom, align="edge", label=g)
        bottom = [x + y for x, y in zip(bottom, heights)]
    ax.set_xlabel("divergence from the consensus (% mismatches, not Kimura-corrected)")
    ax.set_ylabel("% of the genome")
    ax.set_title(f"{args.assembly} masked by {args.version}: approximate repeat landscape")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(targets[1], dpi=150)
    print(f"table: {targets[0]}\nplot:  {targets[1]}")


if __name__ == "__main__":
    main()
