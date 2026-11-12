#!/usr/bin/env python3
"""Do two libraries give the same bases the same label? (cross-tabulation of classes)

On one genome masked by two or more libraries, for every pair A, B:
  - a class matrix in base pairs: of the bases A labels with a class (rows), how many B
    labels with each class (columns), and how many B leaves unmasked; plus a row for the
    bases B masks that A leaves unmasked;
  - agreement: of the bases both libraries mask, the share given the same class;
  - the family pairs that carry most of the disagreement (same bases, different class),
    largest first: the families to inspect or curate.
Classes are grouped as DNA, LINE, LTR, SINE, Unknown and Other (everything else, including
RepeatMasker's own Simple_repeat and Low_complexity).
Each base gets one label per library: that of the highest-scoring annotation covering it (RepeatMasker's
own precedence for overlapping hits; ties go to the longer annotation). RepeatMasker can give a base two
annotations, for example a simple repeat inside a transposon; counting both would count the base twice.

Inputs: <version folder>/masked/<assembly>/<assembly>.bed (made by mask_genomes.py).
Results: <scratch>/comparisons/label_agreement/<assembly>_<A>_vs_<B>.tsv (class matrix),
..._families.tsv (top family pairs) and, with --plot, ..._matrix.png and .pdf. Never overwritten.

Example:
  python scripts/label_agreement.py --assembly YH_YH2f --versions RI_v1 RI_v3              # dry run
  python scripts/label_agreement.py --assembly YH_YH2f --versions RI_v1 RI_v3 --plot --go
"""
from __future__ import annotations

import argparse
import heapq
import itertools
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

GROUPS = ["DNA", "LINE", "LTR", "SINE", "Unknown", "Other"]
UNMASKED = "unmasked"
Row = Tuple[str, int, int, str, str]   # sequence, start, end, family, full label
Scored = Tuple[str, int, int, str, str, int]   # the same, plus RepeatMasker's alignment (SW) score


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly", required=True, help="a genome masked by every version listed")
    parser.add_argument("--versions", nargs="+", required=True, help="two or more versions, e.g. RI_v1 RI_v3")
    parser.add_argument("--top", type=int, default=50, help="family pairs to list per pair of libraries (default 50)")
    parser.add_argument("--plot", action="store_true", help="also draw the class matrix (needs matplotlib)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def group(label: str) -> str:
    cls = c.top_class(label)
    return cls if cls in GROUPS else "Other"


def read_labelled_bed(path: Path) -> List[Scored]:
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3], f[6], int(float(f[8]))))
    return rows


def one_label_per_base(rows: List[Scored]) -> List[Row]:
    """Non-overlapping rows: each masked base keeps only its highest-scoring annotation (ties: the longer one)."""
    by_seq: Dict[str, List[Scored]] = defaultdict(list)
    for r in rows:
        by_seq[r[0]].append(r)
    out: List[Row] = []
    for name, items in by_seq.items():
        items.sort(key=lambda r: r[1])
        points = sorted({p for r in items for p in (r[1], r[2])})
        heap: list = []
        i = 0
        open_seg = None   # [start, end, family, label] of the segment being extended
        for k in range(len(points) - 1):
            x, nxt = points[k], points[k + 1]
            while i < len(items) and items[i][1] <= x:
                _n, start, end, family, label, score = items[i]
                heapq.heappush(heap, (-score, -(end - start), family, end, label))
                i += 1
            while heap and heap[0][3] <= x:
                heapq.heappop(heap)
            if not heap:
                if open_seg:
                    out.append((name, *open_seg))
                    open_seg = None
                continue
            family, label = heap[0][2], heap[0][4]
            if open_seg and open_seg[1] == x and open_seg[2] == family and open_seg[3] == label:
                open_seg[1] = nxt
            else:
                if open_seg:
                    out.append((name, *open_seg))
                open_seg = [x, nxt, family, label]
        if open_seg:
            out.append((name, *open_seg))
    return [tuple(r) for r in out]


def group_sets(rows: List[Row]):
    by_group: Dict[str, list] = defaultdict(list)
    for name, start, end, _family, label in rows:
        by_group[group(label)].append((name, start, end))
    everything = c.merge_by_sequence((r[0], r[1], r[2]) for r in rows)
    return everything, {g: c.merge_by_sequence(spans) for g, spans in by_group.items()}


def class_matrix(a_rows: List[Row], b_rows: List[Row]) -> Dict[Tuple[str, str], int]:
    """(group in A or 'unmasked', group in B or 'unmasked') -> base pairs."""
    a_all, a_groups = group_sets(a_rows)
    b_all, b_groups = group_sets(b_rows)
    matrix: Dict[Tuple[str, str], int] = {}
    for ga, a_spans in a_groups.items():
        for gb, b_spans in b_groups.items():
            matrix[(ga, gb)] = c.intersect_length(a_spans, b_spans)
        matrix[(ga, UNMASKED)] = c.subtract_length(a_spans, b_all)
    for gb, b_spans in b_groups.items():
        matrix[(UNMASKED, gb)] = c.subtract_length(b_spans, a_all)
    return matrix


def agreement(matrix: Dict[Tuple[str, str], int]) -> Tuple[int, int]:
    """(bases masked by both, of which given the same group), counted from the matrix."""
    both = sum(v for (ga, gb), v in matrix.items() if UNMASKED not in (ga, gb))
    same = sum(v for (ga, gb), v in matrix.items() if ga == gb and ga != UNMASKED)
    return both, same


def discordant_families(a_rows: List[Row], b_rows: List[Row]) -> Dict[Tuple[str, str, str, str], int]:
    """(family in A, label in A, family in B, label in B) -> overlapping base pairs where the groups differ."""
    by_seq_a: Dict[str, List[Row]] = defaultdict(list)
    by_seq_b: Dict[str, List[Row]] = defaultdict(list)
    for r in a_rows:
        by_seq_a[r[0]].append(r)
    for r in b_rows:
        by_seq_b[r[0]].append(r)
    pairs: Dict[Tuple[str, str, str, str], int] = defaultdict(int)
    for name, a_list in by_seq_a.items():
        b_list = sorted(by_seq_b.get(name, []), key=lambda r: r[1])
        if not b_list:
            continue
        a_list.sort(key=lambda r: r[1])
        j = 0
        for _n, a_start, a_end, a_fam, a_label in a_list:
            while j < len(b_list) and b_list[j][2] <= a_start:
                j += 1
            k = j
            while k < len(b_list) and b_list[k][1] < a_end:
                _m, b_start, b_end, b_fam, b_label = b_list[k]
                overlap = min(a_end, b_end) - max(a_start, b_start)
                if overlap > 0 and group(a_label) != group(b_label):
                    pairs[(a_fam, a_label, b_fam, b_label)] += overlap
                k += 1
    return pairs


def main(argv=None) -> None:
    args = parse_args(argv)
    if len(set(args.versions)) < 2 or len(set(args.versions)) != len(args.versions):
        c.stop("give two or more different versions")
    beds = {}
    for version in args.versions:
        bed = c.version_folder(version) / "masked" / args.assembly / f"{args.assembly}.bed"
        if not bed.exists():
            c.stop(f"{version} has not masked {args.assembly}: {bed} not found")
        beds[version] = bed
        print(f"{version}: {bed}")
    out_dir = args.scratch / "comparisons" / "label_agreement"
    pairs = list(itertools.combinations(args.versions, 2))
    targets = []
    for a, b in pairs:
        stem = out_dir / f"{args.assembly}_{a}_vs_{b}"
        targets += [Path(f"{stem}.tsv"), Path(f"{stem}_families.tsv")]
        if args.plot:
            targets += [Path(f"{stem}_matrix.png"), Path(f"{stem}_matrix.pdf")]
    print(f"pairs: {', '.join(f'{a} vs {b}' for a, b in pairs)}; results in {out_dir}")
    existing = [p for p in targets if p.exists()]
    if existing:
        c.stop("already exists, never overwritten: " + ", ".join(map(str, existing)))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    genome_bp = sum(c.fasta_lengths(genome).values())
    rows = {}
    for version, bed in beds.items():
        print(f"reading {version} ...", flush=True)
        raw = read_labelled_bed(bed)
        rows[version] = one_label_per_base(raw)
        print(f"  {len(raw):,} annotations -> {len(rows[version]):,} stretches with one label each", flush=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    for a, b in pairs:
        stem = out_dir / f"{args.assembly}_{a}_vs_{b}"
        matrix = class_matrix(rows[a], rows[b])
        row_totals: Dict[str, int] = defaultdict(int)
        for (ga, _gb), v in matrix.items():
            row_totals[ga] += v
        order = GROUPS + [UNMASKED]
        table = []
        for ga in order:
            for gb in order:
                if (ga, gb) in matrix:
                    v = matrix[(ga, gb)]
                    share = 100 * v / row_totals[ga] if row_totals[ga] else 0.0
                    table.append([args.assembly, f"{a} vs {b}", ga, gb, v, f"{100 * v / genome_bp:.3f}", f"{share:.1f}"])
        c.write_tsv(Path(f"{stem}.tsv"), ["assembly", "pair", f"class_in_{a}", f"class_in_{b}", "bp",
                                          "percent_of_genome", "percent_of_row"], table)
        families = discordant_families(rows[a], rows[b])
        top = sorted(families.items(), key=lambda kv: -kv[1])[: args.top]
        c.write_tsv(Path(f"{stem}_families.tsv"),
                    ["assembly", "pair", f"family_in_{a}", f"label_in_{a}", f"family_in_{b}", f"label_in_{b}",
                     "bp", "percent_of_genome"],
                    [[args.assembly, f"{a} vs {b}", fa, la, fb, lb, v, f"{100 * v / genome_bp:.4f}"]
                     for (fa, la, fb, lb), v in top])
        both, same = agreement(matrix)
        direct = c.intersect_length(c.merge_by_sequence((r[0], r[1], r[2]) for r in rows[a]),
                                    c.merge_by_sequence((r[0], r[1], r[2]) for r in rows[b]))
        if both != direct:
            c.stop(f"internal check failed: the class matrix holds {both} bp masked by both, the direct count {direct}")
        print(f"\n{args.assembly}, {a} vs {b}: masked by both {100 * both / genome_bp:.2f}% of the genome; "
              f"same class for {100 * same / both:.1f}% of those bases")
        for ga in GROUPS:
            if row_totals.get(ga):
                parts = [f"{gb} {100 * matrix.get((ga, gb), 0) / row_totals[ga]:.1f}%" for gb in order
                         if matrix.get((ga, gb), 0) / row_totals[ga] >= 0.005]
                print(f"  of {a}'s {ga} ({100 * row_totals[ga] / genome_bp:.2f}% of the genome), {b} calls: "
                      + ", ".join(parts))
        print(f"  largest disagreements (same bases, different class), of {len(families)} family pairs:")
        for (fa, la, fb, lb), v in top[:10]:
            print(f"    {fa} ({la})  vs  {fb} ({lb}): {v:,} bp")
        print(f"  check passed: every base masked by both is counted once ({both:,} bp)")
        print(f"  written: {stem}.tsv and {stem}_families.tsv")

        if args.plot:
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
            except ImportError:
                c.stop("matplotlib is not installed in this environment; the tables were written, the figure was not")
            rows_shown = [g for g in GROUPS if row_totals.get(g)]
            cols_shown = [g for g in order if any((r, g) in matrix for r in rows_shown)]
            grid = [[100 * matrix.get((r, g), 0) / row_totals[r] for g in cols_shown] for r in rows_shown]
            fig, ax = plt.subplots(figsize=(1.1 * len(cols_shown) + 2.4, 0.7 * len(rows_shown) + 1.8))
            image = ax.imshow(grid, cmap="Greys", vmin=0, vmax=100)
            for i, r in enumerate(rows_shown):
                for j, g in enumerate(cols_shown):
                    value = grid[i][j]
                    ax.text(j, i, f"{value:.0f}" if value >= 0.5 else "", ha="center", va="center", fontsize=9,
                            color="white" if value > 55 else "black")
            ax.set_xticks(range(len(cols_shown)))
            ax.set_xticklabels(cols_shown, fontsize=9)
            ax.set_yticks(range(len(rows_shown)))
            ax.set_yticklabels([f"{r} ({100 * row_totals[r] / genome_bp:.1f}%)" for r in rows_shown], fontsize=9)
            ax.set_xlabel(f"label given by {b}")
            ax.set_ylabel(f"label given by {a} (share of the genome)")
            ax.set_title(f"{args.assembly}: of the bases {a} labels with each class, what {b} calls them (%)",
                         fontsize=9)
            fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(Path(f"{stem}_matrix.{ext}"), dpi=200)
            plt.close(fig)
            print(f"  figure: {stem}_matrix.png and .pdf")


if __name__ == "__main__":
    main()
