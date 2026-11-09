#!/usr/bin/env python3
"""Figures for the question 1 pilot (needs matplotlib in ga-repeats).

From one or more comparison folders (compare_libraries.py) and, when present, their
masking summaries (summarize_masking.py) and the masking-overlap tables
(compare_masking.py), writes PNG and PDF figures:

  fig1_shared_clusters_by_threshold   clusters holding both libraries, per identity threshold
  fig2_best_matches                   each family's best match in the other library
  fig3_shared_clusters_by_class       clusters holding both libraries, per class, at --threshold
  fig4_class_composition              percent of each genome masked per class, per library
                                      (genome-wide, after Quah et al. Fig. 7A)
  fig5_masking_overlap                bases masked by both libraries or by only one

Figures go to <scratch>/comparisons/plots/ (or --out). Existing files are never replaced.

Example:
  python scripts/plot_pilot.py --comparisons ~/scratch/ga-repeats/comparisons/RI_v1_vs_RI_v2 \\
                                             ~/scratch/ga-repeats/comparisons/RI_v1_vs_RI_v3 --go
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

COMPARISON_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4"]
CLASS_COLOURS = OrderedDict([("DNA", "#2a78d6"), ("LINE", "#1baf7a"), ("LTR", "#eb6834"),
                             ("SINE", "#eda100"), ("Unknown", "#8c8a84"), ("Other", "#cfcdc4")])
BAR_CLASSES = ["LTR", "DNA", "LINE", "SINE", "Unknown", "ALL"]


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--comparisons", nargs="+", type=Path, required=True,
                        help="comparison folders made by compare_libraries.py")
    parser.add_argument("--threshold", type=float, default=0.80, help="threshold for figures 3 and 4 (default 0.80)")
    parser.add_argument("--overlap-dir", type=Path, default=None,
                        help="folder of compare_masking.py tables (default <scratch>/comparisons/masking_overlap)")
    parser.add_argument("--out", type=Path, default=None, help="figure folder (default <scratch>/comparisons/plots)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as handle:
        return [{k.strip(): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def load_comparison(folder: Path, threshold: float) -> dict:
    meta = json.loads((folder / "meta.json").read_text())
    summary = folder / f"masking_summary_c{threshold:.2f}.tsv"
    return {"label": f"{meta['id_a']} vs {meta['id_b']}", "meta": meta,
            "shared": read_tsv(folder / "shared_unique_by_class.tsv"),
            "similarity": read_tsv(folder / "similarity_by_class.tsv"),
            "masking": read_tsv(summary) if summary.exists() else []}


def shared_percent(rows: List[Dict[str, str]], threshold: str, cls: str):
    for r in rows:
        if r["threshold"] == threshold and r["class"] == cls:
            s, a, b = int(r["clusters_shared"]), int(r["clusters_only_a"]), int(r["clusters_only_b"])
            return 100 * s / (s + a + b) if s + a + b else None
    return None


def class_composition(comparisons: List[dict]) -> "OrderedDict[str, Dict[str, float]]":
    """(library on genome) -> percent masked per class group. A masking that appears in two
    comparison folders (RI_v1 is in both) is counted once."""
    seen = set()
    out: "OrderedDict[str, Dict[str, float]]" = OrderedDict()
    for comp in comparisons:
        for r in comp["masking"]:
            key = f"{r['library_version']} ({r['library_built_from']}) on {r['assembly']}"
            if (key, r["class"]) in seen:
                continue
            seen.add((key, r["class"]))
            group = out.setdefault(key, {k: 0.0 for k in CLASS_COLOURS})
            if r["class"] == "ALL":
                group["_all"] = float(r["percent_masked"])
                continue
            group[r["class"] if r["class"] in CLASS_COLOURS else "Other"] += float(r["percent_masked"])
    return out


def overlap_bars(overlap_dir: Path) -> "OrderedDict[str, tuple]":
    bars: "OrderedDict[str, tuple]" = OrderedDict()
    if not overlap_dir or not overlap_dir.exists():
        return bars
    for table in sorted(overlap_dir.glob("*.tsv")):
        rows = read_tsv(table)
        for pair in OrderedDict.fromkeys(r["pair"] for r in rows):
            a, b = pair.split(" vs ")
            ra = next(r for r in rows if r["pair"] == pair and r["library"] == a and r["class"] == "ALL")
            rb = next(r for r in rows if r["pair"] == pair and r["library"] == b and r["class"] == "ALL")
            genome = int(ra["genome_bp"])
            only_a, only_b = int(ra["masked_bp_not_masked_by_other"]), int(rb["masked_bp_not_masked_by_other"])
            both = int(ra["masked_bp"]) - only_a
            bars[f"{ra['assembly']}: {pair}"] = (100 * both / genome, 100 * only_a / genome,
                                                 100 * only_b / genome, a, b)
    return bars


def main(argv=None) -> None:
    args = parse_args(argv)
    out = args.out or (args.scratch / "comparisons" / "plots")
    overlap_dir = args.overlap_dir or (args.scratch / "comparisons" / "masking_overlap")
    comparisons = [load_comparison(f, args.threshold) for f in args.comparisons]
    t = f"{args.threshold:.2f}"
    names = ["fig1_shared_clusters_by_threshold", "fig2_best_matches", "fig3_shared_clusters_by_class"]
    composition = class_composition(comparisons)
    if composition:
        names.append("fig4_class_composition")
    overlaps = overlap_bars(overlap_dir)
    if overlaps:
        names.append("fig5_masking_overlap")
    print("comparisons: " + "; ".join(comp["label"] for comp in comparisons))
    print(f"masking summaries at {t}: {sum(1 for comp in comparisons if comp['masking'])} of {len(comparisons)}; "
          f"overlap tables: {len(overlaps)} pairs")
    print(f"figures to {out}: " + ", ".join(names) + " (.png and .pdf)")
    clashes = [out / f"{n}.{ext}" for n in names for ext in ("png", "pdf") if (out / f"{n}.{ext}").exists()]
    if clashes:
        c.stop("these figures already exist and are never replaced: " + ", ".join(map(str, clashes)))
    if not args.go:
        print("dry run: nothing was written. Add --go to draw.")
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        c.stop("matplotlib is not installed in this environment; see INSTRUCTIONS.md of repeats_011")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    def save(fig, name):
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(out / f"{name}.{ext}", dpi=200)
        plt.close(fig)
        print(f"written: {out / name}.png and .pdf")

    # 1: shared clusters by threshold
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, comp in enumerate(comparisons):
        thresholds = sorted({r["threshold"] for r in comp["shared"]})
        ys = [shared_percent(comp["shared"], th, "ALL") for th in thresholds]
        ax.plot([float(th) for th in thresholds], ys, marker="os^Dv"[i % 5], linestyle=["-", "--", ":", "-."][i % 4],
                color=COMPARISON_COLOURS[i % 5], linewidth=2, label=comp["label"])
    ax.set_xlabel("cd-hit-est identity threshold")
    ax.set_ylabel("Clusters holding both libraries (%)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    save(fig, "fig1_shared_clusters_by_threshold")

    # 2: best matches
    labels, covering, partial, none = [], [], [], []
    coverage_label = "match covers the coverage threshold"
    for comp in comparisons:
        meta = comp["meta"]
        names = {meta["a"]: f"{meta['a']} ({meta['id_a']})", meta["b"]: f"{meta['b']} ({meta['id_b']})"}
        for r in comp["similarity"]:
            if r["class"] != "ALL":
                continue
            cover_col = next(k for k in r if k.startswith("match_covers_"))
            coverage_label = f"match covers at least {100 * float(cover_col.rsplit('_', 1)[1]):.0f}% of the family"
            n, m, cv = int(r["families"]), int(r["with_a_match"]), int(r[cover_col])
            query, _, library = r["direction"].partition("_in_")
            labels.append(f"{names.get(query, query)} in {names.get(library, library)}")
            covering.append(100 * cv / n)
            partial.append(100 * (m - cv) / n)
            none.append(100 * (n - m) / n)
    fig, ax = plt.subplots(figsize=(9, 0.6 * len(labels) + 1.8))
    y = range(len(labels))
    ax.barh(y, covering, color="#2a78d6", label=coverage_label)
    ax.barh(y, partial, left=covering, color="#1baf7a", label="partial match only")
    ax.barh(y, none, left=[a + b for a, b in zip(covering, partial)], color="#b4b2a9", label="no match")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Families (%)")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=8)
    save(fig, "fig2_best_matches")

    # 3: shared clusters by class
    fig, ax = plt.subplots(figsize=(7, 4))
    width = 0.8 / max(1, len(comparisons))
    for i, comp in enumerate(comparisons):
        values = [shared_percent(comp["shared"], t, cls) or 0 for cls in BAR_CLASSES]
        ax.bar([x + i * width for x in range(len(BAR_CLASSES))], values, width=width,
               color=COMPARISON_COLOURS[i % 5], hatch=["", "//", "..", "xx", "\\\\"][i % 5],
               edgecolor="white", label=comp["label"])
    ax.set_xticks([x + width * (len(comparisons) - 1) / 2 for x in range(len(BAR_CLASSES))])
    ax.set_xticklabels(["All" if cls == "ALL" else cls for cls in BAR_CLASSES])
    ax.set_ylabel(f"Clusters holding both libraries at {t} (%)")
    ax.legend(frameon=False)
    save(fig, "fig3_shared_clusters_by_class")

    # 4: class composition per library and genome
    if composition:
        keys = list(composition)
        fig, ax = plt.subplots(figsize=(8, 0.55 * len(keys) + 1.6))
        left = [0.0] * len(keys)
        for cls, colour in CLASS_COLOURS.items():
            values = [composition[k][cls] for k in keys]
            ax.barh(range(len(keys)), values, left=left, color=colour, label=cls, edgecolor="white")
            left = [a + b for a, b in zip(left, values)]
        for i, k in enumerate(keys):
            if "_all" in composition[k]:
                ax.text(left[i] + 0.4, i, f"{composition[k]['_all']:.1f}% masked", va="center", fontsize=8)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(keys, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Genome masked, by class (%); classes can overlap, so a bar can run past its total")
        ax.set_xlim(0, max(left) + 12)
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=6, fontsize=8)
        save(fig, "fig4_class_composition")

    # 5: masking overlap
    if overlaps:
        keys = list(overlaps)
        fig, ax = plt.subplots(figsize=(8, 0.55 * len(keys) + 1.6))
        both = [overlaps[k][0] for k in keys]
        only_a = [overlaps[k][1] for k in keys]
        only_b = [overlaps[k][2] for k in keys]
        ax.barh(range(len(keys)), both, color="#8c8a84", label="masked by both")
        ax.barh(range(len(keys)), only_a, left=both, color="#2a78d6", label="only by the first library")
        ax.barh(range(len(keys)), only_b, left=[a + b for a, b in zip(both, only_a)], color="#eb6834",
                label="only by the second library")
        for i, k in enumerate(keys):
            ax.text(both[i] + only_a[i] + only_b[i] + 0.4, i,
                    f"{only_a[i]:.2f}% / {only_b[i]:.2f}%", va="center", fontsize=8)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(keys, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Genome (%)")
        ax.set_xlim(0, max(a + b + d for a, b, d in zip(both, only_a, only_b)) + 10)
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8)
        save(fig, "fig5_masking_overlap")


if __name__ == "__main__":
    main()
