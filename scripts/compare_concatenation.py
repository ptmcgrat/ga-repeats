#!/usr/bin/env python3
"""Compare a concatenated library's masking of one genome with a merged library's (the concatenation test).

Decided 2026-11-17: before choosing between merging the libraries (merge_libraries.py, clustering) and
concatenating them (merge_libraries.py --no-clustering, RepeatMasker's scores choosing among near-duplicate
families), mask one genome with both, built from the same sources, and compare:

  coverage        bases masked by each, and by one but not the other (one label per base: the
                  highest-scoring annotation)
  class shares    share of the genome per class group (DNA, LINE, LTR, SINE, Unknown, Other)
  labels          on bases masked by both, the share given the same class group; the largest disagreements
  fragmentation   interspersed repeats only (Simple_repeat and Low_complexity left out): annotations and
                  RepeatMasker elements (fragments sharing an element ID) per Mb of genome, the median
                  annotation length, and annotations per element
  family split    for every cluster of the merged library with two or more members: how the concatenated
                  masking divides the cluster's bases among its members (the top member's share; members
                  holding at least 10%)
  run time        from each masking's run record (minute resolution; the server's load is not recorded,
                  so small differences mean little)

Usage (from the repository root):
  python scripts/compare_concatenation.py --concatenated RI_v8 --merged RI_v5 --assembly DT_004f        # dry run
  python scripts/compare_concatenation.py --concatenated RI_v8 --merged RI_v5 --assembly DT_004f --go

Writes, never overwriting, in <scratch>/comparisons/concatenation/:
  <concatenated>_vs_<merged>_<assembly>_summary.tsv    one row per measure and version
  <concatenated>_vs_<merged>_<assembly>_clusters.tsv   one row per merged cluster with 2+ members
"""
import argparse
import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

SIMPLE = {"Simple_repeat", "Low_complexity"}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--concatenated", required=True, help="a version written by merge_libraries.py --no-clustering, e.g. RI_v8")
    parser.add_argument("--merged", required=True, help="a merged version of the same sources, e.g. RI_v5")
    parser.add_argument("--assembly", required=True, help="a genome masked by both, e.g. DT_004f")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def fragmentation(bed: Path) -> Dict[str, float]:
    """Interspersed annotations (BED rows), elements (distinct sequence and element id) and median annotation length."""
    lengths, elements = [], set()
    with open(bed) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if c.top_class(f[6]) in SIMPLE:
                continue
            lengths.append(int(f[2]) - int(f[1]))
            elements.add((f[0], f[9]))
    return {"annotations": len(lengths), "elements": len(elements),
            "median_length": statistics.median(lengths) if lengths else 0}


def cluster_members(table: Path) -> Dict[str, str]:
    """Family name (without its #label) -> cluster, for the clusters with two or more members."""
    with open(table) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    size = Counter(r["cluster"] for r in rows)
    return {r["family"].split("#")[0]: r["cluster"] for r in rows if size[r["cluster"]] >= 2}


def family_split(rows: List[la.Row], member_of: Dict[str, str]) -> Dict[str, Counter]:
    """cluster -> {member family: bases it holds in the (one label per base) masking}."""
    out: Dict[str, Counter] = defaultdict(Counter)
    for _seq, start, end, family, _label in rows:
        cluster = member_of.get(family)
        if cluster is not None:
            out[cluster][family] += end - start
    return out


def split_summary(split: Dict[str, Counter]) -> Dict[str, float]:
    totals = {cl: sum(v.values()) for cl, v in split.items()}
    total = sum(totals.values())
    under_half = [cl for cl, v in split.items() if max(v.values()) < 0.5 * totals[cl]]
    shared = [cl for cl, v in split.items() if sum(1 for x in v.values() if x >= 0.10 * totals[cl]) >= 2]
    return {"clusters": len(split), "bp": total,
            "top_member_share": sum(max(v.values()) for v in split.values()) / total if total else 0.0,
            "clusters_no_member_half": len(under_half), "bp_no_member_half": sum(totals[cl] for cl in under_half),
            "clusters_two_members_10pct": len(shared)}


def run_minutes(scratch: Path, version: str, assembly: str) -> Optional[int]:
    record = scratch / version / "masking" / assembly / "run_record.json"
    if not record.exists():
        return None
    r = json.loads(record.read_text())
    started = datetime.strptime(r["started"], "%Y-%m-%d %H:%M")
    finished = datetime.strptime(r["finished"], "%Y-%m-%d %H:%M")
    return int((finished - started).total_seconds() // 60)


def sources_of(version: str) -> List[str]:
    record = c.version_folder(version) / "library" / "merge_record.json"
    if not record.exists():
        c.stop(f"{version} has no merge_record.json: it is not a merged or concatenated version")
    return sorted(s["version"] for s in json.loads(record.read_text())["sources"])


def main(argv=None) -> None:
    args = parse_args(argv)
    conc, merged = args.concatenated, args.merged
    record = json.loads((c.version_folder(conc) / "library" / "merge_record.json").read_text()) \
        if (c.version_folder(conc) / "library" / "merge_record.json").exists() else {}
    if record.get("kind") != "concatenated library":
        c.stop(f"{conc} is not a concatenated library (merge_libraries.py --no-clustering)")
    if sources_of(conc) != sources_of(merged):
        c.stop(f"{conc} and {merged} were built from different sources: {sources_of(conc)} and {sources_of(merged)}")
    table = c.version_folder(merged) / "library" / "cluster_table.tsv"
    if not table.exists():
        c.stop(f"{merged} has no cluster table ({table}): it is not a merged library")
    beds = {}
    for v in (conc, merged):
        beds[v] = c.version_folder(v) / "masked" / args.assembly / f"{args.assembly}.bed"
        if not beds[v].exists():
            c.stop(f"{v} has not masked {args.assembly}: {beds[v]} not found")
    out_dir = args.scratch / "comparisons" / "concatenation"
    stem = f"{conc}_vs_{merged}_{args.assembly}"
    targets = [out_dir / f"{stem}_summary.tsv", out_dir / f"{stem}_clusters.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    print(f"{args.assembly}: {conc} (concatenated, sources {', '.join(sources_of(conc))}) against {merged} (merged)")
    print(f"results: {out_dir}/{stem}_*.tsv")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome_bp = sum(c.fasta_lengths(c.genome_path(args.assembly, c.panel_accession(args.assembly))).values())
    print("reading the two maskings ...", flush=True)
    rows = {v: la.one_label_per_base(la.read_labelled_bed(beds[v])) for v in (conc, merged)}
    summary = []
    pct = lambda bp: f"{100 * bp / genome_bp:.3f}"  # noqa: E731
    shares = {}
    for v in (conc, merged):
        groups: Counter = Counter()
        for _s, start, end, _f, label in rows[v]:
            groups[la.group(label)] += end - start
        shares[v] = groups
        summary.append([v, "masked_percent", pct(sum(groups.values()))])
        summary += [[v, f"{g}_percent", pct(groups[g])] for g in la.GROUPS]
    matrix = la.class_matrix(rows[conc], rows[merged])
    both, same = la.agreement(matrix)
    only_conc = sum(n for (a, b), n in matrix.items() if b == la.UNMASKED)
    only_merged = sum(n for (a, b), n in matrix.items() if a == la.UNMASKED)
    pairs = sorted(((n, a, b) for (a, b), n in matrix.items() if la.UNMASKED not in (a, b) and a != b), reverse=True)
    summary += [["both", "masked_by_both_percent", pct(both)], ["both", "same_class_group_percent_of_both", f"{100 * same / both:.2f}" if both else ""],
                [conc, "masked_only_by_this_percent", pct(only_conc)], [merged, "masked_only_by_this_percent", pct(only_merged)]]
    frag = {v: fragmentation(beds[v]) for v in (conc, merged)}
    mb = genome_bp / 1e6
    for v in (conc, merged):
        f = frag[v]
        summary += [[v, "interspersed_annotations_per_Mb", f"{f['annotations'] / mb:.1f}"],
                    [v, "interspersed_elements_per_Mb", f"{f['elements'] / mb:.1f}"],
                    [v, "median_annotation_length", f"{f['median_length']:.0f}"],
                    [v, "annotations_per_element", f"{f['annotations'] / f['elements']:.3f}" if f["elements"] else ""]]
    split = family_split(rows[conc], cluster_members(table))
    ss = split_summary(split)
    summary += [[conc, f"split_{k}", f"{v:.4f}" if isinstance(v, float) else v] for k, v in ss.items()]
    minutes = {v: run_minutes(args.scratch, v, args.assembly) for v in (conc, merged)}
    summary += [[v, "run_minutes", m if m is not None else ""] for v, m in minutes.items()]

    cluster_rows = []
    for cl, v in sorted(split.items(), key=lambda kv: -sum(kv[1].values())):
        total = sum(v.values())
        top, top_bp = max(v.items(), key=lambda kv: kv[1])
        cluster_rows.append([cl, len(v), total, top, f"{top_bp / total:.3f}", sum(1 for x in v.values() if x >= 0.10 * total)])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["version", "measure", "value"], summary)
    c.write_tsv(targets[1], ["cluster", "members_masking", "bp", "top_member", "top_member_share", "members_with_10pct"], cluster_rows)

    print(f"\n{args.assembly} ({genome_bp:,} bp): {conc} (concatenated) against {merged} (merged)")
    print(f"coverage: {conc} {100 * sum(shares[conc].values()) / genome_bp:.2f}%, {merged} {100 * sum(shares[merged].values()) / genome_bp:.2f}%; "
          f"masked only by {conc} {100 * only_conc / genome_bp:.2f}%, only by {merged} {100 * only_merged / genome_bp:.2f}%")
    print(f"class shares ({conc} / {merged}): " + ", ".join(
        f"{g} {100 * shares[conc][g] / genome_bp:.2f}/{100 * shares[merged][g] / genome_bp:.2f}%" for g in la.GROUPS))
    print(f"labels: of the bases masked by both, {100 * same / both:.1f}% get the same class group; largest disagreements "
          f"({conc} vs {merged}): " + "; ".join(f"{a} vs {b} {100 * n / genome_bp:.2f}%" for n, a, b in pairs[:5]))
    print("fragmentation (interspersed repeats): " + "; ".join(
        f"{v}: {frag[v]['annotations'] / mb:.0f} annotations and {frag[v]['elements'] / mb:.0f} elements per Mb, "
        f"median length {frag[v]['median_length']:.0f} bp, {frag[v]['annotations'] / max(1, frag[v]['elements']):.2f} annotations per element"
        for v in (conc, merged)))
    print(f"family split ({merged}'s {ss['clusters']:,} clusters with 2+ members, {ss['bp'] / 1e6:.1f} Mb in {conc}'s masking): "
          f"the top member holds {100 * ss['top_member_share']:.1f}% of their bases; in {ss['clusters_no_member_half']:,} clusters "
          f"({ss['bp_no_member_half'] / 1e6:.1f} Mb) no member holds half; {ss['clusters_two_members_10pct']:,} clusters have "
          f"2+ members with 10% or more")
    fmt = lambda m: f"{m // 60}:{m % 60:02d}" if m is not None else "no run record"  # noqa: E731
    print(f"run time: {conc} {fmt(minutes[conc])}, {merged} {fmt(minutes[merged])} (server load not recorded)")
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
