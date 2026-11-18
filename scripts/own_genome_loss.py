#!/usr/bin/env python3
"""Lineage-specific loss: what a merged library misses on a build genome, against that genome's own library.

A family found in one lineage only stays its own cluster in a merge, but a lineage-specific subfamily that is
close enough to another library's family (80% identity over 80% of both) joins that family's cluster, and its
copies are then searched with the cluster's representative. This measures what that costs on the genome the
family came from: mask a build genome with the merged library (mask_genomes.py) and compare with its masking
by its own library (one label per base: the highest-scoring annotation).

Reports: masked share by each; bases masked by the own library but not the merged one (the loss) and the reverse,
by class; on bases masked by both, the share given the same class group; and the own library's families behind
the loss, each with its cluster in the merged library (its role there, the cluster's size, representative and label).

Usage (from the repository root):
  python scripts/own_genome_loss.py --version RI_v5 --assembly YH_YH2f          # dry run
  python scripts/own_genome_loss.py --version RI_v5 --assembly YH_YH2f --go

Writes, never overwriting, <scratch>/comparisons/own_genome/<version>_<assembly>_lost_families.tsv.
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a merged version, e.g. RI_v5")
    parser.add_argument("--assembly", required=True, help="the build genome of one of its sources, e.g. YH_YH2f")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def lost_by_family(own: List[la.Row], merged: List[la.Row]) -> Counter:
    """(family, label) -> bases of the own library's (one label per base) stretches that the merged masking leaves unmasked."""
    cover = c.merge_by_sequence((s, a, b) for s, a, b, _f, _l in merged)
    by_seq: Dict[str, list] = defaultdict(list)
    for s, a, b, family, label in own:
        by_seq[s].append((a, b, family, label))
    lost: Counter = Counter()
    for seq, items in by_seq.items():
        spans, j = cover.get(seq, []), 0
        for a, b, family, label in sorted(items):
            while j < len(spans) and spans[j][1] <= a:
                j += 1
            covered, k = 0, j
            while k < len(spans) and spans[k][0] < b:
                covered += min(b, spans[k][1]) - max(a, spans[k][0])
                k += 1
            if b - a - covered:
                lost[(family, label)] += b - a - covered
    return lost


def cluster_info(table: Path) -> Dict[str, dict]:
    """Family (without #label) -> its cluster, role, cluster size, representative and cluster label."""
    with open(table) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    size = Counter(r["cluster"] for r in rows)
    rep = {r["cluster"]: r["family"].split("#")[0] for r in rows if r["representative"] == "yes"}
    return {r["family"].split("#")[0]: {"cluster": r["cluster"], "role": "representative" if r["representative"] == "yes" else "member",
                                        "size": size[r["cluster"]], "representative": rep.get(r["cluster"], ""),
                                        "cluster_label": r["cluster_label"]} for r in rows}


def main(argv=None) -> None:
    args = parse_args(argv)
    lib = c.version_folder(args.version) / "library"
    record_file, table = lib / "merge_record.json", lib / "cluster_table.tsv"
    if not record_file.exists() or not table.exists():
        c.stop(f"{args.version} is not a merged library (needs merge_record.json and cluster_table.tsv in {lib})")
    record = json.loads(record_file.read_text())
    source = next((s["version"] for s in record["sources"] if s.get("build_assembly") == args.assembly), None)
    if source is None:
        c.stop(f"{args.assembly} is not the build genome of any source of {args.version}: "
               + ", ".join(f"{s['version']} ({s.get('build_assembly')})" for s in record["sources"]))
    beds = {v: c.version_folder(v) / "masked" / args.assembly / f"{args.assembly}.bed" for v in (args.version, source)}
    for v, bed in beds.items():
        if not bed.exists():
            c.stop(f"{v} has not masked {args.assembly}: {bed} not found")
    out_dir = args.scratch / "comparisons" / "own_genome"
    target = out_dir / f"{args.version}_{args.assembly}_lost_families.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    print(f"{args.assembly}: {args.version} (merged) against {source} (its own library)\nresult: {target}")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome_bp = sum(c.fasta_lengths(c.genome_path(args.assembly, c.panel_accession(args.assembly))).values())
    print("reading the two maskings ...", flush=True)
    own = la.one_label_per_base(la.read_labelled_bed(beds[source]))
    merged = la.one_label_per_base(la.read_labelled_bed(beds[args.version]))
    matrix = la.class_matrix(own, merged)
    both, same = la.agreement(matrix)
    only_own = Counter({a: n for (a, b), n in matrix.items() if b == la.UNMASKED})
    only_merged = Counter({b: n for (a, b), n in matrix.items() if a == la.UNMASKED})
    pct = lambda n: f"{100 * n / genome_bp:.3f}%"  # noqa: E731
    own_bp, merged_bp = sum(b - a for _, a, b, _, _ in own), sum(b - a for _, a, b, _, _ in merged)
    print(f"\n{args.assembly} ({genome_bp:,} bp): own library {source} masks {pct(own_bp)}, merged {args.version} {pct(merged_bp)}")
    print(f"lost (own library masks, merged does not): {pct(sum(only_own.values()))} -- "
          + ", ".join(f"{g} {pct(n)}" for g, n in only_own.most_common()))
    print(f"gained (merged masks, own library does not): {pct(sum(only_merged.values()))} -- "
          + ", ".join(f"{g} {pct(n)}" for g, n in only_merged.most_common()))
    print(f"labels: of the bases masked by both, {100 * same / both:.1f}% get the same class group")

    lost = lost_by_family(own, merged)
    info = cluster_info(table)
    rows, by_role = [], Counter()
    for (family, label), n in lost.most_common():
        i = info.get(family, {})
        by_role[i.get("role", "not in the cluster table")] += n
        rows.append([family, label, la.group(label), n, i.get("cluster", ""), i.get("role", ""), i.get("size", ""),
                     i.get("representative", ""), i.get("cluster_label", "")])
    total = sum(lost.values())
    if total:
        print("the lost bases by the role of their family in the merge: "
              + ", ".join(f"{k} {100 * v / total:.0f}%" for k, v in by_role.most_common()))
        print("largest: " + "; ".join(f"{r[0]} ({r[1]}, {r[5]} of a cluster of {r[6]}) {r[3] / 1e3:.0f} kb" for r in rows[:6]))
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["family", "label", "class", "lost_bp", "cluster", "role", "cluster_size", "representative",
                         "cluster_label"], rows)
    print(f"table: {target}")


if __name__ == "__main__":
    main()
