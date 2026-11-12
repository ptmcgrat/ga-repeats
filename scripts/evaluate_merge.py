#!/usr/bin/env python3
"""What did merging cost, and what did it add? (merge rehearsal)

On one genome masked by a merged version (merge_libraries.py) and by each of its sources, compare
the merged library's masking with the union of the sources' maskings (a base counts if any source
masks it; masking_gain.py's total). Reports, as base pairs and percent of the genome:
  - masked by the merged library, by any source, by both;
  - lost: masked by a source but not by the merged library (what de-duplication cost);
  - gained: masked by the merged library but by no source;
  - lost and gained split by class. Each base takes one label: for the merged library, its
    highest-scoring annotation; for the sources, the highest-scoring annotation among all of them.
The sources are read from the merged version's library/merge_record.json.

Inputs: <version folder>/masked/<assembly>/<assembly>.bed for the merged version and each source.
Result: <scratch>/comparisons/merge/<version>_<assembly>.tsv. Never overwritten.

Example:
  python scripts/evaluate_merge.py --version RI_v4 --assembly DT_004f          # dry run
  python scripts/evaluate_merge.py --version RI_v4 --assembly DT_004f --go
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a merged version, e.g. RI_v4")
    parser.add_argument("--assembly", required=True, help="a genome masked by it and by all its sources, e.g. DT_004f")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def merge_effect(merged_raw: List[la.Scored], source_raws: List[List[la.Scored]]) -> Dict[str, object]:
    """Base-pair counts for the merged masking against the union of the source maskings."""
    merged_rows = la.one_label_per_base(merged_raw)
    union_rows = la.one_label_per_base([r for rows in source_raws for r in rows])
    m_all, m_groups = la.group_sets(merged_rows)
    u_all, u_groups = la.group_sets(union_rows)
    return {
        "merged": c.merged_length(m_all), "union": c.merged_length(u_all), "both": c.intersect_length(m_all, u_all),
        "lost": c.subtract_length(u_all, m_all), "gained": c.subtract_length(m_all, u_all),
        "lost_by_class": {g: c.subtract_length(s, m_all) for g, s in u_groups.items()},
        "gained_by_class": {g: c.subtract_length(s, u_all) for g, s in m_groups.items()},
    }


def main(argv=None) -> None:
    args = parse_args(argv)
    record_file = c.version_folder(args.version) / "library" / "merge_record.json"
    if not record_file.exists():
        c.stop(f"{args.version} is not a merged version: {record_file} not found")
    sources = [s["version"] for s in json.loads(record_file.read_text())["sources"]]
    beds = {}
    for version in [args.version] + sources:
        bed = c.version_folder(version) / "masked" / args.assembly / f"{args.assembly}.bed"
        if not bed.exists():
            c.stop(f"{version} has not masked {args.assembly}: {bed} not found")
        beds[version] = bed
        print(f"{version}: {bed}")
    target = args.scratch / "comparisons" / "merge" / f"{args.version}_{args.assembly}.tsv"
    print(f"merged {args.version} against the union of {', '.join(sources)}; result: {target}")
    if target.exists():
        c.stop(f"{target} already exists, never overwritten")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    genome_bp = sum(c.fasta_lengths(genome).values())
    raws = {}
    for version, bed in beds.items():
        print(f"reading {version} ...", flush=True)
        raws[version] = la.read_labelled_bed(bed)
    effect = merge_effect(raws[args.version], [raws[s] for s in sources])
    pct = lambda bp: f"{100 * bp / genome_bp:.3f}"
    rows = [[args.assembly, args.version, k, "ALL", effect[k], pct(effect[k])]
            for k in ("merged", "union", "both", "lost", "gained")]
    for kind in ("lost", "gained"):
        for g in la.GROUPS:
            bp = effect[f"{kind}_by_class"].get(g, 0)
            rows.append([args.assembly, args.version, kind, g, bp, pct(bp)])
    for s in sources:
        bp = c.merged_length(c.merge_by_sequence((r[0], r[1], r[2]) for r in raws[s]))
        rows.append([args.assembly, args.version, f"masked by {s} alone", "ALL", bp, pct(bp)])
    target.parent.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["assembly", "merged_version", "measure", "class", "bp", "percent_of_genome"], rows)

    print(f"\n{args.assembly}: merged library {pct(effect['merged'])}%, any source {pct(effect['union'])}%, "
          f"both {pct(effect['both'])}%")
    print(f"  lost by merging (a source masks it, the merged library does not): {pct(effect['lost'])}% — "
          + ", ".join(f"{g} {pct(v)}" for g, v in sorted(effect["lost_by_class"].items(), key=lambda kv: -kv[1]) if v))
    print(f"  gained by merging (only the merged library masks it): {pct(effect['gained'])}% — "
          + ", ".join(f"{g} {pct(v)}" for g, v in sorted(effect["gained_by_class"].items(), key=lambda kv: -kv[1]) if v))
    print(f"written: {target}")


if __name__ == "__main__":
    main()
