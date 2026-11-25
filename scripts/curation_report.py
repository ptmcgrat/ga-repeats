#!/usr/bin/env python3
"""What a curation tool did to each family (question 2; a curate_library.py work folder, pilot or full).

Each input family is traced to its outputs (curate_library.source_family: the tools rename hyphens as underscores and
add suffixes for subfamilies and incomplete copies) and reported as
  unchanged  one output with the same sequence as the input
  changed    one output with another sequence (trimmed, extended or rebuilt): length before and after
  split      two or more outputs (subfamilies)
  absent     no output
with its label before and after (MCHelper's Wicker labels in RepeatMasker's form) and TEtrimmer's evaluation. Summaries:
fates; class changes (before -> after, each output counted); Unknown before and after; TEtrimmer's evaluations; the
--watch families in full (default: the composite families at the published tbx2a and mitfa insertions and the
satellite family). It reads only the work folder.

Usage (from the repository root):
  python scripts/curation_report.py --from <scratch>/curation/tetrimmer/RI_v1_MZ_MZ4f_pilot100          # dry run
  python scripts/curation_report.py --from <scratch>/curation/tetrimmer/RI_v1_MZ_MZ4f_pilot100 --go

Writes, never overwriting, <work folder>/report.tsv (one row per input family).
"""
import argparse
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import curate_library as cl  # noqa: E402
import label_agreement as la  # noqa: E402

Record = Tuple[str, str, str, str]      # name, label, extra header fields, sequence


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from", dest="source", type=Path, required=True, help="a curate_library.py work folder")
    parser.add_argument("--watch", nargs="*", default=cl.INCLUDE, help="families to show in full")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def read_records(path: Path) -> List[Record]:
    out: List[Record] = []
    name = label = extra = None
    seq: List[str] = []
    with c.open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, label, extra, "".join(seq).upper()))
                words = line[1:].split()
                name, _, label = words[0].partition("#")
                label, extra, seq = label or "Unknown", " ".join(words[1:]), []
            elif name is not None:
                seq.append(line.strip())
    if name is not None:
        out.append((name, label, extra, "".join(seq).upper()))
    return out


def fates(tool: str, inputs: List[Record], outputs: List[Record]) -> Tuple[Dict[str, dict], List[str]]:
    """Input family -> its fate and outputs; and the outputs no input could be found for."""
    by_norm = {cl.normalize(r[0]): r[0] for r in inputs}
    traced: Dict[str, List[Record]] = defaultdict(list)
    orphans = []
    for r in outputs:
        src = cl.source_family(r[0], by_norm)
        if src is None:
            orphans.append(r[0])
        else:
            traced[src].append(r)
    result = {}
    for name, label, _extra, seq in inputs:
        outs = traced.get(name, [])
        if not outs:
            fate = "absent"
        elif len(outs) > 1:
            fate = "split"
        else:
            fate = "unchanged" if outs[0][3] == seq else "changed"
        result[name] = {"label": label, "length": len(seq), "fate": fate, "outputs": outs,
                        "labels_after": [la.wicker_label(o[1]) if tool == "mchelper" else o[1] for o in outs]}
    return result, orphans


def main(argv=None) -> None:
    args = parse_args(argv)
    rec_file = args.source / "record.json"
    if not rec_file.exists():
        c.stop(f"{rec_file} not found: {args.source} is not a finished curate_library.py run")
    with open(rec_file) as handle:
        run = json.load(handle)
    tool = run["tool"]
    target = args.source / "report.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    print(f"{run['tool_version']} on {run['version']} ({run['families_in']:,} families"
          + (f", pilot of {run['pilot']}" if run.get("pilot") else "") + f"), {run['hours']} h; report: {target}")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to report.")
        return

    inputs = read_records(args.source / "input.fa")
    outputs = read_records(args.source / cl.OUTPUT[tool])
    result, orphans = fates(tool, inputs, outputs)
    fate_n = Counter(r["fate"] for r in result.values())
    print(f"\nfamilies in {len(inputs):,}, out {len(outputs):,}: " + ", ".join(f"{k} {fate_n.get(k, 0)}" for k in
          ("unchanged", "changed", "split", "absent")) + (f"; outputs traced to no input: {len(orphans)}" if orphans else ""))
    changes: Counter = Counter()
    for r in result.values():
        for after in r["labels_after"]:
            changes[(la.group(r["label"]), la.group(after))] += 1
    print("class before -> after (each output counted):")
    for (b, a), n in sorted(changes.items(), key=lambda kv: (kv[0][0] != kv[0][1], -kv[1])):
        print(f"  {b:8s} -> {a:8s} {n:>4}" + ("" if a == b else "   changed"))
    unknown_in = sum(1 for r in inputs if la.group(r[1]) == "Unknown")
    unknown_out = sum(1 for r in result.values() for a in r["labels_after"] if la.group(a) == "Unknown")
    print(f"Unknown: {unknown_in} families in, {unknown_out} outputs")
    changed = [(len(r["outputs"][0][3]), r["length"]) for r in result.values() if r["fate"] == "changed"]
    if changed:
        ratios = [a / b for a, b in changed if b]
        print(f"changed families: length after / before, median {statistics.median(ratios):.2f} "
              f"(range {min(ratios):.2f}-{max(ratios):.2f})")
    if tool == "tetrimmer":
        evals = Counter((o[2].split()[-1] if o[2] else "no evaluation") for r in result.values() for o in r["outputs"])
        print("TEtrimmer's evaluations: " + ", ".join(f"{k} {v}" for k, v in evals.most_common()))
    print("watched families:")
    for name in args.watch:
        r = result.get(name)
        if r is None:
            print(f"  {name}: not in this run")
            continue
        print(f"  {name} ({r['label']}, {r['length']:,} bp): {r['fate']}")
        for o, after in zip(r["outputs"], r["labels_after"]):
            print(f"      {o[0]}  {after}  {len(o[3]):,} bp" + (f"  [{o[2]}]" if o[2] else ""))
    rows = [[name, r["label"], r["length"], r["fate"], len(r["outputs"]), ";".join(o[0] for o in r["outputs"]),
             ";".join(r["labels_after"]), ";".join(str(len(o[3])) for o in r["outputs"]),
             ";".join(o[2] for o in r["outputs"])] for name, r in result.items()]
    c.write_tsv(target, ["family", "label_before", "length_before", "fate", "outputs", "output_names", "labels_after",
                         "lengths_after", "evaluations"], rows)
    if orphans:
        print(f"outputs traced to no input (first 5): {', '.join(orphans[:5])}")
    print(f"table: {target}")


if __name__ == "__main__":
    main()
