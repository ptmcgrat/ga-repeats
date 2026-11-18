#!/usr/bin/env python3
"""Curated-family recovery (question 2): how many curated families does the library recover, and how well?

The measure of the RepeatModeler2 paper (Flynn et al. 2020, PNAS 117:9451): annotate the library's consensus
sequences with a curated library (RepeatMasker), then call each curated family
  perfect    one of our sequences matches it at >95% identity over >95% of its length;
  good       several of our sequences together do (>95% identity, >95% of its length);
  present    one or several match at >80% identity over >80% of its length;
  not found  below that (fragments may still be there).
Identity is 100 - RepeatMasker's "perc div." (mismatches; indels not counted); length is the curated consensus.

The curated families are the Repbase entries of the chosen species (--species; default the 165 cichlid entries:
Cichlidae and Oreochromis_niloticus), taken from the classified Repbase that repbase_match.py built
(<scratch>/repbase/<release>/repbase.fa and repbase_entries.tsv); run repbase_match.py first. Lengths are those of the
sequences searched; where a name occurs twice in the classified Repbase, the first copy is used and the duplicates are
listed (repeats_026: REX1-4_AFC had shown a coverage of 2.93 against the table's length).

Usage (from the repository root):
  python scripts/family_recovery.py --version RI_v5                                  # dry run
  python scripts/family_recovery.py --version RI_v5 --go
  python scripts/family_recovery.py --version RI_v5 --species Danio_rerio --tag zebrafish --go

Writes, never overwriting, in <scratch>/comparisons/family_recovery/: <version>_<tag>_families.tsv (one row per
curated family) and <version>_<tag>_summary.tsv (counts per class and category). Work in
<scratch>/<version>/family_recovery/<tag>/.
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import curated_annotation as ca  # noqa: E402

CATEGORIES = ["perfect", "good", "present", "not found"]
Hit = Tuple[str, str, float, int, int]   # our sequence, curated family, identity, consensus start, consensus end (1-based)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the library version, e.g. RI_v5")
    parser.add_argument("--species", nargs="+", default=["Cichlidae", "Oreochromis_niloticus"],
                        help="Repbase species of the curated families (default: Cichlidae Oreochromis_niloticus)")
    parser.add_argument("--tag", default="cichlid", help="name for this set of curated families (default cichlid)")
    parser.add_argument("--release", default="20181026", help="the classified Repbase release folder (default 20181026)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def parse_hits(lines: Iterable[str]) -> List[Hit]:
    """RepeatMasker .out lines -> (query, repeat, identity, consensus start, consensus end)."""
    out = []
    for line in lines:
        f = line.split()
        if len(f) < 15 or not f[0].isdigit():
            continue
        start, end = (int(f[11]), int(f[12])) if f[8] == "+" else (int(f[13]), int(f[12]))
        out.append((f[4], f[9], 100.0 - float(f[1]), min(start, end), max(start, end)))
    return out


def covered(spans: Iterable[Tuple[int, int]]) -> int:
    """Bases of a consensus covered by 1-based inclusive spans."""
    return c.merged_length(c.merge_by_sequence(("x", a - 1, b) for a, b in spans))


def classify(hits: List[Hit], length: int) -> Dict[str, object]:
    """Flynn et al. 2020's category for one curated family, with the best single sequence of ours."""
    high = [h for h in hits if h[2] > 95]
    mid = [h for h in hits if h[2] > 80]
    by_query: Dict[str, list] = defaultdict(list)
    for h in mid:
        by_query[h[0]].append(h)
    best, best_cov, best_id = "", 0.0, 0.0
    for q, hs in by_query.items():
        cov = covered((h[3], h[4]) for h in hs) / length
        ident = sum(h[2] * (h[4] - h[3] + 1) for h in hs) / sum(h[4] - h[3] + 1 for h in hs)
        if cov > best_cov:
            best, best_cov, best_id = q, cov, ident
    single_95 = max((covered((h[3], h[4]) for h in high if h[0] == q) / length for q in {h[0] for h in high}), default=0.0)
    union_95 = covered((h[3], h[4]) for h in high) / length
    union_80 = covered((h[3], h[4]) for h in mid) / length
    category = ("perfect" if single_95 > 0.95 else "good" if union_95 > 0.95 else "present" if union_80 > 0.80 else "not found")
    return {"category": category, "best": best, "best_identity": best_id, "best_coverage": best_cov,
            "coverage_95": union_95, "coverage_80": union_80, "sequences_80": len(by_query)}


def main(argv=None) -> None:
    args = parse_args(argv)
    lib_dir = args.scratch / "repbase" / args.release
    repbase_fa, entries_file = lib_dir / "repbase.fa", lib_dir / "repbase_entries.tsv"
    for f in (repbase_fa, entries_file):
        if not f.exists():
            c.stop(f"{f} not found: run repbase_match.py first (it builds the classified Repbase)")
    with open(entries_file) as handle:
        entries = [r for r in csv.DictReader(handle, delimiter="\t") if r["species"] in set(args.species)]
    if not entries:
        c.stop(f"no Repbase entries of {', '.join(args.species)} in {entries_file}")
    _, _, library = c.library_files(args.version)
    work = args.scratch / args.version / "family_recovery" / args.tag
    out_dir = args.scratch / "comparisons" / "family_recovery"
    targets = [out_dir / f"{args.version}_{args.tag}_families.tsv", out_dir / f"{args.version}_{args.tag}_summary.tsv"]
    existing = [t for t in targets if t.exists()]
    if existing:
        c.stop("already exists, never overwritten: " + ", ".join(map(str, existing)))
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    pa = max(1, args.threads // c.RMBLAST_CORES_PER_JOB)
    by_type = Counter(r["type"] for r in entries)
    print(f"curated families: {len(entries)} Repbase entries of {', '.join(args.species)} ("
          + ", ".join(f"{k} {v}" for k, v in by_type.most_common()) + f")\nour library: {library}\n"
          f"search in {work}; RepeatMasker -pa {pa}\nresults: {targets[0].name}, {targets[1].name} in {out_dir}")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    started = datetime.now()
    work.mkdir(parents=True)
    wanted = {r["name"] for r in entries}
    curated = work / "curated.fa"
    lengths, duplicates = ca.write_curated(repbase_fa, wanted, curated)
    if duplicates:
        print(f"note: {len(set(duplicates))} entry names occur more than once in {repbase_fa.name}; the first copy is used "
              f"(repeats_026): " + ", ".join(sorted(set(duplicates))[:10]))
    families = [f.name for f in c.read_families(library)]
    names = {f"q{i:06d}": name for i, name in enumerate(families, start=1)}
    sequences, current = {}, None
    with open(library) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split(None, 1)[0]
                sequences[current] = []
            elif current is not None:
                sequences[current].append(line)
    query = work / "query.fa"
    with open(query, "x") as handle:
        for short, name in names.items():
            handle.write(f">{short}\n{''.join(sequences[name])}")
    c.write_tsv(work / "query_names.tsv", ["query", "family"], names.items())
    command = c.container(args.sif, "RepeatMasker", "-lib", str(curated), "-engine", "rmblast", "-pa", str(pa),
                          "-nolow", "-no_is", "-dir", str(work), str(query))
    c.run(command, True, work / "repeatmasker.log")
    out_file = work / "query.fa.out"
    hits = parse_hits(open(out_file)) if out_file.exists() else []
    per_family: Dict[str, List[Hit]] = defaultdict(list)
    for h in hits:
        per_family[h[1]].append(h)

    rows, summary = [], Counter()
    for r in sorted(entries, key=lambda r: (r["type"], r["name"])):
        length = lengths.get(r["name"]) or int(r["length"])   # the sequence searched (repeats_026), else the table
        res = classify(per_family.get(r["name"], []), length)
        summary[(r["type"], res["category"])] += 1
        rows.append([r["name"], r["type"], r["subtype"], length, res["category"], names.get(res["best"], ""),
                     f"{res['best_identity']:.1f}" if res["best"] else "", f"{res['best_coverage']:.3f}",
                     f"{res['coverage_95']:.3f}", f"{res['coverage_80']:.3f}", res["sequences_80"]])
    types = [t for t, _ in by_type.most_common()]
    summary_rows = [[t] + [summary[(t, k)] for k in CATEGORIES] + [by_type[t]] for t in types]
    summary_rows.append(["ALL"] + [sum(summary[(t, k)] for t in types) for k in CATEGORIES] + [len(entries)])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["curated_family", "type", "subtype", "length", "category", "best_family", "best_identity",
                             "best_coverage", "coverage_by_95pct_matches", "coverage_by_80pct_matches", "our_sequences_over_80pct"], rows)
    c.write_tsv(targets[1], ["type"] + CATEGORIES + ["families"], summary_rows)
    with open(work / "record.json", "x") as handle:
        json.dump({"version": args.version, "library": str(library), "library_sha256": c.sha256(library),
                   "species": args.species, "curated_families": len(entries), "command": c.shell_text(command),
                   "started": f"{started:%Y-%m-%d %H:%M}", "finished": f"{datetime.now():%Y-%m-%d %H:%M}"}, handle, indent=2)
        handle.write("\n")
    print(f"\n{args.version} against {len(entries)} curated {args.tag} families (Flynn et al. 2020's categories):")
    print(f"  {'type':10s} " + " ".join(f"{k:>10s}" for k in CATEGORIES) + f" {'families':>9s}")
    for row in summary_rows:
        print(f"  {row[0]:10s} " + " ".join(f"{v:>10d}" for v in row[1:5]) + f" {row[5]:>9d}")
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
