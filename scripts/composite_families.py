#!/usr/bin/env python3
"""Composite families: consensus sequences that contain another class's element (question 2; read-only).

Two of the three published insertions (findings B5) and many curated cichlid families (B11) are best matched by families
labelled with another class. One explanation (C1) is a composite consensus: an LTR element carrying a hAT copy, a LINE
carrying a DNA transposon. This check looks for that in every family. It changes nothing in the library.

It reads the search family_recovery.py already ran: the version's consensus sequences annotated by RepeatMasker with the
curated cichlid Repbase entries (or, with --all-repbase, repbase_match.py's search against all of Repbase). For each
family, every curated entry matching it is joined over its fragments: the stretches of our consensus it covers, the
share of the entry it covers, identity (100 - RepeatMasker's divergence), and the entry's class. Classes are compared
as groups (DNA, LINE, LTR, SINE, Unknown, Other). A family is
  composite  labelled with a class, carrying a stretch of at least --min-length bp (default 300) that matches, at
             --min-identity % or more (default 90), at least --min-entry (default 0.8) of a curated entry of another class;
  relabel    labelled with a class, with another class's entry matching at least --min-family (default 0.8) of the
             family itself at that identity: the family is that element, wrongly labelled;
  unknown    labelled Unknown, with such a match covering most of the entry or of the family: a class it could take;
  partial    another class's match passing identity and length but neither threshold: listed, not flagged.
The positive controls (--controls; default the families at the published mitfa and tbx2a insertions) must come out
composite or relabel; the printout says whether they do. Each family's footprint is the bases it masks in --genome
(default DT_004f; one label per base).

Usage (from the repository root):
  python scripts/composite_families.py --version RI_v5          # dry run
  python scripts/composite_families.py --version RI_v5 --go

Writes, never overwriting, in <scratch>/comparisons/composite/: <version>_<search>_matches.tsv (one row per family and
other-class entry) and <version>_<search>_families.tsv (one row per family with an other-class match).
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

CONTROLS = ["MZ_MZ4f_ltr-1_family-276", "MZ_NMetZ_ltr-1_family-392", "MZ_MZ4f_ltr-1_family-769"]
ORDER = {"relabel": 0, "composite": 1, "unknown": 2, "partial": 3}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the library version, e.g. RI_v5")
    parser.add_argument("--all-repbase", action="store_true", help="use repbase_match.py's search against all of Repbase")
    parser.add_argument("--tag", default="cichlid", help="family_recovery.py's set of curated entries (default cichlid)")
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-length", type=int, default=300)
    parser.add_argument("--min-entry", type=float, default=0.8)
    parser.add_argument("--min-family", type=float, default=0.8)
    parser.add_argument("--genome", default="DT_004f", help="genome masked by the version, for footprints (default DT_004f)")
    parser.add_argument("--controls", nargs="*", default=CONTROLS, help="families that must be flagged")
    parser.add_argument("--release", default="20181026")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def parse_hits(lines: Iterable[str]) -> List[Tuple[str, str, str, float, int, int, int, int]]:
    """RepeatMasker .out -> (query, entry, entry label, identity, query begin, query end, entry start, entry end)."""
    out = []
    for line in lines:
        f = line.split()
        if len(f) < 15 or not f[0].isdigit():
            continue
        es, ee = (int(f[11]), int(f[12])) if f[8] == "+" else (int(f[13]), int(f[12]))
        out.append((f[4], f[9], f[10], 100.0 - float(f[1]), int(f[5]), int(f[6]), min(es, ee), max(es, ee)))
    return out


def span(intervals: Iterable[Tuple[int, int]]) -> int:
    """Bases covered by 1-based inclusive intervals."""
    return c.merged_length(c.merge_by_sequence(("x", a - 1, b) for a, b in intervals))


def join_matches(hits, entry_lengths: Dict[str, int]) -> Dict[Tuple[str, str], dict]:
    """(query, entry) -> label, identity (weighted by bases), query bases covered, entry share covered, query stretches."""
    parts: Dict[Tuple[str, str], list] = defaultdict(list)
    labels = {}
    for q, entry, label, ident, qb, qe, es, ee in hits:
        parts[(q, entry)].append((ident, qb, qe, es, ee))
        labels[(q, entry)] = label
    out = {}
    for key, ps in parts.items():
        bases = sum(p[2] - p[1] + 1 for p in ps)
        out[key] = {"label": labels[key],
                    "identity": sum(p[0] * (p[2] - p[1] + 1) for p in ps) / bases,
                    "query_bp": span((p[1], p[2]) for p in ps),
                    "entry_share": span((p[3], p[4]) for p in ps) / entry_lengths[key[1]] if entry_lengths.get(key[1]) else 0.0,
                    "stretches": sorted((p[1], p[2]) for p in ps)}
    return out


def classify(own_label: str, own_length: int, match: dict, args) -> str:
    """'composite', 'relabel', 'unknown', 'partial' or '' for one family's match to one curated entry."""
    own, other = la.group(own_label), la.group(match["label"])
    if match["identity"] < args.min_identity or match["query_bp"] < args.min_length:
        return ""
    family_share = match["query_bp"] / own_length if own_length else 0.0
    if own == "Unknown":
        return "unknown" if other != "Unknown" and (family_share >= args.min_family or match["entry_share"] >= args.min_entry) else ""
    if other == own or other == "Unknown":
        return ""
    if family_share >= args.min_family:
        return "relabel"
    if match["entry_share"] >= args.min_entry:
        return "composite"
    return "partial"


def fasta_lengths_plain(path: Path) -> Dict[str, int]:
    """Entry name (before any #) -> sequence length, for a plain FASTA."""
    lengths: Dict[str, int] = {}
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split("#", 1)[0].split()[0]
                lengths.setdefault(name, 0)
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def main(argv=None) -> None:
    args = parse_args(argv)
    if args.all_repbase:
        work, search = args.scratch / args.version / "repbase", "allrepbase"
        lib = args.scratch / "repbase" / args.release / "repbase.fa"
    else:
        work, search = args.scratch / args.version / "family_recovery" / args.tag, args.tag
        lib = work / "curated.fa"
    out_file, names_file = work / "query.fa.out", work / "query_names.tsv"
    for f in (out_file, names_file, lib):
        if not f.exists():
            c.stop(f"{f} not found: run {'repbase_match.py' if args.all_repbase else 'family_recovery.py'} first")
    bed = c.version_folder(args.version) / "masked" / args.genome / f"{args.genome}.bed"
    out_dir = args.scratch / "comparisons" / "composite"
    targets = [out_dir / f"{args.version}_{search}_matches.tsv", out_dir / f"{args.version}_{search}_families.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    print(f"{args.version}: families against {'all of Repbase' if args.all_repbase else f'the curated {args.tag} entries'} "
          f"({out_file}); composite: >= {args.min_length} bp at >= {args.min_identity:g}% over >= {100 * args.min_entry:.0f}% of "
          f"an entry of another class; relabel: over >= {100 * args.min_family:.0f}% of the family")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    _, _, library = c.library_files(args.version)
    families = {f.name: f for f in c.read_families(library)}
    with open(names_file) as handle:
        names = {r["query"]: r["family"] for r in csv.DictReader(handle, delimiter="\t")}
    with open(out_file) as handle:
        matches = join_matches(parse_hits(handle), fasta_lengths_plain(lib))
    footprint: Counter = Counter()
    if bed.exists():
        print(f"reading {args.genome}'s masking for footprints ...", flush=True)
        for _s, a, b, fam, _l in la.one_label_per_base(la.read_labelled_bed(bed)):
            footprint[fam] += b - a
    rows, status = [], {}
    for (q, entry), m in sorted(matches.items()):
        name = names.get(q)
        fam = families.get(name)
        if fam is None:
            continue
        kind = classify(fam.label, fam.length, m, args)
        if not kind:
            continue
        rows.append([fam.family, fam.label, fam.length, kind, entry, m["label"], f"{m['identity']:.1f}", m["query_bp"],
                     f"{m['query_bp'] / fam.length:.3f}", f"{m['entry_share']:.3f}",
                     ";".join(f"{a}-{b}" for a, b in m["stretches"]), footprint.get(fam.family, 0)])
        if fam.family not in status or ORDER[kind] < ORDER[status[fam.family][0]]:
            status[fam.family] = (kind, fam)
    fam_rows = []
    for name, (kind, fam) in sorted(status.items(), key=lambda kv: (ORDER[kv[1][0]], -footprint.get(kv[0], 0))):
        origin = c.origin(name)
        others = sorted({r[5] for r in rows if r[0] == name})
        fam_rows.append([name, fam.label, fam.length, origin, kind, ",".join(others), footprint.get(name, 0)])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["family", "label", "length", "kind", "entry", "entry_label", "identity", "family_bp_matched",
                             "family_share", "entry_share", "family_stretches", "family_footprint_bp"], rows)
    c.write_tsv(targets[1], ["family", "label", "length", "origin", "status", "other_class_entries", "footprint_bp"], fam_rows)

    genome_bp = sum(c.fasta_lengths(c.genome_path(args.genome, c.panel_accession(args.genome))).values()) if bed.exists() else 0
    pct = lambda n: f"{100 * n / genome_bp:.2f}% of {args.genome}" if genome_bp else "footprint not measured"  # noqa: E731
    print(f"\n{len(families):,} families; with another class's curated entry inside: {len(status):,}")
    for kind in ("relabel", "composite", "unknown", "partial"):
        group = [r for r in fam_rows if r[4] == kind]
        if group:
            print(f"  {kind:9s} {len(group):>4}  ({', '.join(f'{k} {v}' for k, v in Counter(la.group(r[1]) for r in group).most_common())}; "
                  f"LTR pipeline {sum(1 for r in group if r[3] == 'LTR pipeline')}); masking {pct(sum(r[6] for r in group))}")
    top = [r for r in fam_rows if r[4] in ("composite", "relabel")][:12]
    if top:
        print("largest flagged families by footprint:")
        for r in top:
            print(f"  {r[0]:30s} {r[1]:18s} {r[4]:9s} contains {r[5]}  ({r[6] / 1e3:,.0f} kb)")
    print("positive controls: " + "; ".join(
        f"{f} -> {status[f][0] if f in status else 'not flagged'}" for f in args.controls))
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
