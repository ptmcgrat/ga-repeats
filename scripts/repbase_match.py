#!/usr/bin/env python3
"""Our families against curated Repbase: the label check, and transposon support for the gene filter.

1. builds a classified Repbase library once, in <scratch>/repbase/<release>/, from the RepeatMasker
   edition's two halves joined by entry name: RMRBSeqs.embl (the sequences) and RMRBMeta.embl (Type,
   SubType, Species). repbase.fa has headers NAME#Type/SubType; repbase_entries.tsv lists every entry's
   class, species and length. Never overwritten;
2. searches a version's families against it with RepeatMasker, as compare_libraries.py does;
3. one row per family: its best Repbase match (highest score), that entry's class and species, how much
   of the family the match covers, whether the top classes agree, and whether the family has Repbase
   transposon support: its best match is a transposon class (TE_CLASSES) covering at least --support
   (default 0.50) of it. A summary by our class gives agreement among families matched over --coverage
   (default 0.80) of their length.

Inputs: --repbase-dir holding RMRBSeqs.embl and RMRBMeta.embl (default ~/Repbase/20181026).
Results: <scratch>/comparisons/repbase/<version>_best_match.tsv and _by_class.tsv; the search runs in
<scratch>/<version>/repbase/. Never overwritten.

Example:
  python scripts/repbase_match.py --version RI_v5          # dry run
  python scripts/repbase_match.py --version RI_v5 --go
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import compare_libraries as cl  # noqa: E402

TE_CLASSES = ("DNA", "LINE", "LTR", "SINE", "RC", "Retroposon", "PLE")
DEFAULT_REPBASE = Path.home() / "Repbase" / "20181026"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the library version to check, e.g. RI_v5")
    parser.add_argument("--repbase-dir", type=Path, default=DEFAULT_REPBASE,
                        help="folder with RMRBSeqs.embl and RMRBMeta.embl (default ~/Repbase/20181026)")
    parser.add_argument("--support", type=float, default=0.50,
                        help="coverage of the family a transposon match needs to count as support (default 0.50)")
    parser.add_argument("--coverage", type=float, default=0.80,
                        help="coverage for the label-agreement summary (default 0.80)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def embl_entries(lines: Iterable[str]) -> Iterable[Tuple[str, Dict[str, str], str]]:
    """(name, annotations, sequence) per entry of a RepeatMasker EMBL file (either half)."""
    name, notes, seq, in_seq = None, {}, [], False
    for line in lines:
        if line.startswith("ID "):
            name, notes, seq, in_seq = line.split()[1], {}, [], False
        elif line.startswith("//"):
            if name is not None:
                yield name, notes, "".join(seq).upper()
            name, in_seq = None, False
        elif line.startswith("SQ "):
            in_seq = True
        elif in_seq:
            seq.append("".join(ch for ch in line if ch.isalpha()))
        elif line.startswith("CC ") and ":" in line:
            key, _, value = line[2:].strip().partition(":")
            if key in ("Type", "SubType", "Species"):
                notes[key] = value.strip()


def build_repbase(seqs_lines: Iterable[str], meta_lines: Iterable[str]) -> List[Tuple[str, str, str, str, str]]:
    """(name, type, subtype, species, sequence) for every entry present in both halves."""
    meta = {name: notes for name, notes, _seq in embl_entries(meta_lines)}
    out, missing = [], []
    for name, _notes, seq in embl_entries(seqs_lines):
        notes = meta.get(name)
        if notes is None:
            missing.append(name)
            continue
        out.append((name, notes.get("Type") or "Unknown", notes.get("SubType", ""), notes.get("Species", ""), seq))
    if missing:
        c.stop(f"{len(missing)} sequence entries have no metadata entry, e.g. {', '.join(missing[:3])}")
    if len(out) != len(meta):
        c.stop(f"the metadata has {len(meta)} entries but only {len(out)} have sequences")
    return out


def label_of(entry_type: str, subtype: str) -> str:
    return f"{entry_type}/{subtype}" if subtype else entry_type


def summarise(rows: List[list], coverage: float) -> List[list]:
    """By our top class: families, with a match, matched over the coverage, of those same top class."""
    by: Dict[str, List[list]] = defaultdict(list)
    for r in rows:
        by[r[2]].append(r)
        by["ALL"].append(r)
    out = []
    for cls in sorted(k for k in by if k != "ALL") + ["ALL"]:
        group = by[cls]
        matched = [r for r in group if r[4]]
        well = [r for r in matched if float(r[10]) >= coverage]
        same = [r for r in well if r[12] == "yes"]
        out.append([cls, len(group), len(matched), len(well), len(same),
                    f"{100 * len(same) / len(well):.1f}" if well else ""])
    return out


def main(argv=None) -> None:
    args = parse_args(argv)
    seqs_file, meta_file = args.repbase_dir / "RMRBSeqs.embl", args.repbase_dir / "RMRBMeta.embl"
    for f in (seqs_file, meta_file):
        if not f.exists():
            c.stop(f"{f} not found")
    release = args.repbase_dir.name
    lib_dir = args.scratch / "repbase" / release
    repbase_fa = lib_dir / "repbase.fa"
    _, _, library = c.library_files(args.version)
    work = args.scratch / args.version / "repbase"
    out_dir = args.scratch / "comparisons" / "repbase"
    targets = [out_dir / f"{args.version}_best_match.tsv", out_dir / f"{args.version}_by_class.tsv"]
    pa = max(1, args.threads // c.RMBLAST_CORES_PER_JOB)
    print(f"Repbase: {seqs_file} and {meta_file}\nclassified library: {repbase_fa}"
          f"{' (exists, reused)' if repbase_fa.exists() else ' (built now)'}\nfamilies: {library}\n"
          f"search in {work}; RepeatMasker -pa {pa}\nresults: {targets[0].name}, {targets[1].name} in {out_dir}")
    existing = [p for p in targets if p.exists()]
    if existing:
        c.stop("already exists, never overwritten: " + ", ".join(map(str, existing)))
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    started = datetime.now()
    if not repbase_fa.exists():
        with open(seqs_file) as s, open(meta_file) as m:
            entries = build_repbase(s, m)
        lib_dir.mkdir(parents=True, exist_ok=True)
        with open(repbase_fa, "x") as handle:
            for name, t, st, _sp, seq in entries:
                handle.write(f">{name}#{label_of(t, st)}\n" + "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)))
        c.write_tsv(lib_dir / "repbase_entries.tsv", ["name", "type", "subtype", "species", "length"],
                    [[n, t, st, sp, len(seq)] for n, t, st, sp, seq in entries])
        with open(lib_dir / "repbase_record.json", "x") as handle:
            json.dump({"RMRBSeqs.embl": {"path": str(seqs_file), "sha256": c.sha256(seqs_file)},
                       "RMRBMeta.embl": {"path": str(meta_file), "sha256": c.sha256(meta_file)},
                       "entries": len(entries), "by_type": dict(Counter(e[1] for e in entries).most_common()),
                       "built": f"{datetime.now():%Y-%m-%d %H:%M}"}, handle, indent=2)
            handle.write("\n")
        print(f"classified Repbase: {len(entries):,} entries; "
              + ", ".join(f"{t} {n:,}" for t, n in Counter(e[1] for e in entries).most_common(8)))
    species = {}
    with open(lib_dir / "repbase_entries.tsv") as handle:
        next(handle)
        for line in handle:
            f = line.rstrip("\n").split("\t")
            species[f[0]] = f[3]

    families = {f.name: f for f in c.read_families(library)}
    names = {f"q{i:06d}": name for i, name in enumerate(families, start=1)}
    work.mkdir(parents=True)
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
    c.run(c.container(args.sif, "RepeatMasker", "-lib", str(repbase_fa), "-engine", "rmblast", "-pa", str(pa),
                      "-nolow", "-no_is", "-dir", str(work), str(query)), True, work / "repeatmasker.log")

    base = cl.best_matches(work / "query.fa.out", names, families, args.coverage)
    rows = []
    for name, label, cls, length, rep, rep_label, score, div, cov, _well in base:
        rep_cls = c.top_class(rep_label) if rep else ""
        support = "yes" if rep and rep_cls in TE_CLASSES and float(cov) >= args.support else "no"
        same = ("yes" if rep_cls == cls else "no") if rep else ""
        rows.append([name, label, cls, length, rep, rep_label, rep_cls, species.get(rep, ""), score, div,
                     cov, support, same])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["family", "label", "class", "length", "repbase_entry", "repbase_label", "repbase_class",
                             "repbase_species", "score", "divergence", "coverage_of_family", "transposon_support",
                             "same_class"], rows)
    summary = summarise(rows, args.coverage)
    c.write_tsv(targets[1], ["class", "families", "with_a_match", f"matched_over_{args.coverage:.2f}",
                             "of_those_same_class", "percent_same_class"], summary)
    print(f"\n{args.version}: {len(rows):,} families; with a Repbase match {sum(1 for r in rows if r[4]):,}; "
          f"with transposon support {sum(1 for r in rows if r[11] == 'yes'):,}")
    print(f"labels, among families matched over {args.coverage:.0%} of their length:")
    for cls, n, matched, well, same, pct in summary:
        print(f"  {cls:12s} {n:>6,} families, {well:>5,} well matched, same class {same:>5,} ({pct}%)" if well else
              f"  {cls:12s} {n:>6,} families, none well matched")
    print(f"written: {targets[0].name}, {targets[1].name} in {out_dir}; run time "
          f"{int((datetime.now() - started).total_seconds() // 60)} min")


if __name__ == "__main__":
    main()
