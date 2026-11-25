#!/usr/bin/env python3
"""Import a library curated by MCHelper or TEtrimmer as a new version (question 2).

Reads the final library from a curate_library.py work folder and writes a version whose library follows our
conventions, so mask_genomes.py and every check treat it like our own versions:
  names    the tool's family names (they derive from ours and keep the build prefix), with characters other than
           letters, digits, '_', '-' and '.' replaced by '_', and made unique where needed
  labels   RepeatMasker's form. TEtrimmer's labels are RepeatMasker's already (its classifier is RepeatClassifier);
           only the first word of each header is kept and its evaluation (Perfect, Good, ...) goes to the table.
           MCHelper writes Wicker's classification, mapped by order: CLASSI/LTR -> LTR, CLASSI/LINE -> LINE,
           CLASSI/SINE -> SINE, CLASSI/DIRS -> LTR/DIRS, CLASSI/PLE -> PLE, CLASSI/TRIM or LARD -> LTR,
           CLASSII/TIR or MITE -> DNA, CLASSII/HELITRON -> RC/Helitron, CLASSII/MAVERICK -> DNA/Maverick,
           CLASSII/CRYPTON -> DNA/Crypton; its superfamily follows ('GYPSY' -> 'Gypsy', 'TC1-MARINER' -> 'TcMar', ...).
           A class with no order, an order not in this list, or 'unclassified' becomes Unknown, never a guess.
The table keeps every family's original name and label. A pilot run is refused unless --pilot-ok.

Writes <version folder>/library/<version>-families.prefixed.fa, import_table.tsv and curation_record.json (the tool,
its version and environment, the source version, the run's own record, the label mapping used, counts). The version
folder must not exist yet. Work in <scratch>/<version>/import/.

Usage (from the repository root):
  python scripts/import_curated.py --from <scratch>/curation/tetrimmer/RI_v1_MZ_MZ4f --version RI_v10          # dry run
  python scripts/import_curated.py --from <scratch>/curation/tetrimmer/RI_v1_MZ_MZ4f --version RI_v10 --go
"""
import argparse
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import curate_library as cl  # noqa: E402
import label_agreement as la  # noqa: E402

ORDERS, SUPERFAMILIES = la.ORDERS, la.SUPERFAMILIES


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from", dest="source", type=Path, required=True, help="a curate_library.py work folder")
    parser.add_argument("--version", required=True, help="the new version, e.g. RI_v10")
    parser.add_argument("--pilot-ok", action="store_true", help="allow importing a pilot run")
    c.add_common_args(parser)
    return parser.parse_args(argv)


mchelper_label = la.wicker_label


def parse_header(line: str) -> Tuple[str, str, str]:
    """'>name#label  extra fields' -> (name, label, the extra fields joined by spaces)."""
    words = line[1:].split()
    name, _, label = words[0].partition("#") if words else ("", "", "")
    return name, label or "Unknown", " ".join(words[1:])


def clean_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", name)


def convert(lines, tool: str) -> Tuple[List[Tuple[str, str, str]], List[list]]:
    """Tool FASTA lines -> records (new name, new label, sequence) and table rows."""
    records, rows, seen = [], [], Counter()
    name = label = extra = None
    seq: List[str] = []

    def flush():
        if name is None:
            return
        base = clean_name(name)
        seen[base] += 1
        new = base if seen[base] == 1 else f"{base}_dup{seen[base]}"
        new_label = mchelper_label(label) if tool == "mchelper" else (label or "Unknown")
        evaluation = extra if tool == "tetrimmer" else ("incomplete" if name.endswith("_inc") else "")
        sequence = "".join(seq)
        records.append((new, new_label, sequence))
        rows.append([new, new_label, name, label, evaluation, len(sequence)])
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            name, label, extra = parse_header(line)
            seq = []
        elif name is not None:
            seq.append(line.strip())
    flush()
    return records, rows


def main(argv=None) -> None:
    args = parse_args(argv)
    record_file = args.source / "record.json"
    if not record_file.exists():
        c.stop(f"{record_file} not found: {args.source} is not a finished curate_library.py run")
    with open(record_file) as handle:
        run = json.load(handle)
    if run.get("pilot") and not args.pilot_ok:
        c.stop(f"{args.source} is a pilot ({run['pilot']} families); add --pilot-ok to import it anyway")
    tool = run["tool"]
    tool_output = args.source / cl.OUTPUT[tool]
    target = c.version_folder(args.version)
    if target.exists():
        c.stop(f"{target} already exists; a version is never changed")
    work = args.scratch / args.version / "import"
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    print(f"{run['tool_version']} on {run['version']} against {run['assembly']} ({run['families_in']:,} families in, "
          f"{run['families_out']:,} out) -> {args.version}\nfrom: {tool_output}\nto:   {target / 'library'}")
    if not args.go:
        print("dry run: nothing was written. Add --go to import.")
        return

    with open(tool_output) as handle:
        records, rows = convert(handle, tool)
    work.mkdir(parents=True)
    library = work / f"{args.version}-families.prefixed.fa"
    with open(library, "x") as out:
        for new, label, seq in records:
            out.write(f">{new}#{label}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    c.write_tsv(work / "import_table.tsv", ["family", "label", "tool_name", "tool_label", "evaluation", "length"], rows)
    mapping = Counter(f"{r[3]} -> {r[1]}" for r in rows if r[3] != r[1])
    by_class = Counter(la.group(r[1]) for r in rows)
    record = {"version": args.version, "kind": "curated library", "tool": tool, "tool_version": run["tool_version"],
              "environment": run["environment"], "source_version": run["version"], "assembly": run["assembly"],
              "tool_output": str(tool_output), "tool_output_sha256": c.sha256(tool_output), "run": run,
              "counts": {"families": len(records), "by_class": dict(by_class.most_common()),
                         "renamed": sum(1 for r in rows if r[0] != r[2]),
                         "unknown_after_mapping": sum(1 for r in rows if r[1] == "Unknown")},
              "label_mapping": dict(mapping.most_common()),
              "written": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(work / "curation_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    for name in (library.name, "import_table.tsv", "curation_record.json"):
        c.copy_new(work / name, target / "library" / name)
    print(f"\n{len(records):,} families ({', '.join(f'{k} {v}' for k, v in by_class.most_common())}); "
          f"labels mapped: {sum(mapping.values()):,}; renamed: {record['counts']['renamed']:,}")
    for k, v in mapping.most_common(8):
        print(f"  {v:>5}  {k}")
    print(f"written: {target / 'library'}/\nnext: python scripts/mask_genomes.py --version {args.version} --assemblies AC_fAstCal1_2")


if __name__ == "__main__":
    main()
