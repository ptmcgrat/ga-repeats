#!/usr/bin/env python3
"""A version's library plus an added library, combined without clustering (question 2; Quah et al. 2025's procedure).

Quah et al. combined their de novo libraries with Repbase-derived RepeatMasker libraries and masked with the result. This
writes a new version whose library is the version's library followed by the added entries (for instance the Repbase
lineage subset of repbase_lineage.py), both unchanged: where several families match a locus, RepeatMasker keeps the
best-scoring one. Family names must not collide (the Repbase_ prefix of repbase_lineage.py guarantees it).

Writes <version folder>/library/<new>-families.prefixed.fa and combine_record.json (the two sources and their sha256,
the added library's own record where one lies beside it, counts). The new version folder must not exist yet. Work in
<scratch>/<new>/combine/.

Usage (from the repository root):
  python scripts/add_to_library.py --version RI_v8 --add <file> --new RI_v9          # dry run
  python scripts/add_to_library.py --version RI_v8 --add <file> --new RI_v9 --go
"""
import argparse
import json
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the version whose library comes first, e.g. RI_v8")
    parser.add_argument("--add", type=Path, required=True, help="the library to add (FASTA, name#class headers)")
    parser.add_argument("--new", required=True, help="the new version, e.g. RI_v9")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def names_of(path: Path):
    with c.open_text(path) as handle:
        return [line[1:].split()[0].split("#", 1)[0] for line in handle if line.startswith(">")]


def main(argv=None) -> None:
    args = parse_args(argv)
    if not args.add.exists():
        c.stop(f"{args.add} not found")
    _, _, base = c.library_files(args.version)
    target = c.version_folder(args.new)
    if target.exists():
        c.stop(f"{target} already exists; a version is never changed")
    work = args.scratch / args.new / "combine"
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    base_names, added_names = names_of(base), names_of(args.add)
    clashes = sorted(set(base_names) & set(added_names))
    dups = [n for n, k in Counter(added_names).items() if k > 1]
    if clashes or dups:
        c.stop("names would collide: " + ", ".join((clashes + dups)[:10]))
    print(f"{args.new} = {args.version} ({len(base_names):,} families) + {args.add.name} ({len(added_names):,} entries), "
          f"combined without clustering\nto: {target / 'library'}")
    if not args.go:
        print("dry run: nothing was written. Add --go to write the version.")
        return

    work.mkdir(parents=True)
    library = work / f"{args.new}-families.prefixed.fa"
    with open(library, "x") as out:
        for src in (base, args.add):
            with c.open_text(src) as handle:
                for line in handle:
                    out.write(line if line.endswith("\n") else line + "\n")
    fams = c.read_families(library)
    added_record = args.add.parent / "record.json"
    record = {"version": args.new, "kind": "combined library, no clustering", "base_version": args.version,
              "base_library": str(base), "base_sha256": c.sha256(base), "added_library": str(args.add),
              "added_sha256": c.sha256(args.add),
              "added_record": json.loads(added_record.read_text()) if added_record.exists() else None,
              "counts": {"base": len(base_names), "added": len(added_names), "families": len(fams),
                         "by_class": dict(Counter(la.group(f.label) for f in fams).most_common())},
              "written": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(work / "combine_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    for name in (library.name, "combine_record.json"):
        c.copy_new(work / name, target / "library" / name)
    print(f"\n{len(fams):,} families ({', '.join(f'{k} {v}' for k, v in record['counts']['by_class'].items())})")
    print(f"written: {target / 'library'}/\nnext: python scripts/mask_genomes.py --version {args.new} "
          f"--assemblies AC_fAstCal1_2 --go")


if __name__ == "__main__":
    main()
