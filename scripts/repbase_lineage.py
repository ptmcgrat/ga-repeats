#!/usr/bin/env python3
"""Repbase entries of a species' lineage, as RepeatMasker's -species chooses them (question 2).

Quah et al. 2025 combined their de novo libraries with "Repbase-derived RepeatMasker libraries". The usual way to derive
such a library, and the rule RepeatMasker itself applies with -species, is to keep the Repbase entries assigned to the
species or to one of its ancestral clades; entries of sister taxa (another cichlid species, zebrafish) are left out
(spec section 6, 2026-11-24). This applies that rule to our Repbase release (20181026, classified by repbase_match.py
into <scratch>/repbase/<release>/repbase.fa and repbase_entries.tsv):
  lineage   read from the printout of RepeatMasker 4.1.0's own taxonomy query, run read-only and saved to a file:
              perl $RMDIR/util/queryTaxonomyDatabase.pl -species "Astatotilapia calliptera" > <file>
            (the RepeatMasker in the GenomeSize environment; the taxonomy file it read is recorded by its sha256)
  matching  an entry's species or clade name (underscores read as spaces) is compared, ignoring case, with each line of
            the lineage and with each word of each line (a few lines carry two names, e.g. "Gnathostomata vertebrate")
  names     each selected entry is written as Repbase_<name>#<Type/SubType>, so that masking reports tell Repbase
            families from ours; characters other than letters, digits, '_', '-' and '.' become '_'
It prints the matched names with their entry counts, and the largest groups left out, for review.

Usage (from the repository root):
  python scripts/repbase_lineage.py --lineage-file <file> --taxonomy-sha256 <sha256 of taxonomy.dat>          # dry run
  python scripts/repbase_lineage.py --lineage-file <file> --taxonomy-sha256 <sha256 of taxonomy.dat> --go

Writes, never overwriting, <scratch>/repbase/<release>/lineage_<species>/: repbase_lineage.fa, selection.tsv (one row per
selected entry) and record.json.
"""
import argparse
import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import List, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lineage-file", type=Path, required=True, help="saved printout of queryTaxonomyDatabase.pl")
    parser.add_argument("--taxonomy-sha256", required=True, help="sha256 of the taxonomy.dat the query read")
    parser.add_argument("--release", default="20181026")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def parse_lineage(text: str) -> Tuple[str, List[str]]:
    """queryTaxonomyDatabase.pl's printout -> (species, lineage lines from the species up to the root)."""
    species, lineage, inside = "", [], False
    for line in text.splitlines():
        if line.startswith("Species ="):
            species = line.split("=", 1)[1].strip()
        elif line.startswith("Lineage ="):
            inside = True
            lineage.append(line.split("=", 1)[1].strip())
        elif inside:
            if not line.strip():
                break
            lineage.append(line.strip())
    if not species or not lineage:
        raise ValueError("no 'Species =' and 'Lineage =' lines: not a queryTaxonomyDatabase.pl printout")
    return species, lineage


def lineage_terms(lineage: List[str]) -> Set[str]:
    terms = set()
    for line in lineage:
        terms.add(line.lower())
        terms.update(w.lower() for w in line.split())
    return terms


def in_lineage(repbase_species: str, terms: Set[str]) -> bool:
    return repbase_species.replace("_", " ").strip().lower() in terms


def clean(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", name)


def main(argv=None) -> None:
    args = parse_args(argv)
    lib_dir = args.scratch / "repbase" / args.release
    repbase_fa, entries_file = lib_dir / "repbase.fa", lib_dir / "repbase_entries.tsv"
    for f in (args.lineage_file, repbase_fa, entries_file):
        if not f.exists():
            c.stop(f"{f} not found" + ("" if f == args.lineage_file else ": run repbase_match.py first"))
    try:
        species, lineage = parse_lineage(args.lineage_file.read_text())
    except ValueError as error:
        c.stop(f"{args.lineage_file}: {error}")
    slug = clean(species.replace(" ", "_"))
    out_dir = lib_dir / f"lineage_{slug}"
    if out_dir.exists():
        c.stop(f"{out_dir} exists; it is never overwritten")
    with open(entries_file, newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    terms = lineage_terms(lineage)
    chosen = [e for e in entries if in_lineage(e["species"], terms)]
    matched = Counter(e["species"] for e in chosen)
    left_out = Counter(e["species"] for e in entries if not in_lineage(e["species"], terms))
    print(f"{species}: lineage of {len(lineage)} levels ({lineage[-1]} ... {lineage[0]}); Repbase {args.release}: "
          f"{len(entries):,} entries")
    print(f"selected {len(chosen):,} entries under {len(matched)} names: "
          + ", ".join(f"{k} {v}" for k, v in matched.most_common()))
    print("largest groups left out: " + ", ".join(f"{k} {v}" for k, v in left_out.most_common(8)))
    print(f"to: {out_dir}")
    if not args.go:
        print("dry run: nothing was written. Add --go to write the subset.")
        return

    wanted = {e["name"] for e in chosen}
    out_dir.mkdir(parents=True)
    written, renamed, keep = {}, 0, False
    with open(repbase_fa) as src, open(out_dir / "repbase_lineage.fa", "x") as dst:
        for line in src:
            if line.startswith(">"):
                name, _, label = line[1:].strip().partition("#")
                keep = name in wanted and name not in written
                if keep:
                    new = "Repbase_" + clean(name)
                    renamed += new != "Repbase_" + name
                    written[name] = (new, label or "Unknown")
                    dst.write(f">{new}#{label or 'Unknown'}\n")
            elif keep:
                dst.write(line)
    by_entry = {e["name"]: e for e in chosen}
    c.write_tsv(out_dir / "selection.tsv", ["repbase_name", "library_name", "label", "species", "length"],
                [[n, new, lab, by_entry[n]["species"], by_entry[n].get("length", "")] for n, (new, lab) in written.items()])
    classes = Counter(la.group(lab) for _new, lab in written.values())
    record = {"species": species, "lineage": lineage, "rule": "Repbase entries whose species or clade is in the lineage "
              "(RepeatMasker's -species rule)", "repbase_release": args.release,
              "repbase_fa_sha256": c.sha256(repbase_fa), "lineage_file": str(args.lineage_file),
              "lineage_file_sha256": c.sha256(args.lineage_file), "taxonomy_dat_sha256": args.taxonomy_sha256,
              "taxonomy_source": "RepeatMasker 4.1.0 (GenomeSize environment), util/queryTaxonomyDatabase.pl, read-only",
              "entries": len(written), "names_matched": dict(matched.most_common()), "by_class": dict(classes.most_common()),
              "renamed_characters": renamed, "written": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(out_dir / "record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    missing = len(wanted) - len(written)
    print(f"\nwritten {len(written):,} entries ({', '.join(f'{k} {v}' for k, v in classes.most_common())})"
          + (f"; {missing} selected names were not found in {repbase_fa.name}" if missing else ""))
    print(f"subset: {out_dir / 'repbase_lineage.fa'}")
    print(f"next: python scripts/add_to_library.py --version RI_v8 --add {out_dir / 'repbase_lineage.fa'} --new RI_v9")


if __name__ == "__main__":
    main()
