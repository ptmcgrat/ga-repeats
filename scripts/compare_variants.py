#!/usr/bin/env python3
"""Compare library variants side by side: raw, curated by a tool, corrected by our rules (question 2).

For each version in --versions (e.g. RI_v1 RI_v9 RI_v10) it gathers what the checks have already written and prints
one table, a column per version:
  run          the tool and its version, and run time (curation_record.json of an imported version)
  library      families; share labelled Unknown; median consensus length; families by class group
  masking      share of --assembly masked, and share of the masked bases labelled Unknown   (mask_genomes.py)
  curated      of the bases the curated library masks: the share the version masks too, and the share it gives the
               curated class, overall and per curated class                               (curated_annotation.py)
  recovery     curated families perfect / good / present / not found                     (family_recovery.py)
  composites   relabel / composite / unknown / partial families, and the composites' footprint (composite_families.py)
  Ensembl      for a reference genome: of Ensembl's masked bases, the share the version masks, and the reverse
               (ensembl_masking.py)
A check not yet run for a version shows as '-', and the printout names the command that fills it.

Usage (from the repository root):
  python scripts/compare_variants.py --versions RI_v1 RI_v9 RI_v10          # dry run: lists what is there
  python scripts/compare_variants.py --versions RI_v1 RI_v9 RI_v10 --go

Writes, never overwriting, <scratch>/comparisons/variants/<versions joined by '_'>_<assembly>.tsv.
"""
import argparse
import csv
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

CLASSES = ["DNA", "LINE", "LTR", "SINE", "Other", "Unknown"]
TOOL_NAMES = {"mchelper": "MCHelper", "tetrimmer": "TEtrimmer"}
CURATED = ["DNA", "LINE", "LTR", "SINE", "Other"]


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--versions", nargs="+", required=True)
    parser.add_argument("--assembly", default="AC_fAstCal1_2",
                        help="genome all variants are masked on (default fAstCal1.2, spec section 6, 2026-11-23)")
    parser.add_argument("--tag", default="cichlid", help="the curated set (family_recovery.py, curated_annotation.py)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def pct(n: float, d: float) -> Optional[float]:
    return 100.0 * n / d if d else None


def read_tsv(path: Path) -> List[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def curated_totals(bed: Path) -> Counter:
    totals: Counter = Counter()
    for _s, a, b, _f, label in la.one_label_per_base(la.read_labelled_bed(bed)):
        totals[la.group(label)] += b - a
    return totals


def metrics(version: str, args, genome_bp: int, cur_totals: Optional[Counter], missing: List[str]) -> Dict[str, object]:
    m: Dict[str, object] = {}
    lib_dir = c.version_folder(version) / "library"
    rec = lib_dir / "curation_record.json"
    if rec.exists():
        with open(rec) as handle:
            r = json.load(handle)
        number = re.search(r"\d+(\.\d+)+", r["tool_version"])
        m["run"] = f"{TOOL_NAMES.get(r['tool'], r['tool'])} {number.group(0) if number else ''}".strip()
        m["source version"] = r["source_version"]
        m["run hours"] = r["run"].get("hours")
    else:
        m["run"], m["source version"], m["run hours"] = "as built", "-", None
    _, _, library = c.library_files(version)
    fams = c.read_families(library)
    groups = Counter(la.group(f.label) for f in fams)
    m["families"] = len(fams)
    m["families Unknown %"] = pct(groups["Unknown"], len(fams))
    m["median length bp"] = int(statistics.median(f.length for f in fams)) if fams else None
    for g in CLASSES:
        m[f"families {g}"] = groups.get(g, 0)
    bed = c.version_folder(version) / "masked" / args.assembly / f"{args.assembly}.bed"
    if bed.exists():
        rows = la.one_label_per_base(la.read_labelled_bed(bed))
        masked = sum(b - a for _s, a, b, _f, _l in rows)
        unknown = sum(b - a for _s, a, b, _f, l in rows if la.group(l) == "Unknown")
        m[f"masked % of {args.assembly}"] = pct(masked, genome_bp)
        m["masked bases Unknown %"] = pct(unknown, masked)
    else:
        missing.append(f"python scripts/mask_genomes.py --version {version} --assemblies {args.assembly} --go")
    classes = args.scratch / "comparisons" / "curated_annotation" / f"{version}_{args.tag}_{args.assembly}_classes.tsv"
    if classes.exists() and cur_totals:
        table = Counter()
        for r in read_tsv(classes):
            table[(r["curated_class"], r["version_class"])] += int(r["bp"])
        both = sum(table.values())
        m["curated masked by version %"] = pct(both, sum(cur_totals.values()))
        m["curated same class %"] = pct(sum(v for (a, b), v in table.items() if a == b), both)
        for g in CURATED:
            row_total = sum(v for (a, _b), v in table.items() if a == g)
            m[f"curated {g} same class %"] = pct(table.get((g, g), 0), row_total)
    else:
        missing.append(f"python scripts/curated_annotation.py --version {version} --assembly {args.assembly} --go")
    summary = args.scratch / "comparisons" / "family_recovery" / f"{version}_{args.tag}_summary.tsv"
    if summary.exists():
        all_row = next((r for r in read_tsv(summary) if r["type"] == "ALL"), None)
        if all_row:
            m["recovery perfect/good/present/not found"] = "/".join(all_row[k] for k in ("perfect", "good", "present", "not found"))
    else:
        missing.append(f"python scripts/family_recovery.py --version {version} --go")
    ens = args.scratch / "comparisons" / "ensembl" / f"{version}_{args.assembly}.tsv"
    if ens.exists():
        values = {r["metric"]: r["value"] for r in read_tsv(ens)}
        for k in ("Ensembl-masked bases masked by version %", "version-masked bases masked by Ensembl %"):
            if k in values:
                m[k] = float(values[k])
    elif c.is_reference(args.assembly):
        missing.append(f"python scripts/ensembl_masking.py --version {version} --assembly {args.assembly} --go")
    comp = args.scratch / "comparisons" / "composite" / f"{version}_{args.tag}_families.tsv"
    if comp.exists():
        rows = read_tsv(comp)
        status = Counter(r["status"] for r in rows)
        for s in ("relabel", "composite", "unknown", "partial"):
            m[f"{s} families"] = status.get(s, 0)
        m["composite footprint Mb"] = sum(int(r["footprint_bp"]) for r in rows if r["status"] == "composite") / 1e6
    else:
        missing.append(f"python scripts/composite_families.py --version {version} --go (after family_recovery.py)")
    return m


def fmt(v: object) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.2f}"
    return str(v)


def main(argv=None) -> None:
    args = parse_args(argv)
    for v in args.versions:
        if not c.version_folder(v).exists():
            c.stop(f"{c.version_folder(v)} not found")
    out_dir = args.scratch / "comparisons" / "variants"
    target = out_dir / f"{'_'.join(args.versions)}_{args.assembly}.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    print(f"comparing {', '.join(args.versions)} on {args.assembly}; table: {target}")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to compare.")
        return
    genome_bp = sum(c.fasta_lengths(c.genome_path(args.assembly, c.panel_accession(args.assembly))).values())
    cur_bed = args.scratch / "curated_annotation" / args.tag / args.assembly / f"{args.assembly}.curated.bed"
    cur_totals = curated_totals(cur_bed) if cur_bed.exists() else None
    missing: List[str] = []
    if cur_totals is None:
        missing.append(f"python scripts/curated_annotation.py --version <any> --assembly {args.assembly} --go "
                       f"(it makes the curated masking)")
    results = {v: metrics(v, args, genome_bp, cur_totals, missing) for v in args.versions}
    keys: List[str] = []
    for m in results.values():
        keys += [k for k in m if k not in keys]
    width = max(len(k) for k in keys)
    print(f"\n{'':{width}s}  " + "  ".join(f"{v:>16s}" for v in args.versions))
    for k in keys:
        print(f"{k:{width}s}  " + "  ".join(f"{fmt(results[v].get(k))[:16]:>16s}" for v in args.versions))
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["metric"] + args.versions, [[k] + [fmt(results[v].get(k)) for v in args.versions] for k in keys])
    if missing:
        print("\nnot yet run ('-' above):")
        for line in dict.fromkeys(missing):
            print(f"  {line}")
    print(f"\ntable: {target}")


if __name__ == "__main__":
    main()
