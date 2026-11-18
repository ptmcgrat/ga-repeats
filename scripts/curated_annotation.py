#!/usr/bin/env python3
"""Annotation against the curated cichlid library, genome-wide (question 2).

The known insertions test three sites; this tests every copy of the curated families in a genome. It masks one genome
(default DT_004f) with the curated cichlid Repbase entries alone (the Repbase entries of --species, default Cichlidae and
Oreochromis_niloticus: 165 entries, 70% DNA; taken from the classified Repbase repbase_match.py built), with the same
RepeatMasker settings as our maskings (-engine rmblast -no_is), then compares base by base (one label per base: the
highest-scoring annotation) with the version's masking:
  sensitivity  of the bases the curated library masks, the share the version masks too, per curated class
  class        of the bases both mask, the share the version gives the same class group, per curated class, and the
               classes it gives instead
  families     the version's families that most often give curated bases another class, with their origin (LTR pipeline
               or RepeatScout/RECON) and, where composite_families.py has run, its verdict on them
The curated library is partial, so bases it does not mask say nothing about the version.

The curated masking does not depend on the version: it is kept in <scratch>/curated_annotation/<tag>/<assembly>/ and
reused when it is complete.

Usage (from the repository root):
  python scripts/curated_annotation.py --version RI_v5 --assembly DT_004f          # dry run
  python scripts/curated_annotation.py --version RI_v5 --assembly DT_004f --go

Writes, never overwriting, in <scratch>/comparisons/curated_annotation/: <version>_<tag>_<assembly>_classes.tsv
(curated class by the version's class, in bases) and <version>_<tag>_<assembly>_families.tsv (the version's families on
curated bases of another class).
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the version to compare, e.g. RI_v5")
    parser.add_argument("--assembly", default="DT_004f", help="a genome the version has masked (default DT_004f)")
    parser.add_argument("--species", nargs="+", default=["Cichlidae", "Oreochromis_niloticus"])
    parser.add_argument("--tag", default="cichlid")
    parser.add_argument("--release", default="20181026")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def write_curated(repbase_fa: Path, wanted: set, dst: Path) -> Tuple[Dict[str, int], List[str]]:
    """Write the wanted entries (first copy of each name) to dst; return their lengths and any duplicated names."""
    lengths: Dict[str, int] = {}
    duplicates, keep, name = [], False, None
    with open(repbase_fa) as src, open(dst, "x") as out:
        for line in src:
            if line.startswith(">"):
                name = line[1:].split("#", 1)[0].strip()
                keep = name in wanted and name not in lengths
                if name in wanted and name in lengths:
                    duplicates.append(name)
                if keep:
                    lengths[name] = 0
            elif keep:
                lengths[name] += len(line.strip())
            if keep:
                out.write(line)
    return lengths, duplicates


def overlaps(a: List[la.Row], b: List[la.Row]) -> Iterator[Tuple[la.Row, la.Row, int]]:
    """Pairs of (one label per base) rows of a and b that overlap, with the overlap's length."""
    by_a: Dict[str, list] = defaultdict(list)
    by_b: Dict[str, list] = defaultdict(list)
    for r in a:
        by_a[r[0]].append(r)
    for r in b:
        by_b[r[0]].append(r)
    for seq, rows_a in by_a.items():
        rows_b = sorted(by_b.get(seq, []), key=lambda r: r[1])
        j = 0
        for ra in sorted(rows_a, key=lambda r: r[1]):
            while j < len(rows_b) and rows_b[j][2] <= ra[1]:
                j += 1
            k = j
            while k < len(rows_b) and rows_b[k][1] < ra[2]:
                yield ra, rows_b[k], min(ra[2], rows_b[k][2]) - max(ra[1], rows_b[k][1])
                k += 1


def main(argv=None) -> None:
    args = parse_args(argv)
    lib_dir = args.scratch / "repbase" / args.release
    repbase_fa, entries_file = lib_dir / "repbase.fa", lib_dir / "repbase_entries.tsv"
    for f in (repbase_fa, entries_file):
        if not f.exists():
            c.stop(f"{f} not found: run repbase_match.py first (it builds the classified Repbase)")
    ours_bed = c.version_folder(args.version) / "masked" / args.assembly / f"{args.assembly}.bed"
    if not ours_bed.exists():
        c.stop(f"{args.version} has not masked {args.assembly}: {ours_bed} not found")
    work = args.scratch / "curated_annotation" / args.tag / args.assembly
    curated_bed = work / f"{args.assembly}.curated.bed"
    out_dir = args.scratch / "comparisons" / "curated_annotation"
    stem = f"{args.version}_{args.tag}_{args.assembly}"
    targets = [out_dir / f"{stem}_classes.tsv", out_dir / f"{stem}_families.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    reuse = curated_bed.exists()
    if work.exists() and not reuse:
        c.stop(f"{work} exists from an unfinished attempt; clear it with: rm -r {work}")
    with open(entries_file) as handle:
        wanted = {r["name"] for r in csv.DictReader(handle, delimiter="\t") if r["species"] in set(args.species)}
    pa = max(1, args.threads // c.RMBLAST_CORES_PER_JOB)
    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    fasta = work / f"{args.assembly}.fna"
    command = c.container(args.sif, "RepeatMasker", "-lib", str(work / "curated.fa"), "-engine", "rmblast", "-pa", str(pa),
                          "-no_is", "-dir", str(work), str(fasta))
    print(f"{args.assembly}: masked with the {len(wanted)} curated entries of {', '.join(args.species)}, then compared with "
          f"{args.version}\ncurated masking: {curated_bed}" + (" (exists, reused)" if reuse else ""))
    if not reuse:
        c.run(command, False, work / "repeatmasker.log")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    if not reuse:
        work.mkdir(parents=True)
        started = datetime.now()
        lengths, duplicates = write_curated(repbase_fa, wanted, work / "curated.fa")
        if duplicates:
            print(f"note: {len(duplicates)} entry names occur more than once in {repbase_fa.name}; the first copy is used: "
                  + ", ".join(sorted(set(duplicates))[:10]))
        c.decompress(genome, fasta)
        seq_lengths = c.fasta_lengths(fasta)
        c.run(command, True, work / "repeatmasker.log")
        out_file = work / f"{args.assembly}.fna.out"
        if not out_file.exists():
            c.stop(f"RepeatMasker wrote no {out_file.name}; see {work / 'repeatmasker.log'}")
        bed_rows = [c.hit_to_bed(h) for h in c.parse_out(out_file)]
        problems = c.check_bed(bed_rows, seq_lengths)
        if problems:
            c.stop("BED check failed:\n  " + "\n  ".join(problems))
        c.write_tsv(curated_bed, ["#chrom", "start", "end", "family", "score", "strand", "class_family", "divergence",
                                  "sw_score", "rm_id"], bed_rows)
        fasta.unlink()
        with open(work / "record.json", "x") as handle:
            json.dump({"assembly": args.assembly, "species": args.species, "entries": len(lengths),
                       "command": c.shell_text(command), "started": f"{started:%Y-%m-%d %H:%M}",
                       "finished": f"{datetime.now():%Y-%m-%d %H:%M}"}, handle, indent=2)
            handle.write("\n")

    print("reading the two maskings ...", flush=True)
    genome_bp = sum(c.fasta_lengths(genome).values())
    cur = la.one_label_per_base(la.read_labelled_bed(curated_bed))
    ours = la.one_label_per_base(la.read_labelled_bed(ours_bed))
    cur_total: Counter = Counter()
    for _s, a, b, _f, label in cur:
        cur_total[la.group(label)] += b - a
    table: Counter = Counter()          # (curated class, our class) -> bases
    wrong: Dict[Tuple[str, str], Counter] = defaultdict(Counter)   # (our family, our label) -> {curated class: bases}
    for ra, rb, n in overlaps(cur, ours):
        cg, og = la.group(ra[4]), la.group(rb[4])
        table[(cg, og)] += n
        if og != cg:
            wrong[(rb[3], rb[4])][cg] += n
    composite = {}
    comp_file = args.scratch / "comparisons" / "composite" / f"{args.version}_{args.tag}_families.tsv"
    if comp_file.exists():
        with open(comp_file) as handle:
            composite = {r["family"]: r["status"] for r in csv.DictReader(handle, delimiter="\t")}
    rows = [[cg, og, n] for (cg, og), n in sorted(table.items(), key=lambda kv: -kv[1])]
    fam_rows = []
    for (fam, label), by in sorted(wrong.items(), key=lambda kv: -sum(kv[1].values())):
        fam_rows.append([fam, label, "LTR pipeline" if "_ltr-1_" in fam else "RepeatScout/RECON", sum(by.values()),
                         ";".join(f"{k} {v}" for k, v in by.most_common()), composite.get(fam, "not run" if not composite else "")])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["curated_class", "version_class", "bp"], rows)
    c.write_tsv(targets[1], ["family", "label", "origin", "bp_on_curated_other_class", "curated_classes", "composite_status"], fam_rows)

    total = sum(cur_total.values())
    print(f"\n{args.assembly} ({genome_bp:,} bp): the curated library masks {100 * total / genome_bp:.2f}% "
          f"({', '.join(f'{k} {100 * v / genome_bp:.2f}%' for k, v in cur_total.most_common())})")
    print(f"{'curated class':14s} {'curated bp':>12s} {'masked by us':>13s} {'same class':>11s}   our class instead")
    for cg, n in cur_total.most_common():
        both = sum(v for (a, b), v in table.items() if a == cg and b != la.UNMASKED)
        same = table.get((cg, cg), 0)
        instead = Counter({b: v for (a, b), v in table.items() if a == cg and b not in (cg, la.UNMASKED)})
        print(f"{cg:14s} {n:>12,} {100 * both / n:>12.1f}% {100 * same / both if both else 0:>10.1f}%   "
              + ", ".join(f"{k} {100 * v / both:.1f}%" for k, v in instead.most_common(3)))
    print("families giving curated bases another class (largest first):")
    for r in fam_rows[:12]:
        print(f"  {r[0]:30s} {r[1]:20s} {r[2]:18s} {r[3] / 1e3:8,.0f} kb on {r[4]}" + (f"  [{r[5]}]" if r[5] else ""))
    lp = sum(r[3] for r in fam_rows if r[2] == "LTR pipeline")
    allw = sum(r[3] for r in fam_rows)
    if allw:
        print(f"of the curated bases given another class, {100 * lp / allw:.0f}% are given it by LTR-pipeline families")
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
