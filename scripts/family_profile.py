#!/usr/bin/env python3
"""Profile one repeat family: is it dispersed copies or tandem arrays, where does it sit, what does its consensus hold?

Written for the Unknown family behind the spike at 3-4% divergence in RI_v5's landscape of DT_004f
(MZ_MZ4f_rnd-5_family-9801: 12.9 Mb in 102 annotations, longer than its 13.4 kb consensus, as RepeatMasker reports joined
copies of a tandemly repeated consensus); it works for any family.

For each assembly (masked by the version):
  annotations  number, total bases, lengths (median, largest), divergence (median), and how many are longer than the
               consensus (joined tandem copies)
  arrays       annotations of the family no more than --join bases apart (default 1,000), merged: how many, their
               lengths, the share of the family's bases in arrays of 10 kb or more, and the largest arrays with position
  chromosomes  bases per chromosome (named from sequence_report.jsonl) and where along each chromosome (fraction of its
               length) the bases sit
and once for the family:
  consensus    length and label; the consensus aligned to itself (blastn in the container): off-diagonal matches mean
               an internally repeated unit, and their most common offset is the period
  Repbase      its best Repbase match, from repbase_match.py's table, if there is one

Usage (from the repository root):
  python scripts/family_profile.py --version RI_v5 --family MZ_MZ4f_rnd-5_family-9801 --assemblies DT_004f MZ_MZ4f          # dry run
  python scripts/family_profile.py --version RI_v5 --family MZ_MZ4f_rnd-5_family-9801 --assemblies DT_004f MZ_MZ4f --go

Writes, never overwriting, <scratch>/comparisons/family_profile/<version>_<family>_arrays.tsv (every array).
"""
import argparse
import csv
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

Ann = Tuple[str, int, int, float]   # sequence, start, end (0-based half-open), divergence


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the library version, e.g. RI_v5")
    parser.add_argument("--family", required=True, help="the family name without its #class, e.g. MZ_MZ4f_rnd-5_family-9801")
    parser.add_argument("--assemblies", nargs="+", required=True, help="assemblies masked by the version")
    parser.add_argument("--join", type=int, default=1000, help="annotations at most this far apart form one array (default 1000)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def family_annotations(lines, family: str) -> List[Ann]:
    out = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[3] == family:
            out.append((f[0], int(f[1]), int(f[2]), float(f[7])))
    return sorted(out)


def arrays(anns: List[Ann], join: int) -> List[Tuple[str, int, int, int, int]]:
    """(sequence, start, end, annotations, bases annotated) for runs of annotations at most join bases apart."""
    out: List[list] = []
    for seq, a, b, _d in anns:
        if out and out[-1][0] == seq and a - out[-1][2] <= join:
            out[-1][2] = max(out[-1][2], b)
            out[-1][3] += 1
            out[-1][4] += b - a
        else:
            out.append([seq, a, b, 1, b - a])
    return [tuple(x) for x in out]


def period(self_hits: List[Tuple[int, int, int, int]], length: int) -> Tuple[int, int]:
    """From the consensus aligned to itself ((qstart, qend, sstart, send), 1-based): the most common offset between
    off-diagonal matches (0 if none) and the consensus bases they cover."""
    offsets: Counter = Counter()
    covered = set()
    for qs, qe, ss, se in self_hits:
        if (qs, qe) == (ss, se) or (qe - qs + 1) >= 0.95 * length:
            continue
        offsets[abs(qs - ss)] += qe - qs + 1
        covered.update(range(qs, qe + 1))
    return (offsets.most_common(1)[0][0] if offsets else 0), len(covered)


def main(argv=None) -> None:
    args = parse_args(argv)
    _, _, library = c.library_files(args.version)
    families = {f.family: f for f in c.read_families(library)}
    if args.family not in families:
        c.stop(f"{args.family} is not a family of {args.version} ({library})")
    fam = families[args.family]
    beds = {}
    for assembly in args.assemblies:
        beds[assembly] = c.version_folder(args.version) / "masked" / assembly / f"{assembly}.bed"
        if not beds[assembly].exists():
            c.stop(f"{args.version} has not masked {assembly}: {beds[assembly]} not found")
    out_dir = args.scratch / "comparisons" / "family_profile"
    target = out_dir / f"{args.version}_{args.family}_arrays.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    work = args.scratch / args.version / "family_profile" / args.family
    print(f"{args.family} ({fam.label}, consensus {fam.length:,} bp) of {args.version} in " + ", ".join(args.assemblies))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    # the consensus against itself
    work.mkdir(parents=True, exist_ok=True)
    consensus = work / "consensus.fa"
    seq, current = [], None
    with open(library) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split("#", 1)[0].split()[0]
            elif current == args.family:
                seq.append(line.strip())
    consensus.write_text(f">{args.family}\n{''.join(seq)}\n")
    self_out = work / "self_hits.tsv"
    c.run(c.container(args.sif, "blastn", "-query", str(consensus), "-subject", str(consensus), "-evalue", "1e-10",
                      "-outfmt", "6 qstart qend sstart send", "-out", str(self_out)), True, work / "blastn.log")
    self_hits = [tuple(int(x) for x in line.split("\t")[:4]) for line in open(self_out) if line.strip()]
    per, covered = period(self_hits, fam.length)
    print(f"\nconsensus {fam.length:,} bp, {fam.label}: "
          + (f"internally repeated: off-diagonal matches cover {covered:,} bp, most common offset {per:,} bp"
             if per else "no internal repeat found"))
    repbase = args.scratch / "comparisons" / "repbase" / f"{args.version}_best_match.tsv"
    if repbase.exists():
        with open(repbase) as handle:
            row = next((r for r in csv.DictReader(handle, delimiter="\t") if r["family"].split("#")[0] == args.family), None)
        print("best Repbase match: " + (f"{row.get('repbase_entry')} ({row.get('repbase_label')}), covering "
                                         f"{row.get('coverage_of_family')} of the family" if row and row.get("repbase_entry") else "none"))

    table = []
    for assembly in args.assemblies:
        genome = c.genome_path(assembly, c.panel_accession(assembly))
        names, lengths = {}, {}
        report = genome.parent / "sequence_report.jsonl"
        if report.exists():
            for line in open(report):
                if line.strip():
                    d = json.loads(line)
                    for key in ("genbankAccession", "refseqAccession"):
                        if d.get(key):
                            names[d[key]] = str(d.get("chrName") or d[key])
                            lengths[d[key]] = int(d.get("length", 0))
        with open(beds[assembly]) as handle:
            anns = family_annotations(handle, args.family)
        if not anns:
            print(f"\n{assembly}: no annotations of {args.family}")
            continue
        sizes = [b - a for _, a, b, _ in anns]
        total = sum(sizes)
        arr = arrays(anns, args.join)
        in_big = sum(x[4] for x in arr if x[2] - x[1] >= 10000)
        print(f"\n{assembly}: {len(anns):,} annotations, {total:,} bp; length median {statistics.median(sizes):,.0f} bp, largest "
              f"{max(sizes):,} bp; {sum(1 for s in sizes if s > fam.length):,} longer than the consensus; divergence median "
              f"{statistics.median(d for *_, d in anns):.1f}%")
        print(f"  arrays (annotations at most {args.join:,} bp apart): {len(arr):,}; {100 * in_big / total:.0f}% of the bases in "
              f"arrays of 10 kb or more; largest: " + "; ".join(
                  f"{names.get(s, s)} {a / 1e6:.2f}-{b / 1e6:.2f} Mb ({(b - a) / 1e3:.0f} kb"
                  + (f", at {100 * (a + b) / 2 / lengths[s]:.0f}% of the chromosome" if lengths.get(s) else "") + ")"
                  for s, a, b, _n, _bp in sorted(arr, key=lambda x: -(x[2] - x[1]))[:5]))
        by_chrom: Dict[str, int] = defaultdict(int)
        for s, a, b, _ in anns:
            by_chrom[names.get(s, s)] += b - a
        print("  by chromosome: " + ", ".join(f"{k} {v / 1e6:.2f} Mb" for k, v in sorted(by_chrom.items(), key=lambda kv: -kv[1])[:8]))
        for s, a, b, n, bp in arr:
            table.append([assembly, names.get(s, s), s, a, b, b - a, n, bp,
                          f"{100 * (a + b) / 2 / lengths[s]:.1f}" if lengths.get(s) else ""])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["assembly", "chromosome", "seq", "start", "end", "length", "annotations", "bp_annotated",
                         "position_percent_of_chromosome"], table)
    print(f"table: {target}")


if __name__ == "__main__":
    main()
