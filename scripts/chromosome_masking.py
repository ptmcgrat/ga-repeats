#!/usr/bin/env python3
"""Masking per chromosome, LG3 in particular (question 2, briefing C2.7).

Kumar et al. 2025 describe chromosome 3 as highly repetitive and as a place where the M. zebra reference differs
from other assemblies; in the P1 download its length varies 1.8-fold across assemblies (58.7-105.3 Mb). This asks
whether LG3's masked share is high, and whether it tracks those length differences (a question, not a finding).

For each assembly masked by the version: every chromosome's length, its gaps (runs of N, counted in the genome), masked
bases (any annotation) and class shares (one label per base: the highest-scoring annotation); shares are given over the
sequence without gaps (repeats_025: DT_004f's LG3 is 27% gap, which made it look under-masked), and the raw share too.
Sequences not assigned to a chromosome are pooled as "unplaced"; the mitochondrion appears as MT where an assembly
lists it. Chromosome names come from NCBI's sequence_report.jsonl beside the genome (fields chrName, genbankAccession,
refseqAccession, length, role); the script stops if the report names none of the masked sequences. LG3 is the
chromosome named 3, LG3 or chr3.

Usage (from the repository root):
  python scripts/chromosome_masking.py --version RI_v5 --assemblies DT_004f DT_005m MZ_NMK2024a AC_fAstCal68          # dry run
  python scripts/chromosome_masking.py --version RI_v5 --assemblies DT_004f DT_005m MZ_NMK2024a AC_fAstCal68 --go

Writes, never overwriting, <scratch>/comparisons/chromosomes/<version>_chromosomes.tsv (one row per assembly and
chromosome).
"""
import argparse
import gzip
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

LG3 = re.compile(r"^(LG|chr)?0*3$", re.IGNORECASE)
UNPLACED = "unplaced"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a version that has masked the assemblies, e.g. RI_v5")
    parser.add_argument("--assemblies", nargs="+", required=True, help="assemblies masked by it")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def chromosome_names(lines: Iterable[str]) -> Tuple[Dict[str, Tuple[str, int]], Dict[str, Tuple[str, int]]]:
    """From sequence_report.jsonl: (GenBank accession -> (chromosome or 'unplaced', length), the same by RefSeq accession)."""
    genbank, refseq = {}, {}
    for line in lines:
        if not line.strip():
            continue
        d = json.loads(line)
        name = str(d.get("chrName") or "")
        if not name or name.lower() in ("un", "na") or d.get("role", "assembled-molecule") != "assembled-molecule":
            name = UNPLACED
        value = (name, int(d.get("length", 0)))
        if d.get("genbankAccession"):
            genbank[d["genbankAccession"]] = value
        if d.get("refseqAccession"):
            refseq[d["refseqAccession"]] = value
    return genbank, refseq


def gap_counts(lines: Iterable[str]) -> Counter:
    """Sequence name -> number of N (either case), from FASTA lines."""
    gaps: Counter = Counter()
    name = None
    for line in lines:
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif name is not None:
            gaps[name] += line.count("N") + line.count("n")
    return gaps


def per_chromosome(names: Dict[str, Tuple[str, int]], raw: Iterable[la.Scored], rows: Iterable[la.Row],
                   gaps: Counter = None) -> Dict[str, Counter]:
    """names: each genome sequence once -> (chromosome, length). Per chromosome: length, gaps, masked (any annotation)
    and bases per class group (one label per base)."""
    out: Dict[str, Counter] = defaultdict(Counter)
    for seq, (chrom, length) in names.items():
        out[chrom]["length"] += length
        out[chrom]["gaps"] += (gaps or {}).get(seq, 0)
    for seq, spans in c.merge_by_sequence((r[0], r[1], r[2]) for r in raw).items():
        out[names.get(seq, (UNPLACED, 0))[0]]["masked"] += sum(b - a for a, b in spans)
    for seq, a, b, _family, label in rows:
        out[names.get(seq, (UNPLACED, 0))[0]][la.group(label)] += b - a
    return out


def chromosome_order(name: str):
    digits = re.sub(r"\D", "", name)
    return (name == UNPLACED, int(digits) if digits else 999, name)


def main(argv=None) -> None:
    args = parse_args(argv)
    plan = []
    for assembly in args.assemblies:
        genome = c.genome_path(assembly, c.panel_accession(assembly))
        report = genome.parent / "sequence_report.jsonl"
        bed = c.version_folder(args.version) / "masked" / assembly / f"{assembly}.bed"
        for f, what in ((report, "sequence report"), (bed, f"{args.version}'s masking")):
            if not f.exists():
                c.stop(f"{assembly}: no {what} ({f})")
        plan.append((assembly, report, bed))
    out_dir = args.scratch / "comparisons" / "chromosomes"
    target = out_dir / f"{args.version}_chromosomes.tsv"
    if target.exists():
        with open(target) as handle:
            old = "gap_bp" not in handle.readline().split("\t")
        c.stop(f"{target} exists" + (f" from the earlier version (without gaps); remove it with: rm {target}" if old
                                     else "; results are never overwritten"))
    print(f"{args.version}: masking per chromosome in " + ", ".join(a for a, _, _ in plan) + f"\nresult: {target}")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    table, lg3 = [], []
    for assembly, report, bed in plan:
        with open(report) as handle:
            genbank, refseq = chromosome_names(handle)
        raw = la.read_labelled_bed(bed)
        masked_seqs = {r[0] for r in raw}
        names = max((genbank, refseq), key=lambda m: len(masked_seqs & set(m)))   # the accessions this genome uses
        if not masked_seqs & set(names):
            c.stop(f"{assembly}: the sequence report names none of the masked sequences; check its fields ({report})")
        print(f"{assembly}: counting gaps ...", flush=True)
        with gzip.open(c.genome_path(assembly, c.panel_accession(assembly)), "rt") as handle:
            gaps = gap_counts(handle)
        stats = per_chromosome(names, raw, la.one_label_per_base(raw), gaps)
        total_len = sum(s["length"] - s["gaps"] for s in stats.values())
        total_masked = sum(s["masked"] for s in stats.values())
        placed = {}
        for chrom in sorted(stats, key=chromosome_order):
            s = stats[chrom]
            real = s["length"] - s["gaps"]
            if real <= 0:
                continue
            share = 100 * s["masked"] / real
            if chrom != UNPLACED:
                placed[chrom] = share
            table.append([assembly, chrom, s["length"], s["gaps"], s["masked"], f"{100 * s['masked'] / s['length']:.2f}",
                          f"{share:.2f}"] + [f"{100 * s[g] / real:.2f}" for g in la.GROUPS])
        three = next((k for k in placed if LG3.match(k)), None)
        gap_total = sum(s["gaps"] for s in stats.values())
        line = (f"{assembly}: gaps {100 * gap_total / (total_len + gap_total):.2f}%; genome {100 * total_masked / total_len:.2f}% "
                f"of non-gap sequence masked; chromosomes "
                f"{min(placed.values()):.1f}-{max(placed.values()):.1f}% (median {statistics.median(placed.values()):.1f}%)"
                if placed else f"{assembly}: no chromosomes named in the report")
        if three:
            rank = sorted(placed.values(), reverse=True).index(placed[three]) + 1
            real3 = stats[three]["length"] - stats[three]["gaps"]
            line += (f"; LG3 ({three}) {stats[three]['length'] / 1e6:.1f} Mb ({stats[three]['gaps'] / 1e6:.1f} Mb gap), "
                     f"{placed[three]:.1f}% of non-gap masked, rank {rank} of {len(placed)}")
            lg3.append((assembly, real3, placed[three], stats[three]))
        print(line)
    if lg3:
        print("\nLG3 across assemblies (non-gap length, masked share and class shares of the non-gap sequence):")
        for assembly, length, p, s in sorted(lg3, key=lambda x: x[1]):
            print(f"  {assembly:<14} {length / 1e6:6.1f} Mb  {p:5.1f}%  " + ", ".join(f"{g} {100 * s[g] / length:.1f}" for g in la.GROUPS))
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["assembly", "chromosome", "length", "gap_bp", "masked_bp", "masked_percent_of_length",
                         "masked_percent_of_non_gap"] + [f"{g}_percent_of_non_gap" for g in la.GROUPS], table)
    print(f"table: {target}")


if __name__ == "__main__":
    main()
