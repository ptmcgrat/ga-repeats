#!/usr/bin/env python3
"""Which bases a library masks, against Ensembl's published masking of the same assembly (question 2).

Ensembl's soft-masked genome (fetch_ensembl.py; lower-case bases are those Ensembl's repeat pipeline masked) is read as
intervals. Its sequence names (chromosome numbers, GenBank accessions) are mapped to NCBI's RefSeq names through NCBI's
assembly report (assembly_report.txt: sequence name, assigned chromosome for chromosomes only, GenBank and RefSeq
accessions), and every mapping is confirmed by sequence length; a sequence that cannot be mapped so is left out and
listed. Against the version's masking of the same assembly (one label per base) it reports the bases masked by Ensembl,
by the version, by both and by only one, and, for each of the version's class groups, the share of its bases Ensembl
also masks. Ensembl's masking carries no classes: this measures over- and under-masking, not labels.

Usage (from the repository root):
  python scripts/ensembl_masking.py --version RI_v9 --assembly AC_fAstCal1_2          # dry run
  python scripts/ensembl_masking.py --version RI_v9 --assembly AC_fAstCal1_2 --go

Writes, never overwriting, <scratch>/comparisons/ensembl/<version>_<assembly>.tsv (metric, value) and
<version>_<assembly>_classes.tsv.
"""
import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import fetch_ensembl as fe  # noqa: E402
import label_agreement as la  # noqa: E402

LOWER = re.compile(r"[a-z]+")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--assembly", default="AC_fAstCal1_2")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def report_aliases(lines) -> Dict[str, str]:
    """NCBI assembly report -> every name a sequence goes by -> its RefSeq accession. The assigned chromosome name
    is used only for the chromosome itself (role assembled-molecule), never for scaffolds placed on it."""
    aliases: Dict[str, str] = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 7 or f[6] in ("", "na"):
            continue
        refseq = f[6]
        names = [f[0], f[4], f[6]] + ([f[2]] if f[1] == "assembled-molecule" else [])
        for n in names:
            if n and n != "na":
                for k in {n, n.split(".")[0], re.sub(r"^(LG|chr)", "", n)}:
                    aliases.setdefault(k, refseq)
    return aliases


def softmasked_intervals(lines) -> Iterator[Tuple[str, int, List[Tuple[int, int]]]]:
    """A soft-masked FASTA -> (name, length, lower-case intervals, 0-based half-open), one record at a time."""
    name, pos, runs = None, 0, []
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                yield name, pos, runs
            name, pos, runs = line[1:].split()[0], 0, []
            continue
        seq = line.strip()
        for m in LOWER.finditer(seq):
            a, b = pos + m.start(), pos + m.end()
            if runs and runs[-1][1] == a:
                runs[-1] = (runs[-1][0], b)
            else:
                runs.append((a, b))
        pos += len(seq)
    if name is not None:
        yield name, pos, runs


def map_names(ensembl_lengths: Dict[str, int], aliases: Dict[str, str], ncbi_lengths: Dict[str, int]) -> Tuple[Dict[str, str], List[str]]:
    """Ensembl name -> RefSeq name, by alias then by unique length; every mapping confirmed by length."""
    mapping, used = {}, set()
    for name, length in ensembl_lengths.items():
        for key in (name, name.split(".")[0]):
            target = aliases.get(key)
            if target and ncbi_lengths.get(target) == length and target not in used:
                mapping[name] = target
                used.add(target)
                break
    by_length: Dict[int, List[str]] = defaultdict(list)
    for n, length in ncbi_lengths.items():
        if n not in used:
            by_length[length].append(n)
    for name, length in ensembl_lengths.items():
        if name not in mapping and len(by_length.get(length, [])) == 1:
            mapping[name] = by_length[length].pop()
    return mapping, [n for n in ensembl_lengths if n not in mapping]


def overlap(intervals: List[Tuple[int, int]], rows) -> Iterator[Tuple[tuple, int]]:
    """For each (sorted, non-overlapping) row of one sequence, its bases inside the sorted intervals."""
    j = 0
    for r in sorted(rows, key=lambda r: r[1]):
        while j < len(intervals) and intervals[j][1] <= r[1]:
            j += 1
        k, n = j, 0
        while k < len(intervals) and intervals[k][0] < r[2]:
            n += min(r[2], intervals[k][1]) - max(r[1], intervals[k][0])
            k += 1
        yield r, n


def main(argv=None) -> None:
    args = parse_args(argv)
    if not c.is_reference(args.assembly):
        c.stop(f"{args.assembly} is not a reference genome: Ensembl's masking is fetched for references only")
    folder = c.reference_folder(args.assembly)
    soft, report = folder / fe.LOCAL, folder / "assembly_report.txt"
    bed = c.version_folder(args.version) / "masked" / args.assembly / f"{args.assembly}.bed"
    for f in (soft, report, bed):
        if not f.exists():
            c.stop(f"{f} not found")
    out_dir = args.scratch / "comparisons" / "ensembl"
    targets = [out_dir / f"{args.version}_{args.assembly}.tsv", out_dir / f"{args.version}_{args.assembly}_classes.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    print(f"{args.version} on {args.assembly} against Ensembl's soft-masking ({soft.name})")
    if not args.go:
        print("dry run: nothing was read or written. Add --go to compare.")
        return

    print("reading Ensembl's soft-masked genome ...", flush=True)
    ens: Dict[str, List[Tuple[int, int]]] = {}
    ens_len: Dict[str, int] = {}
    with c.open_text(soft) as handle:
        for name, length, runs in softmasked_intervals(handle):
            ens[name], ens_len[name] = runs, length
    ncbi_len = c.fasta_lengths(c.genome_path(args.assembly, c.panel_accession(args.assembly)))
    with c.open_text(report) as handle:
        aliases = report_aliases(handle)
    mapping, unmapped = map_names(ens_len, aliases, ncbi_len)
    E = {mapping[n]: runs for n, runs in ens.items() if n in mapping}
    genome_bp = sum(ncbi_len.values())
    mapped_bp = sum(ncbi_len[t] for t in E)
    print(f"Ensembl sequences {len(ens):,}: mapped {len(mapping):,} ({100 * mapped_bp / genome_bp:.2f}% of the assembly)"
          + (f"; not mapped {len(unmapped)}: {', '.join(unmapped[:5])}" if unmapped else ""))
    ours = la.one_label_per_base(la.read_labelled_bed(bed))
    by_seq = defaultdict(list)
    for r in ours:
        by_seq[r[0]].append(r)
    e_bp = sum(b - a for runs in E.values() for a, b in runs)
    v_bp = both = 0
    cls_bp: Counter = Counter()
    cls_both: Counter = Counter()
    for seq, rows in by_seq.items():
        if seq not in E:
            continue
        for r, n in overlap(E[seq], rows):
            g = la.group(r[4])
            v_bp += r[2] - r[1]
            both += n
            cls_bp[g] += r[2] - r[1]
            cls_both[g] += n
    pct = lambda n, d: 100.0 * n / d if d else 0.0  # noqa: E731
    metrics = [("assembly bp (mapped)", mapped_bp), ("Ensembl masked bp", e_bp), ("version masked bp", v_bp),
               ("both bp", both), ("version only bp", v_bp - both), ("Ensembl only bp", e_bp - both),
               ("Ensembl masked %", round(pct(e_bp, mapped_bp), 2)), ("version masked %", round(pct(v_bp, mapped_bp), 2)),
               ("Ensembl-masked bases masked by version %", round(pct(both, e_bp), 2)),
               ("version-masked bases masked by Ensembl %", round(pct(both, v_bp), 2)),
               ("Ensembl sequences not mapped", len(unmapped))]
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["metric", "value"], metrics)
    c.write_tsv(targets[1], ["class", "version_bp", "also_masked_by_Ensembl_bp", "share_pct"],
                [[g, cls_bp[g], cls_both[g], f"{pct(cls_both[g], cls_bp[g]):.2f}"] for g, _ in cls_bp.most_common()])
    print(f"\nof the mapped assembly: Ensembl masks {pct(e_bp, mapped_bp):.2f}%, {args.version} {pct(v_bp, mapped_bp):.2f}%; "
          f"both {pct(both, mapped_bp):.2f}%")
    print(f"of Ensembl's masked bases, {args.version} masks {pct(both, e_bp):.2f}%; of {args.version}'s, Ensembl masks "
          f"{pct(both, v_bp):.2f}%")
    print("share of each of our classes that Ensembl also masks: "
          + ", ".join(f"{g} {pct(cls_both[g], cls_bp[g]):.1f}%" for g, _ in cls_bp.most_common()))
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
