#!/usr/bin/env python3
"""TE composition of an assembly's masking, beside published compositions of the same assembly (question 2).

Like for like (spec section 6, 2026-11-23: the reference is published data, on the same assembly, with the same
definitions). The version's masking of --assembly (one label per base: each masked base counted once, under its
highest-scoring annotation) is grouped as RepeatMasker's summary table groups it:
  TEs      DNA, LINE (Penelope/PLE counted with the LINEs, as RepeatMasker's summary does), LTR, SINE, RC (helitrons),
           Retroposon, Unknown
  not TEs  Satellite, Simple_repeat, Low_complexity, small RNA (rRNA, tRNA, snRNA, scRNA, srpRNA), other
Each class is given as a share of the whole assembly and of its non-N bases (published shares often do not say which).
Beside it:
  NCBI     NCBI's own RepeatMasker run on the assembly (<scratch>/references/<id>/ncbi_rm.out.gz), where fetched
  published  references/published_te_composition.tsv: the values, their source and wording (for fAstCal1.2: Quah et al.
           2025, whose "TEs" are the seven TE classes above, Unknown included)
Families tagged as gene-derived (RI_v6 and later, '_hostgene-' in the name) stay in their class and are also shown on a
line of their own.

Usage (from the repository root):
  python scripts/te_composition.py --version RI_v5 --assembly AC_fAstCal1_2          # dry run
  python scripts/te_composition.py --version RI_v5 --assembly AC_fAstCal1_2 --go

Writes, never overwriting, <scratch>/comparisons/te_composition/<version>_<assembly>.tsv.
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

TE_CLASSES = ["DNA", "LINE", "LTR", "SINE", "RC", "Retroposon", "Unknown"]
OTHER_CLASSES = ["Satellite", "Simple_repeat", "Low_complexity", "small RNA", "other"]
SMALL_RNA = {"rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "RNA"}
PUBLISHED = Path(__file__).resolve().parent.parent / "references" / "published_te_composition.tsv"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True)
    parser.add_argument("--assembly", default="AC_fAstCal1_2")
    parser.add_argument("--published", type=Path, default=PUBLISHED)
    c.add_common_args(parser)
    return parser.parse_args(argv)


def category(label: str) -> str:
    """A RepeatMasker label in the classes of RepeatMasker's summary table."""
    top = label.split("/", 1)[0].rstrip("?")
    if top == "PLE":
        return "LINE"
    if top in TE_CLASSES or top in ("Satellite", "Simple_repeat", "Low_complexity"):
        return top
    if top in SMALL_RNA:
        return "small RNA"
    return "other"


def composition(rows) -> Counter:
    """One-label-per-base rows -> bases per class, plus 'TEs' and 'gene-derived (tagged)'."""
    out: Counter = Counter()
    for _s, a, b, family, label in rows:
        cat = category(label)
        out[cat] += b - a
        if cat in TE_CLASSES:
            out["TEs"] += b - a
        if "_hostgene-" in family:
            out["gene-derived (tagged)"] += b - a
    return out


def n_bases(path: Path) -> int:
    n = 0
    with c.open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                n += line.count("N") + line.count("n")
    return n


def ncbi_rows(path: Path) -> List[la.Row]:
    scored = [(h.query, h.begin - 1, h.end, h.repeat, h.label, h.score) for h in c.parse_out(path)]
    return la.one_label_per_base(scored)


def published_values(path: Path, assembly: str) -> Dict[str, str]:
    if not path.exists():
        return {}
    with open(path, newline="") as handle:
        return {r["class"]: r["percent"] for r in csv.DictReader(handle, delimiter="\t") if r["assembly"] == assembly}


def main(argv=None) -> None:
    args = parse_args(argv)
    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    bed = c.version_folder(args.version) / "masked" / args.assembly / f"{args.assembly}.bed"
    for f in (genome, bed):
        if not f.exists():
            c.stop(f"{f} not found" + (f": run mask_genomes.py --version {args.version} --assemblies {args.assembly}"
                                       if f == bed else ""))
    target = args.scratch / "comparisons" / "te_composition" / f"{args.version}_{args.assembly}.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    ncbi = c.reference_folder(args.assembly) / "ncbi_rm.out.gz" if c.is_reference(args.assembly) else None
    published = published_values(args.published, args.assembly)
    print(f"{args.version} on {args.assembly}" + (f"; beside NCBI's RepeatMasker run" if ncbi and ncbi.exists() else "")
          + (f"; beside the published composition ({len(published)} values)" if published else "; no published values"))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to compute.")
        return

    total = sum(c.fasta_lengths(genome).values())
    non_n = total - n_bases(genome)
    ours = composition(la.one_label_per_base(la.read_labelled_bed(bed)))
    theirs: Optional[Counter] = composition(ncbi_rows(ncbi)) if ncbi and ncbi.exists() else None
    pct = lambda n, d: 100.0 * n / d if d else 0.0  # noqa: E731
    order = ["TEs"] + TE_CLASSES + ["SINE+RC+Retroposon"] + OTHER_CLASSES + ["gene-derived (tagged)"]
    for counts in [ours] + ([theirs] if theirs else []):
        counts["SINE+RC+Retroposon"] = counts["SINE"] + counts["RC"] + counts["Retroposon"]
    rows = []
    print(f"\n{args.assembly}: {total:,} bp, of which {total - non_n:,} N ({pct(total - non_n, total):.2f}%)")
    head = f"{'class':24s} {'ours, % of all':>15s} {'% of non-N':>11s}" + (f" {'NCBI, % all':>12s}" if theirs else "") \
        + f" {'published %':>12s} {'ours - pub.':>12s}"
    print(head)
    for k in order:
        pub = published.get(k, "")
        diff = ""
        try:
            diff = f"{pct(ours[k], total) - float(pub):+.2f}"
        except ValueError:
            pass
        row = [k, f"{pct(ours[k], total):.2f}", f"{pct(ours[k], non_n):.2f}",
               f"{pct(theirs[k], total):.2f}" if theirs else "", pub, diff]
        rows.append(row)
        if k == "Satellite":
            print("  -- not TEs --")
        print(f"{k:24s} {row[1]:>15s} {row[2]:>11s}" + (f" {row[3]:>12s}" if theirs else "") + f" {pub or '-':>12s} {diff:>12s}")
    target.parent.mkdir(parents=True, exist_ok=True)
    c.write_tsv(target, ["class", "ours_pct_of_assembly", "ours_pct_of_non_N", "ncbi_pct_of_assembly", "published_pct",
                         "ours_minus_published"], rows)
    if published:
        print(f"\npublished: {args.published} (values, source and wording)")
    print(f"table: {target}")


if __name__ == "__main__":
    main()
