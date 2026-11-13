#!/usr/bin/env python3
"""How much protein-coding sequence does a library mask? (false positives; question 2)

On an annotated genome masked by a version: the coding sequence (the CDS features of the genome's
annotation.gff3.gz, merged across each gene's transcripts) that the library masks, by class, by
family and per gene. A transposable-element library should mask little coding sequence. Simple
repeats in low-complexity protein stretches are expected; a family that masks coding sequence in
many genes is likely gene-derived (a false family) or carries a captured gene fragment.

Genes whose annotated product names a transposable-element protein (see TE_PRODUCT below:
transposase, reverse transcriptase, integrase, ...) are counted separately: masking those is
expected, not a false positive. Each masked base keeps one label (its highest-scoring annotation).

Inputs: annotation.gff3.gz next to the genome (<data root>/Genomes/<assembly>/<accession>/), and
<version folder>/masked/<assembly>/<assembly>.bed.
Results in <scratch>/comparisons/genes/: <version>_<assembly>_classes.tsv, _families.tsv and
_genes.tsv. Never overwritten.

Example:
  python scripts/gene_overlap.py --version RI_v5 --assembly DT_005m          # dry run
  python scripts/gene_overlap.py --version RI_v5 --assembly DT_005m --go
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from urllib.parse import unquote

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

TE_PRODUCT = re.compile(
    r"transpos|retrotranspos|reverse transcriptase|integrase|\bgag\b|\bpol\b|piggybac|tigger|\bpogo\b|mariner|"
    r"\bhat\b|harbinger|helitron|endogenous retrovir|line-1|\bl1re", re.IGNORECASE)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a version that has masked the assembly, e.g. RI_v5")
    parser.add_argument("--assembly", required=True, help="an annotated panel assembly, e.g. DT_005m")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def attributes(text: str) -> Dict[str, str]:
    out = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = unquote(value)
    return out


def read_cds(lines: Iterable[str]) -> Dict[str, dict]:
    """gene -> {"product", "te_like", "spans": {sequence: merged [(start, end)]}} from GFF3 CDS lines (0-based, half-open)."""
    raw: Dict[str, dict] = {}
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = attributes(f[8])
        gene = a.get("gene") or a.get("Parent") or a.get("ID")
        entry = raw.setdefault(gene, {"product": a.get("product", ""), "intervals": []})
        if not entry["product"] and a.get("product"):
            entry["product"] = a["product"]
        entry["intervals"].append((f[0], int(f[3]) - 1, int(f[4])))
    genes = {}
    for gene, entry in raw.items():
        genes[gene] = {"product": entry["product"], "te_like": bool(TE_PRODUCT.search(entry["product"])),
                       "spans": c.merge_by_sequence(entry["intervals"])}
    return genes


def index_rows(rows: List[la.Row]) -> Dict[str, Tuple[List[la.Row], List[int]]]:
    """Per sequence: the masked stretches sorted by start, and their ends (non-overlapping, so also sorted)."""
    by_seq: Dict[str, List[la.Row]] = defaultdict(list)
    for r in rows:
        by_seq[r[0]].append(r)
    index = {}
    for name, items in by_seq.items():
        items.sort(key=lambda r: r[1])
        index[name] = (items, [r[2] for r in items])
    return index


def overlaps(spans: Dict[str, list], index) -> Iterable[Tuple[str, int, int, str, str, int]]:
    """(sequence, start, end, family, label, overlap bp) for each masked stretch meeting the spans."""
    for name, targets in spans.items():
        if name not in index:
            continue
        items, ends = index[name]
        for start, end in targets:
            k = bisect.bisect_right(ends, start)
            while k < len(items) and items[k][1] < end:
                bp = min(end, items[k][2]) - max(start, items[k][1])
                if bp > 0:
                    yield (name, start, end, items[k][3], items[k][4], bp)
                k += 1


def summarise(genes: Dict[str, dict], rows: List[la.Row]):
    """(class totals per gene set, per-family, per-gene) for the masked rows (one label per base)."""
    index = index_rows(rows)
    sets = {"all genes": [g for g in genes], "genes without TE products": [g for g, v in genes.items() if not v["te_like"]],
            "genes with TE products": [g for g, v in genes.items() if v["te_like"]]}
    classes = {}
    for name, members in sets.items():
        union = c.merge_by_sequence((s, a, b) for g in members for s, spans in genes[g]["spans"].items() for a, b in spans)
        by_class: Dict[str, int] = defaultdict(int)
        for _s, _a, _b, _fam, label, bp in overlaps(union, index):
            by_class[la.group(label)] += bp
        classes[name] = (c.merged_length(union), dict(by_class))
    families: Dict[Tuple[str, str], dict] = {}
    per_gene = {}
    for gene, v in genes.items():
        fam_bp: Dict[Tuple[str, str], int] = defaultdict(int)
        for _s, _a, _b, fam, label, bp in overlaps(v["spans"], index):
            fam_bp[(fam, label)] += bp
        cds_bp = sum(e - s for spans in v["spans"].values() for s, e in spans)
        masked = sum(fam_bp.values())
        top = max(fam_bp.items(), key=lambda kv: kv[1])[0] if fam_bp else ("", "")
        per_gene[gene] = (v["product"], v["te_like"], cds_bp, masked, top)
        for key, bp in fam_bp.items():
            entry = families.setdefault(key, {"bp": 0, "genes": [], "te_genes": 0})
            if v["te_like"]:
                entry["te_genes"] += 1
            else:
                entry["bp"] += bp
                entry["genes"].append(gene)
    return classes, families, per_gene


def main(argv=None) -> None:
    args = parse_args(argv)
    genome = c.genome_path(args.assembly, c.panel_accession(args.assembly))
    annotation = genome.parent / "annotation.gff3.gz"
    bed = c.version_folder(args.version) / "masked" / args.assembly / f"{args.assembly}.bed"
    out_dir = args.scratch / "comparisons" / "genes"
    targets = [out_dir / f"{args.version}_{args.assembly}_{part}.tsv" for part in ("classes", "families", "genes")]
    print(f"annotation: {annotation}\nmasking: {bed}\nresults: {out_dir}/{args.version}_{args.assembly}_*.tsv")
    if not annotation.exists():
        c.stop(f"{args.assembly} has no annotation: {annotation} not found")
    if not bed.exists():
        c.stop(f"{args.version} has not masked {args.assembly}: {bed} not found")
    existing = [p for p in targets if p.exists()]
    if existing:
        c.stop("already exists, never overwritten: " + ", ".join(map(str, existing)))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    with gzip.open(annotation, "rt") as handle:
        genes = read_cds(handle)
    if not genes:
        c.stop(f"no CDS features in {annotation}")
    lengths = c.fasta_lengths(genome)
    cds_bp = {s: sum(e - a for a, e in spans) for v in genes.values() for s, spans in v["spans"].items()}
    known = sum(bp for s, bp in cds_bp.items() if s in lengths)
    if known < 0.95 * sum(cds_bp.values()):
        c.stop("the annotation's sequence names do not match the genome's: e.g. annotation "
               + ", ".join(list(cds_bp)[:3]) + "; genome " + ", ".join(list(lengths)[:3]))
    print(f"{len(genes):,} genes with CDS; reading the masking ...", flush=True)
    rows = la.one_label_per_base(la.read_labelled_bed(bed))
    classes, families, per_gene = summarise(genes, rows)

    out_dir.mkdir(parents=True, exist_ok=True)
    table = []
    for name, (total, by_class) in classes.items():
        masked = sum(by_class.values())
        table.append([args.assembly, args.version, name, "ALL", total, masked, f"{100 * masked / total:.3f}" if total else ""])
        for g in la.GROUPS:
            bp = by_class.get(g, 0)
            table.append([args.assembly, args.version, name, g, total, bp, f"{100 * bp / total:.3f}" if total else ""])
    c.write_tsv(targets[0], ["assembly", "version", "gene_set", "class", "cds_bp", "masked_bp", "percent_of_cds"], table)
    fam_rows = sorted(families.items(), key=lambda kv: (-len(kv[1]["genes"]), -kv[1]["bp"]))
    c.write_tsv(targets[1], ["family", "label", "class", "cds_bp_masked_in_other_genes", "other_genes_hit",
                             "te_product_genes_hit", "example_genes"],
                [[fam, label, la.group(label), v["bp"], len(v["genes"]), v["te_genes"], ",".join(sorted(v["genes"])[:8])]
                 for (fam, label), v in fam_rows])
    gene_rows = sorted(((g, *v) for g, v in per_gene.items() if v[3]), key=lambda r: -r[4] / r[3])
    c.write_tsv(targets[2], ["gene", "product", "te_product", "cds_bp", "masked_bp", "percent_masked", "top_family", "top_label"],
                [[g, prod, "yes" if te else "no", cds, masked, f"{100 * masked / cds:.1f}", top[0], top[1]]
                 for g, prod, te, cds, masked, top in gene_rows])

    print()
    for name, (total, by_class) in classes.items():
        masked = sum(by_class.values())
        n = sum(1 for v in genes.values() if (name == "all genes" or (v["te_like"] == (name == "genes with TE products"))))
        print(f"{name} ({n:,}): {total:,} bp of CDS, {100 * masked / total:.2f}% masked — "
              + ", ".join(f"{g} {100 * by_class.get(g, 0) / total:.2f}%" for g in la.GROUPS if by_class.get(g)))
    other = [r for r in gene_rows if not r[2]]
    half = sum(1 for r in other if r[4] >= 0.5 * r[3])
    print(f"genes without TE products with half or more of their CDS masked: {half:,} of {len(genes) - sum(v['te_like'] for v in genes.values()):,}")
    print("families masking CDS in the most genes without TE products:")
    for (fam, label), v in [kv for kv in fam_rows if kv[1]["genes"]][:10]:
        print(f"  {fam} ({label}): {len(v['genes'])} genes, {v['bp']:,} bp")
    print(f"written: {', '.join(p.name for p in targets)} in {out_dir}")


if __name__ == "__main__":
    main()
