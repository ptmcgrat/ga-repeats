#!/usr/bin/env python3
"""The known insertions (question 2, briefing C2.4): what a version's masking finds at rx1 and mitfa.

Carleton et al. 2020 (Molecular Ecology) and Quah et al. 2025 describe transposon insertions at two loci studied
here (tbx2a and fhl2b wait until their genes are identified in the annotations):

  rx1    about 2.5 kb 5' of the coding sequence. Most Lake Malawi species carry a 413 bp deletion there (short and
         medium palette, M. zebra among them); long-palette species carry an insertion of 268-421 bp that mostly
         matches Rex1-5_AFC (LINE/Rex-Babar); A. calliptera and some Rhamphochromis have neither.
  mitfa  intron 1. A 1,408 bp insertion: a 522 bp Rex1 match (LINE/Rex-Babar) plus 50-70 bp fragments of L1, Mariner
         and Copia. Common (20 of 74 individuals homozygous, 12 heterozygous) but absent from A. calliptera and nearly
         absent from mbuna.

So M. zebra and A. calliptera are negative controls at both loci, and which other panel species carry the
insertions is not known: the check reports what is there in every annotated genome, so carriers stand out by a
longer region and a Rex-Babar block.

For each assembly (it needs annotation.gff3.gz and the version's masking), for every gene named rx1 or mitfa:
  rx1    the --upstream bases (default 5,000) 5' of the coding start of each transcript (distinct windows only);
  mitfa  intron 1 of each transcript (between its first and second exon; distinct introns only);
and in each region every repeat annotation (one label per base: the highest-scoring), its class and length,
the masked bases, and the bases labelled LINE/Rex-Babar.

Usage (from the repository root):
  python scripts/known_insertions.py --version RI_v5 --assemblies DT_005m MZ_NMK2024a AC_fAstCal68          # dry run
  python scripts/known_insertions.py --version RI_v5 --assemblies DT_005m MZ_NMK2024a AC_fAstCal68 --go

Writes, never overwriting, <scratch>/comparisons/known_insertions/<version>_known_insertions.tsv
(one row per region) and <version>_known_insertions_annotations.tsv (one row per annotation in a region).
"""
import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import gene_overlap as go  # noqa: E402
import label_agreement as la  # noqa: E402

GENES = {"rx1": "upstream", "mitfa": "intron 1"}
REX = "Rex-Babar"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a version that has masked the assemblies, e.g. RI_v5")
    parser.add_argument("--assemblies", nargs="+", required=True, help="annotated assemblies, e.g. DT_005m MZ_NMK2024a")
    parser.add_argument("--upstream", type=int, default=5000, help="bases 5' of rx1's coding start to report (default 5000)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def read_gene_models(lines: Iterable[str], wanted: Iterable[str]) -> List[dict]:
    """Genes whose Name (or gene) attribute is in wanted: seq, strand, and per transcript its exons and CDS (0-based, half-open)."""
    wanted = {w.lower() for w in wanted}
    genes: Dict[str, dict] = {}
    transcript_gene: Dict[str, str] = {}
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = go.attributes(f[8])
        if f[2] == "gene":
            name = (a.get("Name") or a.get("gene") or "").lower()
            if name in wanted:
                genes[a.get("ID", "")] = {"name": name, "id": a.get("ID", ""), "seq": f[0], "strand": f[6], "transcripts": {}}
        elif a.get("Parent") in genes and f[2] in ("mRNA", "transcript"):
            transcript_gene[a.get("ID", "")] = a["Parent"]
            genes[a["Parent"]]["transcripts"][a.get("ID", "")] = {"exons": [], "cds": []}
        elif f[2] in ("exon", "CDS") and a.get("Parent") in transcript_gene:
            gene = genes[transcript_gene[a["Parent"]]]
            gene["transcripts"][a["Parent"]]["exons" if f[2] == "exon" else "cds"].append((int(f[3]) - 1, int(f[4])))
    return list(genes.values())


def intron1(exons: List[Tuple[int, int]], strand: str) -> Optional[Tuple[int, int]]:
    """Intron 1 in transcript order: between the first and second exon."""
    if len(exons) < 2:
        return None
    ordered = sorted(exons) if strand == "+" else sorted(exons, reverse=True)
    first, second = ordered[0], ordered[1]
    return (first[1], second[0]) if strand == "+" else (second[1], first[0])


def upstream(cds: List[Tuple[int, int]], strand: str, length: int) -> Optional[Tuple[int, int]]:
    """The length bases 5' of the coding start."""
    if not cds:
        return None
    if strand == "+":
        start = min(s for s, _ in cds)
        return max(0, start - length), start
    end = max(e for _, e in cds)
    return end, end + length


def regions(gene: dict, upstream_length: int) -> List[Tuple[int, int]]:
    kind = GENES[gene["name"]]
    found = set()
    for t in gene["transcripts"].values():
        r = upstream(t["cds"], gene["strand"], upstream_length) if kind == "upstream" else intron1(t["exons"], gene["strand"])
        if r is not None and r[1] > r[0]:
            found.add(r)
    return sorted(found)


def annotations_in(rows: List[la.Row], seq: str, start: int, end: int) -> List[Tuple[int, int, str, str]]:
    """(start, end, family, label) of the (one label per base) stretches inside [start, end) of seq, clipped."""
    out = []
    for s, a, b, family, label in rows:
        if s == seq and b > start and a < end:
            out.append((max(a, start), min(b, end), family, label))
    return sorted(out)


def main(argv=None) -> None:
    args = parse_args(argv)
    plan = []
    for assembly in args.assemblies:
        genome = c.genome_path(assembly, c.panel_accession(assembly))
        gff = genome.parent / "annotation.gff3.gz"
        bed = c.version_folder(args.version) / "masked" / assembly / f"{assembly}.bed"
        for f, what in ((gff, "annotation"), (bed, f"{args.version}'s masking")):
            if not f.exists():
                c.stop(f"{assembly}: no {what} ({f})")
        plan.append((assembly, gff, bed))
    out_dir = args.scratch / "comparisons" / "known_insertions"
    targets = [out_dir / f"{args.version}_known_insertions.tsv", out_dir / f"{args.version}_known_insertions_annotations.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    print(f"{args.version}: rx1 ({args.upstream:,} bp upstream of the coding start) and mitfa (intron 1) in "
          + ", ".join(a for a, _, _ in plan))
    if not args.go:
        print("dry run: nothing was read or written. Add --go to run.")
        return

    region_rows, annotation_rows = [], []
    for assembly, gff, bed in plan:
        with gzip.open(gff, "rt") as handle:
            genes = read_gene_models(handle, GENES)
        rows = la.one_label_per_base(la.read_labelled_bed(bed))
        names = Counter(g["name"] for g in genes)
        for missing in sorted(set(GENES) - set(names)):
            print(f"{assembly}: no gene named {missing}")
        for gene in sorted(genes, key=lambda g: g["name"]):
            for start, end in regions(gene, args.upstream):
                hits = annotations_in(rows, gene["seq"], start, end)
                masked = sum(b - a for a, b, _, _ in hits)
                rex = sum(b - a for a, b, _, label in hits if REX in label)
                by_class = Counter()
                for a, b, _, label in hits:
                    by_class[la.group(label)] += b - a
                region_rows.append([assembly, gene["name"], gene["id"], GENES[gene["name"]], gene["seq"], start, end, gene["strand"],
                                    end - start, masked, rex, "; ".join(f"{k} {v}" for k, v in by_class.most_common())])
                for a, b, family, label in hits:
                    rel = (a - start, b - start) if gene["strand"] == "+" else (end - b, end - a)
                    annotation_rows.append([assembly, gene["name"], GENES[gene["name"]], gene["seq"], start, end, family, label,
                                            la.group(label), b - a, rel[0], rel[1]])
                print(f"{assembly} {gene['name']} {GENES[gene['name']]} {gene['seq']}:{start + 1:,}-{end:,} ({gene['strand']}): "
                      f"{end - start:,} bp, masked {masked:,} bp ({100 * masked / (end - start):.0f}%), Rex-Babar {rex:,} bp; "
                      + (", ".join(f"{label} {b - a}" for a, b, _, label in hits) or "no repeats"))
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["assembly", "gene", "gene_id", "region", "seq", "start", "end", "strand", "length", "masked_bp",
                             "rex_babar_bp", "by_class"], region_rows)
    c.write_tsv(targets[1], ["assembly", "gene", "region", "seq", "region_start", "region_end", "family", "label", "class", "bp",
                             "start_in_region", "end_in_region"], annotation_rows)
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
