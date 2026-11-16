#!/usr/bin/env python3
"""Remove gene-derived families from a library (gene filter, option C; decided 2026-11-13).

Evidence A, from gene_overlap.py's tables on the annotated genomes (--genomes): on every one of them,
  the family is the main masker of at least --a-host-genes genes of known multi-copy gene families
  (HOST_FAMILIES; half or more of their CDS masked), or masks CDS in at least --a-genes genes.
Evidence B, from blastx: the family's consensus matches a host protein (evalue <= --b-evalue) with a
  bitscore at least --b-ratio times its best match to a transposon protein (RepeatPeps; 0 when none).
  Host proteins: the longest protein of each gene of --proteins-from, translated from its genome and
  annotation.gff3.gz, leaving out genes whose product names a transposon protein.
Decision (2026-11-13): flagged by A and B: removed. Flagged by one: removed unless repbase_match.py
  found Repbase transposon support. Every flagged family is listed with its evidence.

Runs:
  no --go            the plan only; nothing is run.
  --go               runs the searches in <scratch>/<new version>/filter/ (reused if already there), writes
                     filter_table.tsv there and prints the counts. No version is written.
  --go --make-version  also writes the new version: <new>-families.prefixed.fa (the library without the
                     removed families), filter_table.tsv and filter_record.json in its library/.
Needs: gene_overlap.py tables for --version on each genome, and repbase_match.py's table for --version.

Example:
  python scripts/gene_filter.py --version RI_v5 --new-version RI_v6                    # plan
  python scripts/gene_filter.py --version RI_v5 --new-version RI_v6 --go               # report
  python scripts/gene_filter.py --version RI_v5 --new-version RI_v6 --go --make-version
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import gene_overlap as go  # noqa: E402

HOST_FAMILIES = re.compile(
    r"NACHT|NLR family|NLRC|NLRP|LRR and PYD|tripartite motif|TRIM\d|immunoglobulin|novel immune-type receptor|"
    r"T cell receptor|MHC|major histocompatibility|olfactory|odorant|taste receptor|vomeronasal|protocadherin|"
    r"zinc finger|histone|tubulin", re.IGNORECASE)
DEFAULT_TE_PROTEINS = Path.home() / ".RepeatModelerCache" / "PROT-Dfam_4.0" / "repeatpeps.lib"
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {a + b + c_: AMINO[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
          for k, c_ in enumerate(BASES)}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the library to filter, e.g. RI_v5")
    parser.add_argument("--new-version", required=True, help="the filtered version to write, e.g. RI_v6")
    parser.add_argument("--genomes", nargs="+", default=["DT_005m", "MZ_NMK2024a", "AC_fAstCal68"],
                        help="annotated genomes with gene_overlap.py tables (default: the three)")
    parser.add_argument("--proteins-from", default="MZ_NMK2024a", help="annotated genome for host proteins")
    parser.add_argument("--te-proteins", type=Path, default=DEFAULT_TE_PROTEINS,
                        help="transposon proteins, FASTA (default RepeatModeler's cached Dfam 4.0 RepeatPeps)")
    parser.add_argument("--a-host-genes", type=int, default=5)
    parser.add_argument("--a-genes", type=int, default=20)
    parser.add_argument("--b-evalue", type=float, default=1e-10)
    parser.add_argument("--b-ratio", type=float, default=2.0)
    parser.add_argument("--make-version", action="store_true", help="with --go: write the new version")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(nt: str) -> str:
    nt = nt.upper()
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def read_transcripts(lines: Iterable[str]) -> Dict[str, dict]:
    """transcript -> {sequence, strand, gene, product, segments [(start0, end, phase)]} from GFF3 CDS lines."""
    tx: Dict[str, dict] = {}
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = go.attributes(f[8])
        key = a.get("Parent") or a.get("ID")
        entry = tx.setdefault(key, {"sequence": f[0], "strand": f[6], "gene": a.get("gene") or key,
                                    "product": a.get("product", ""), "segments": []})
        entry["segments"].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7].isdigit() else 0))
    return tx


def protein_of(t: dict, genome: Dict[str, str]) -> str:
    segments = sorted(t["segments"], key=lambda s: s[0], reverse=(t["strand"] == "-"))
    chrom = genome[t["sequence"]]
    parts = [chrom[s:e] for s, e, _ in segments]
    nt = "".join(revcomp(p) for p in parts) if t["strand"] == "-" else "".join(parts)
    protein = translate(nt[segments[0][2]:])
    return protein[:-1] if protein.endswith("*") else protein


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with c.open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def host_proteins(transcripts: Dict[str, dict], genome: Dict[str, str]) -> Tuple[Dict[str, Tuple[str, str]], Counter]:
    """gene -> (product, longest protein), leaving out transposon products and proteins with internal stops."""
    best: Dict[str, Tuple[str, str]] = {}
    skipped = Counter()
    for t in transcripts.values():
        if go.TE_PRODUCT.search(t["product"]):
            skipped["transposon product"] += 1
            continue
        if t["sequence"] not in genome:
            skipped["sequence not in genome"] += 1
            continue
        protein = protein_of(t, genome)
        if "*" in protein or len(protein) < 30:
            skipped["internal stop or under 30 aa"] += 1
            continue
        if t["gene"] not in best or len(protein) > len(best[t["gene"]][1]):
            best[t["gene"]] = (t["product"], protein)
    return best, skipped


def evidence_a(genes_tables: Dict[str, List[dict]], families_tables: Dict[str, List[dict]],
               host_genes: int, many_genes: int) -> Dict[str, dict]:
    """family -> {"flag", per-genome host-family genes it mainly masks, per-genome genes it masks}."""
    host_main: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for genome, rows in genes_tables.items():
        for r in rows:
            if r["te_product"] == "no" and float(r["percent_masked"]) >= 50 and r["top_family"] \
                    and HOST_FAMILIES.search(r["product"]):
                host_main[r["top_family"]][genome] += 1
    hit: Dict[str, Dict[str, int]] = defaultdict(dict)
    for genome, rows in families_tables.items():
        for r in rows:
            hit[r["family"]][genome] = int(r["other_genes_hit"])
    genomes = list(genes_tables)
    out = {}
    for fam in set(host_main) | set(hit):
        hm = [host_main[fam].get(g, 0) for g in genomes]
        gh = [hit[fam].get(g, 0) for g in genomes]
        out[fam] = {"flag": all(x >= host_genes for x in hm) or all(x >= many_genes for x in gh),
                    "host_main": hm, "genes_hit": gh}
    return out


def best_hits(lines: Iterable[str]) -> Dict[str, Tuple[str, float, float]]:
    """query -> (subject, evalue, bitscore) of its highest-bitscore hit, from blast -outfmt '6 qseqid sseqid evalue bitscore'."""
    best: Dict[str, Tuple[str, float, float]] = {}
    for line in lines:
        f = line.split()
        if len(f) < 4:
            continue
        q, s, e, b = f[0], f[1], float(f[2]), float(f[3])
        if q not in best or b > best[q][2]:
            best[q] = (s, e, b)
    return best


def evidence_b(host: Dict[str, Tuple[str, float, float]], te: Dict[str, Tuple[str, float, float]],
               evalue: float, ratio: float) -> Set[str]:
    return {q for q, (_s, e, b) in host.items() if e <= evalue and b >= ratio * (te[q][2] if q in te else 0.0)}


def decide(in_a: bool, in_b: bool, repbase_support: bool) -> Tuple[str, str]:
    if in_a and in_b:
        return "remove", "A and B"
    if in_a or in_b:
        line = "A only" if in_a else "B only"
        return ("keep", f"{line}; Repbase transposon support") if repbase_support else ("remove", f"{line}; no Repbase support")
    return "keep", ""


def read_tsv(path: Path) -> List[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(argv=None) -> None:
    args = parse_args(argv)
    if args.new_version == args.version:
        c.stop("the filtered library must be a new version")
    _, _, library = c.library_files(args.version)
    genes_dir = args.scratch / "comparisons" / "genes"
    tables = {g: (genes_dir / f"{args.version}_{g}_genes.tsv", genes_dir / f"{args.version}_{g}_families.tsv")
              for g in args.genomes}
    repbase_table = args.scratch / "comparisons" / "repbase" / f"{args.version}_best_match.tsv"
    genome = c.genome_path(args.proteins_from, c.panel_accession(args.proteins_from))
    annotation = genome.parent / "annotation.gff3.gz"
    work = args.scratch / args.new_version / "filter"
    target = c.version_folder(args.new_version)
    print(f"library: {library}\nevidence A: gene_overlap.py tables for {', '.join(args.genomes)}\n"
          f"evidence B: host proteins of {args.proteins_from}; transposon proteins {args.te_proteins}\n"
          f"Repbase support: {repbase_table}\nwork (reused): {work}\n"
          + (f"new version: {target}" if args.make_version else "report only; add --make-version to write the version"))
    missing = [str(p) for pair in tables.values() for p in pair if not p.exists()]
    missing += [str(p) for p in (repbase_table, args.te_proteins, annotation) if not p.exists()]
    if missing:
        c.stop("missing: " + ", ".join(missing))
    if args.make_version and target.exists():
        c.stop(f"{target} already exists; a version is never changed")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    work.mkdir(parents=True, exist_ok=True)
    families = {f.name: f for f in c.read_families(library)}
    names = {f"q{i:06d}": name for i, name in enumerate(families, start=1)}
    short_of = {v: k for k, v in names.items()}
    query, host_fa = work / "families.fa", work / "host_proteins.fa"
    if not query.exists():
        seqs = read_fasta(library)
        with open(query, "x") as handle:
            for short, name in names.items():
                handle.write(f">{short}\n{seqs[name]}\n")
        c.write_tsv(work / "query_names.tsv", ["query", "family"], names.items())
    if not host_fa.exists():
        print(f"translating the CDS of {args.proteins_from} ...", flush=True)
        with gzip.open(annotation, "rt") as handle:
            transcripts = read_transcripts(handle)
        proteins, skipped = host_proteins(transcripts, read_fasta(genome))
        with open(host_fa, "x") as handle:
            for gene, (_product, protein) in proteins.items():
                handle.write(f">{gene}\n{protein}\n")
        c.write_tsv(work / "host_proteins.tsv", ["gene", "product", "length"],
                    [[g, p, len(s)] for g, (p, s) in proteins.items()])
        print(f"  {len(proteins):,} host proteins; left out: " + ", ".join(f"{k} {v:,}" for k, v in skipped.items()))
    for fasta, db in ((host_fa, work / "host_db"), (args.te_proteins, work / "te_db")):
        if not Path(str(db) + ".pdb").exists() and not Path(str(db) + ".phr").exists():
            c.run(c.container(args.sif, "makeblastdb", "-in", str(fasta), "-dbtype", "prot", "-out", str(db)),
                  True, work / "makeblastdb.log")
    results = {}
    for tag, db in (("host", work / "host_db"), ("te", work / "te_db")):
        out = work / f"blastx_{tag}.tsv"
        if not out.exists():
            print(f"blastx against the {'host' if tag == 'host' else 'transposon'} proteins ...", flush=True)
            c.run(c.container(args.sif, "blastx", "-query", str(query), "-db", str(db), "-evalue", "1e-5",
                              "-outfmt", "6 qseqid sseqid evalue bitscore", "-max_target_seqs", "5", "-max_hsps", "1",
                              "-num_threads", str(args.threads), "-out", str(out)), True, work / f"blastx_{tag}.log")
        with open(out) as handle:
            results[tag] = best_hits(handle)

    genes_tables = {g: read_tsv(t[0]) for g, t in tables.items()}
    families_tables = {g: read_tsv(t[1]) for g, t in tables.items()}
    a = evidence_a(genes_tables, families_tables, args.a_host_genes, args.a_genes)
    b = evidence_b(results["host"], results["te"], args.b_evalue, args.b_ratio)
    repbase = {r["family"]: r for r in read_tsv(repbase_table)}
    rows, removed = [], []
    for name, fam in families.items():
        ev_a = a.get(fam.family, {})   # gene tables name families without their #class suffix
        in_a = ev_a.get("flag", False)
        in_b = short_of[name] in b
        if not (in_a or in_b):
            continue
        rb = repbase.get(name, {})
        decision, reason = decide(in_a, in_b, rb.get("transposon_support") == "yes")
        host_hit, te_hit = results["host"].get(short_of[name]), results["te"].get(short_of[name])
        rows.append([name, fam.label, fam.top_class, fam.length, "yes" if in_a else "no",
                     "/".join(map(str, ev_a.get("host_main", []))), "/".join(map(str, ev_a.get("genes_hit", []))),
                     "yes" if in_b else "no", host_hit[0] if host_hit else "", f"{host_hit[2]:.0f}" if host_hit else "",
                     te_hit[0] if te_hit else "", f"{te_hit[2]:.0f}" if te_hit else "",
                     rb.get("transposon_support", ""), rb.get("repbase_entry", ""), rb.get("repbase_label", ""),
                     rb.get("coverage_of_family", ""), decision, reason])
        if decision == "remove":
            removed.append(name)
    header = ["family", "label", "class", "length", "evidence_A", "host_family_genes_mainly_masked", "genes_masked",
              "evidence_B", "best_host_protein", "host_bitscore", "best_te_protein", "te_bitscore",
              "repbase_transposon_support", "repbase_entry", "repbase_label", "repbase_coverage", "decision", "reason"]
    report = work / "filter_table.tsv"
    if not report.exists():
        c.write_tsv(report, header, rows)
    n_a, n_b = sum(1 for r in rows if r[4] == "yes"), sum(1 for r in rows if r[7] == "yes")
    both = sum(1 for r in rows if r[4] == r[7] == "yes")
    reasons = Counter(r[17] for r in rows)
    print(f"\n{args.version}: {len(families):,} families; evidence A {n_a}, evidence B {n_b}, both {both}")
    for reason, n in reasons.most_common():
        print(f"  {reason}: {n}")
    print(f"removed: {len(removed)} families ({', '.join(f'{k} {v}' for k, v in Counter(families[f].top_class for f in removed).most_common())})")
    print(f"table: {report}")
    if not args.make_version:
        print("report only: no version written. Review the table, then add --make-version.")
        return

    library_dir = target / "library"
    library_dir.mkdir(parents=True)
    seqs = read_fasta(library)
    removed_set = set(removed)
    kept = [n for n in families if n not in removed_set]
    out_lib = library_dir / f"{args.new_version}-families.prefixed.fa"
    with open(out_lib, "x") as handle:
        for name in kept:
            s = seqs[name]
            handle.write(f">{name}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
    c.write_tsv(library_dir / "filter_table.tsv", header, rows)
    record = {"version": args.new_version, "kind": "filtered library", "source_version": args.version,
              "source_library": {"path": str(library), "sha256": c.sha256(library)},
              "te_proteins": {"path": str(args.te_proteins), "sha256": c.sha256(args.te_proteins)},
              "host_proteins_from": args.proteins_from, "annotation": str(annotation),
              "gene_tables": {g: [str(p) for p in t] for g, t in tables.items()}, "repbase_table": str(repbase_table),
              "thresholds": {"a_host_genes": args.a_host_genes, "a_genes": args.a_genes, "b_evalue": args.b_evalue,
                             "b_ratio": args.b_ratio},
              "host_families_pattern": HOST_FAMILIES.pattern, "rule": "A and B: remove; one line: remove unless Repbase transposon support",
              "counts": {"families_in": len(families), "removed": len(removed), "kept": len(kept),
                         "evidence_A": n_a, "evidence_B": n_b, "both": both, "by_reason": dict(reasons)},
              "written": f"{datetime.now():%Y-%m-%d %H:%M}",
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(library_dir / "filter_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"written: {library_dir}/ ({out_lib.name}, filter_table.tsv, filter_record.json); {len(kept):,} families kept")


if __name__ == "__main__":
    main()
