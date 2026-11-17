#!/usr/bin/env python3
"""Tag gene-derived families in a library (gene filter, option C; decided 2026-11-13; tags instead of removal 2026-11-17).

A family judged gene-derived keeps its sequence and class but is renamed <family>_hostgene-<CATEGORY>#<class>,
so masking and class shares stay as they were and any analysis can include or exclude tagged families.
CATEGORY (CATEGORIES below): NLR, TRIM, IMMUNE, OR, TAAR, TASTE, PCDH, ZNF, HIST, TUB, OTHER (a named host gene
outside these), MULTI (no protein match, but masks coding sequence in many named genes).

Evidence A, from gene_overlap.py's tables on the annotated genomes (--genomes): on every one of them,
  the family is the main masker of at least --a-host-genes genes of known multi-copy gene families
  (HOST_FAMILIES; half or more of their CDS masked), or masks CDS in at least --a-genes named genes
  (uncharacterized gene models, mostly built from transposons, are not counted; repeats_020).
Evidence B, from blastx: the family's consensus matches a host protein (evalue <= --b-evalue) with a
  bitscore at least --b-ratio times its best match to a transposon protein (RepeatPeps; 0 when none).
  Host proteins: the longest protein of each gene of --proteins-from, translated from its genome and
  annotation.gff3.gz, keeping only named proteins of host genes: products called "uncharacterized" (mostly
  gene models built from transposons), products naming a transposon protein, and proteins that blastp
  matches to RepeatPeps (evalue <= --host-te-evalue; domesticated transposon genes such as ZBED or NYNRIN)
  are left out. How much of each family the host match covers is recorded (repeats_019; the first version
  kept uncharacterized and transposon-derived host proteins and flagged real transposons).
  Evidence B is strong when the host match covers at least --b-min-coverage of the family, or the host
  gene belongs to a known multi-copy gene family (HOST_FAMILIES); otherwise weak, and weak B does not count.
Decision (2026-11-13, refined 2026-11-14): flagged by A and strong B: tagged. Flagged by one: tagged unless
  repbase_match.py found Repbase transposon support. Weak B alone: not tagged, listed (a family that is mostly
  transposon with a small host-like stretch). Every flagged family is listed with its evidence. The report also
  counts what the alternative rule (strong B by coverage only, without the gene-family exception) would tag.
The report shows how much the tagged families mask: of the CDS of each annotated genome, and of --impact-genome
  (default DT_004f).

Runs:
  no --go            the plan only; nothing is run.
  --go               runs the searches in <scratch>/<new version>/filter/ (reused if already there), writes
                     filter_report.tsv there and prints the counts. No version is written.
  --go --make-version  also writes the new version: <new>-families.prefixed.fa (every family, the tagged ones
                     renamed), tag_table.tsv and filter_record.json in its library/.
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
import label_agreement as la  # noqa: E402

HOST_FAMILIES = re.compile(
    r"NACHT|NLR family|NLRC|NLRP|LRR and PYD|tripartite motif|TRIM\d|immunoglobulin|novel immune-type receptor|"
    r"T cell receptor|MHC|major histocompatibility|olfactory|odorant|taste receptor|vomeronasal|protocadherin|"
    r"zinc finger|\bhistone H|tubulin (alpha|beta)|trace amine", re.IGNORECASE)
CATEGORIES = [   # tag categories (decided 2026-11-17), checked in this order
    ("NLR", re.compile(r"NACHT|NLR family|NLRC|NLRP|LRR and PYD", re.I)),
    ("TRIM", re.compile(r"tripartite motif|TRIM\d", re.I)),
    ("IMMUNE", re.compile(r"immunoglobulin|novel immune-type receptor|T cell receptor|MHC|major histocompatibility", re.I)),
    ("OR", re.compile(r"olfactory|odorant", re.I)),
    ("TAAR", re.compile(r"trace amine", re.I)),
    ("TASTE", re.compile(r"taste receptor|vomeronasal", re.I)),
    ("PCDH", re.compile(r"protocadherin", re.I)),
    ("ZNF", re.compile(r"zinc finger", re.I)),
    ("HIST", re.compile(r"\bhistone H", re.I)),
    ("TUB", re.compile(r"tubulin (alpha|beta)", re.I)),
]


def category_of(product: str) -> str:
    for name, pattern in CATEGORIES:
        if pattern.search(product or ""):
            return name
    return "OTHER"


def tagged_name(family: str, category: str, label: str) -> str:
    return f"{family}_hostgene-{category}#{label}"


def tag_category(b_product: str, b_strength: str, a: dict) -> str:
    """The best host protein's category if it names a known family, else the family A counted, else what is left."""
    cat_b = category_of(b_product) if b_strength == "strong" else ""
    cat_a = a.get("category", "") if a.get("flag") else ""
    for candidate in (cat_b if cat_b != "OTHER" else "", cat_a if cat_a != "MULTI" else "", cat_b, cat_a):
        if candidate:
            return candidate
    return "OTHER"
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
    parser.add_argument("--b-min-coverage", type=float, default=0.30,
                        help="host match coverage of the family for strong evidence B (default 0.30)")
    parser.add_argument("--impact-genome", default="DT_004f",
                        help="genome masked by --version on which to measure how much tagged families mask (default DT_004f)")
    parser.add_argument("--host-te-evalue", type=float, default=1e-5,
                        help="a host protein matching RepeatPeps at this evalue or lower is left out (default 1e-5)")
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


UNCHARACTERIZED = re.compile(r"uncharacterized", re.IGNORECASE)
OLD_METHOD_FILES = ("host_proteins.fa", "host_proteins.tsv", "blastx_host.tsv", "filter_table.tsv")


def host_proteins(transcripts: Dict[str, dict], genome: Dict[str, str]) -> Tuple[Dict[str, Tuple[str, str]], Counter]:
    """gene -> (product, longest protein) for named host genes: leaves out uncharacterized and transposon
    products, and proteins with internal stops or under 30 aa."""
    best: Dict[str, Tuple[str, str]] = {}
    skipped = Counter()
    for t in transcripts.values():
        if UNCHARACTERIZED.search(t["product"]) or not t["product"]:
            skipped["uncharacterized or unnamed"] += 1
            continue
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
    host_cats: Dict[str, Counter] = defaultdict(Counter)
    for genome, rows in genes_tables.items():
        for r in rows:
            if r["te_product"] == "no" and float(r["percent_masked"]) >= 50 and r["top_family"] \
                    and HOST_FAMILIES.search(r["product"]):
                host_main[r["top_family"]][genome] += 1
                host_cats[r["top_family"]][category_of(r["product"])] += 1
    hit: Dict[str, Dict[str, int]] = defaultdict(dict)
    for genome, rows in families_tables.items():
        for r in rows:
            if "named_genes_hit" not in r:
                c.stop(f"the gene table for {genome} has no named_genes_hit column: run gene_overlap.py (repeats_020) again")
            hit[r["family"]][genome] = int(r["named_genes_hit"])
    genomes = list(genes_tables)
    out = {}
    for fam in set(host_main) | set(hit):
        hm = [host_main[fam].get(g, 0) for g in genomes]
        gh = [hit[fam].get(g, 0) for g in genomes]
        by_families, by_named = all(x >= host_genes for x in hm), all(x >= many_genes for x in gh)
        out[fam] = {"flag": by_families or by_named, "host_main": hm, "genes_hit": gh,
                    "category": host_cats[fam].most_common(1)[0][0] if host_cats.get(fam) else "MULTI",
                    "by": "both" if by_families and by_named else "gene families" if by_families else "named genes" if by_named else ""}
    return out


def best_hits(lines: Iterable[str]) -> Dict[str, tuple]:
    """query -> (subject, evalue, bitscore, coverage of the query or None) of its highest-bitscore hit, from blast
    -outfmt '6 qseqid sseqid evalue bitscore [qstart qend qlen]'."""
    best: Dict[str, tuple] = {}
    for line in lines:
        f = line.split()
        if len(f) < 4:
            continue
        q, s, e, b = f[0], f[1], float(f[2]), float(f[3])
        cov = (abs(int(f[5]) - int(f[4])) + 1) / int(f[6]) if len(f) >= 7 else None
        if q not in best or b > best[q][2]:
            best[q] = (s, e, b, cov)
    return best


def stale_files(work: Path) -> List[Path]:
    """Files left by the first version of this filter, whose host protein set was not cleaned."""
    return [work / n for n in OLD_METHOD_FILES if (work / n).exists()] + sorted(work.glob("host_db.*"))


def evidence_b(host: Dict[str, tuple], te: Dict[str, tuple], evalue: float, ratio: float) -> Set[str]:
    return {q for q, h in host.items() if h[1] <= evalue and h[2] >= ratio * (te[q][2] if q in te else 0.0)}


def b_strength(hit: tuple, product: str, min_coverage: float) -> str:
    """'strong' if the host match covers enough of the family or the host gene is a known multi-copy family, else 'weak'."""
    covered = hit[3] is not None and hit[3] >= min_coverage
    return "strong" if covered or HOST_FAMILIES.search(product or "") else "weak"


def decide(in_a: bool, b: str, repbase_support: bool) -> Tuple[str, str]:
    """b: 'strong', 'weak' or '' (no evidence B)."""
    in_b = b == "strong"
    if in_a and in_b:
        return "tag", "A and B"
    if in_a or in_b:
        line = "A only" if in_a else "B only"
        return ("keep", f"{line}; Repbase transposon support") if repbase_support else ("tag", f"{line}; no Repbase support")
    if b == "weak":
        return "keep", "weak B only (small host match, not a gene family): not tagged"
    return "keep", ""


def report_is_stale(report: Path) -> bool:
    """A report written by an earlier rule (without the tag_category column)."""
    if not report.exists():
        return False
    with open(report) as handle:
        return "tag_category" not in handle.readline().rstrip("\n").split("\t")


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
    stale = stale_files(work) + ([work / "filter_report.tsv"] if report_is_stale(work / "filter_report.tsv") else [])
    if stale:
        c.stop("the work folder holds results of an earlier filter method; remove them with:\n"
               "  rm " + " ".join(str(p) for p in stale))
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    work.mkdir(parents=True, exist_ok=True)
    families = {f.name: f for f in c.read_families(library)}
    names = {f"q{i:06d}": name for i, name in enumerate(families, start=1)}
    short_of = {v: k for k, v in names.items()}
    query, named_fa, host_fa = work / "families.fa", work / "host_named.fa", work / "host_clean.fa"
    if not query.exists():
        seqs = read_fasta(library)
        with open(query, "x") as handle:
            for short, name in names.items():
                handle.write(f">{short}\n{seqs[name]}\n")
        c.write_tsv(work / "query_names.tsv", ["query", "family"], names.items())
    def make_db(fasta: Path, db: Path) -> None:
        if not Path(str(db) + ".pdb").exists() and not Path(str(db) + ".phr").exists():
            c.run(c.container(args.sif, "makeblastdb", "-in", str(fasta), "-dbtype", "prot", "-out", str(db)),
                  True, work / "makeblastdb.log")

    make_db(args.te_proteins, work / "te_db")
    products: Dict[str, str] = {}
    if not named_fa.exists():
        print(f"translating the CDS of {args.proteins_from} ...", flush=True)
        with gzip.open(annotation, "rt") as handle:
            transcripts = read_transcripts(handle)
        proteins, skipped = host_proteins(transcripts, read_fasta(genome))
        with open(named_fa, "x") as handle:
            for gene, (_product, protein) in proteins.items():
                handle.write(f">{gene}\n{protein}\n")
        c.write_tsv(work / "host_named.tsv", ["gene", "product", "length"], [[g, p, len(s)] for g, (p, s) in proteins.items()])
        print(f"  {len(proteins):,} named host proteins; left out: " + ", ".join(f"{k} {v:,}" for k, v in skipped.items()))
    products = {r["gene"]: r["product"] for r in read_tsv(work / "host_named.tsv")}
    te_like_out = work / "blastp_host_te.tsv"
    if not te_like_out.exists():
        print("blastp of the host proteins against the transposon proteins ...", flush=True)
        c.run(c.container(args.sif, "blastp", "-query", str(named_fa), "-db", str(work / "te_db"), "-evalue",
                          f"{args.host_te_evalue:g}", "-outfmt", "6 qseqid sseqid evalue bitscore", "-max_target_seqs", "1",
                          "-max_hsps", "1", "-num_threads", str(args.threads), "-out", str(te_like_out)), True,
              work / "blastp_host_te.log")
    with open(te_like_out) as handle:
        te_like = {g: h for g, h in best_hits(handle).items() if h[1] <= args.host_te_evalue}
    if not host_fa.exists():
        named = read_fasta(named_fa)
        with open(host_fa, "x") as handle:
            for gene, protein in named.items():
                if gene not in te_like:
                    handle.write(f">{gene}\n{protein}\n")
        c.write_tsv(work / "host_left_out_te_like.tsv", ["gene", "product", "repeatpeps_match", "evalue"],
                    [[g, products.get(g, ""), h[0], f"{h[1]:g}"] for g, h in te_like.items()])
        print(f"  {len(te_like):,} named host proteins match a transposon protein and are left out; {len(named) - len(te_like):,} kept")
    make_db(host_fa, work / "host_clean_db")
    results = {}
    for tag, db, fmt in (("host_clean", work / "host_clean_db", "6 qseqid sseqid evalue bitscore qstart qend qlen"),
                         ("te", work / "te_db", "6 qseqid sseqid evalue bitscore")):
        out = work / f"blastx_{tag}.tsv"
        if not out.exists():
            print(f"blastx against the {'cleaned host' if tag == 'host_clean' else 'transposon'} proteins ...", flush=True)
            c.run(c.container(args.sif, "blastx", "-query", str(query), "-db", str(db), "-evalue", "1e-5",
                              "-outfmt", fmt, "-max_target_seqs", "5", "-max_hsps", "1",
                              "-num_threads", str(args.threads), "-out", str(out)), True, work / f"blastx_{tag}.log")
        with open(out) as handle:
            results["host" if tag == "host_clean" else "te"] = best_hits(handle)

    genes_tables = {g: read_tsv(t[0]) for g, t in tables.items()}
    families_tables = {g: read_tsv(t[1]) for g, t in tables.items()}
    a = evidence_a(genes_tables, families_tables, args.a_host_genes, args.a_genes)
    b = evidence_b(results["host"], results["te"], args.b_evalue, args.b_ratio)
    repbase = {r["family"]: r for r in read_tsv(repbase_table)}
    header = ["family", "label", "class", "length", "evidence_A", "evidence_A_by", "host_family_genes_mainly_masked",
              "named_genes_masked", "evidence_B", "evidence_B_strength", "best_host_protein", "best_host_product",
              "host_bitscore", "host_match_coverage_of_family", "best_te_protein", "te_bitscore",
              "repbase_transposon_support", "repbase_entry", "repbase_label", "repbase_coverage",
              "tag_category", "tagged_name", "decision", "reason"]
    rows, tagged = [], []
    tagged_alt = 0   # what the alternative rule (strong B by coverage only) would tag
    for name, fam in families.items():
        ev_a = a.get(fam.family, {})   # gene tables name families without their #class suffix
        in_a = ev_a.get("flag", False)
        host_hit, te_hit = results["host"].get(short_of[name]), results["te"].get(short_of[name])
        strength = b_strength(host_hit, products.get(host_hit[0], ""), args.b_min_coverage) if short_of[name] in b else ""
        if not (in_a or strength):
            continue
        rb = repbase.get(name, {})
        decision, reason = decide(in_a, strength, rb.get("transposon_support") == "yes")
        covered = bool(strength) and host_hit[3] is not None and host_hit[3] >= args.b_min_coverage
        alt, _ = decide(in_a, ("strong" if covered else "weak") if strength else "", rb.get("transposon_support") == "yes")
        tagged_alt += alt == "tag"
        product = products.get(host_hit[0], "") if host_hit else ""
        category = tag_category(product, strength, ev_a) if decision == "tag" else ""
        rows.append({"family": name, "label": fam.label, "class": fam.top_class, "length": fam.length,
                     "evidence_A": "yes" if in_a else "no", "evidence_A_by": ev_a.get("by", ""),
                     "host_family_genes_mainly_masked": "/".join(map(str, ev_a.get("host_main", []))),
                     "named_genes_masked": "/".join(map(str, ev_a.get("genes_hit", []))),
                     "evidence_B": "yes" if strength else "no", "evidence_B_strength": strength,
                     "best_host_protein": host_hit[0] if host_hit else "",
                     "best_host_product": products.get(host_hit[0], "") if host_hit else "",
                     "host_bitscore": f"{host_hit[2]:.0f}" if host_hit else "",
                     "host_match_coverage_of_family": f"{host_hit[3]:.2f}" if host_hit and host_hit[3] is not None else "",
                     "best_te_protein": te_hit[0] if te_hit else "", "te_bitscore": f"{te_hit[2]:.0f}" if te_hit else "",
                     "repbase_transposon_support": rb.get("transposon_support", ""), "repbase_entry": rb.get("repbase_entry", ""),
                     "repbase_label": rb.get("repbase_label", ""), "repbase_coverage": rb.get("coverage_of_family", ""),
                     "tag_category": category, "tagged_name": tagged_name(fam.family, category, fam.label) if category else "",
                     "decision": decision, "reason": reason})
        if decision == "tag":
            tagged.append(name)
    table_rows = [[r[h] for h in header] for r in rows]
    report = work / "filter_report.tsv"
    if not report.exists():
        c.write_tsv(report, header, table_rows)
    n_a = sum(1 for r in rows if r["evidence_A"] == "yes")
    n_strong = sum(1 for r in rows if r["evidence_B_strength"] == "strong")
    n_weak = sum(1 for r in rows if r["evidence_B_strength"] == "weak")
    both = sum(1 for r in rows if r["evidence_A"] == "yes" and r["evidence_B_strength"] == "strong")
    a_by = Counter(r["evidence_A_by"] for r in rows if r["evidence_A"] == "yes")
    print(f"\n{args.version}: {len(families):,} families; evidence A {n_a} (" + ", ".join(f"by {k} {v}" for k, v in a_by.most_common())
          + f"); evidence B strong {n_strong}, weak {n_weak}; A and strong B {both}")
    for reason, n in Counter(r["reason"] for r in rows).most_common():
        print(f"  {reason}: {n}")
    print(f"tagged: {len(tagged)} families; by class: {', '.join(f'{k} {v}' for k, v in Counter(families[f].top_class for f in tagged).most_common())}")
    print("  by category: " + ", ".join(f"{k} {v}" for k, v in Counter(r["tag_category"] for r in rows if r["decision"] == "tag").most_common()))
    print(f"  the alternative rule (strong B by coverage only, no gene-family exception) would tag {tagged_alt}")

    # how much the tagged families mask
    bare = {families[f].family for f in tagged}
    impact = {}
    for g in args.genomes:
        classes = read_tsv(genes_dir / f"{args.version}_{g}_classes.tsv")
        total = next(int(r["cds_bp"]) for r in classes if r["gene_set"] == "genes without TE products" and r["class"] == "ALL")
        masked = next(int(r["masked_bp"]) for r in classes if r["gene_set"] == "genes without TE products" and r["class"] == "ALL")
        bp = sum(int(r["cds_bp_masked_in_other_genes"]) for r in families_tables[g] if r["family"] in bare)
        impact[g] = {"cds_bp": total, "cds_masked_bp": masked, "cds_bp_masked_by_tagged": bp}
        print(f"  {g}: tagged families mask {100 * bp / total:.2f}% of CDS, of the {100 * masked / total:.2f}% masked by all")
    bed = c.version_folder(args.version) / "masked" / args.impact_genome / f"{args.impact_genome}.bed"
    impact_fasta = c.genome_path(args.impact_genome, c.panel_accession(args.impact_genome))
    if bed.exists() and impact_fasta.exists():
        print(f"reading the masking of {args.impact_genome} ...", flush=True)
        stretches = la.one_label_per_base(la.read_labelled_bed(bed))
        genome_bp = sum(c.fasta_lengths(impact_fasta).values())
        all_bp = sum(r[2] - r[1] for r in stretches)
        rem_bp = sum(r[2] - r[1] for r in stretches if r[3] in bare)
        impact[args.impact_genome] = {"genome_bp": genome_bp, "masked_bp": all_bp, "masked_by_tagged_bp": rem_bp}
        print(f"  {args.impact_genome}: tagged families mask {100 * rem_bp / genome_bp:.2f}% of the genome, "
              f"of the {100 * all_bp / genome_bp:.2f}% masked by all")
    else:
        print(f"  (no masking of {args.impact_genome} by {args.version}, or no genome file; genome masking not estimated)")
    print(f"table: {report}")
    if not args.make_version:
        print("report only: no version written. Review the table, then add --make-version.")
        return

    library_dir = target / "library"
    library_dir.mkdir(parents=True)
    seqs = read_fasta(library)
    new_names = {r["family"]: r["tagged_name"] for r in rows if r["decision"] == "tag"}
    out_lib = library_dir / f"{args.new_version}-families.prefixed.fa"
    with open(out_lib, "x") as handle:
        for name in families:
            s = seqs[name]
            handle.write(f">{new_names.get(name, name)}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
    c.write_tsv(library_dir / "tag_table.tsv", header, table_rows)
    record = {"version": args.new_version, "kind": "tagged library", "source_version": args.version,
              "source_library": {"path": str(library), "sha256": c.sha256(library)},
              "te_proteins": {"path": str(args.te_proteins), "sha256": c.sha256(args.te_proteins)},
              "host_proteins_from": args.proteins_from, "annotation": str(annotation),
              "gene_tables": {g: [str(p) for p in t] for g, t in tables.items()}, "repbase_table": str(repbase_table),
              "thresholds": {"a_host_genes": args.a_host_genes, "a_genes": args.a_genes, "b_evalue": args.b_evalue,
                             "b_ratio": args.b_ratio, "host_te_evalue": args.host_te_evalue, "b_min_coverage": args.b_min_coverage},
              "estimated_cost": impact,
              "host_set": "named proteins; uncharacterized and transposon products left out; RepeatPeps matches left out (blastp)",
              "host_families_pattern": HOST_FAMILIES.pattern, "rule": "A and strong B: tag; one line: tag unless Repbase transposon support; weak B alone: not tagged, listed",
              "tag_format": "<family>_hostgene-<CATEGORY>#<class>", "categories": [n for n, _ in CATEGORIES] + ["OTHER", "MULTI"],
              "counts": {"families": len(families), "tagged": len(tagged), "tagged_by_alternative_rule": tagged_alt, "evidence_A": n_a,
                         "evidence_B_strong": n_strong, "evidence_B_weak": n_weak, "both": both,
                         "by_reason": dict(Counter(r["reason"] for r in rows))},
              "written": f"{datetime.now():%Y-%m-%d %H:%M}",
              "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()}
    with open(library_dir / "filter_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"written: {library_dir}/ ({out_lib.name}, tag_table.tsv, filter_record.json); {len(families):,} families, {len(tagged)} tagged")


if __name__ == "__main__":
    main()
