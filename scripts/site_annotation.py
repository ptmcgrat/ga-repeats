#!/usr/bin/env python3
"""The insertion sites annotated with Repbase: an independent annotation beside ours (question 2).

insertion_sites.py found the published insertions and our library's annotation of them (findings B5): MC_MC1m's mitfa
insertion masked as LTR where Carleton et al. 2020 report mainly Rex1 (LINE), MZ_NMK2024a's tbx2a insertion masked as
LTR/Gypsy where they report hAT-8_AFC (DNA). This annotates the same sequences with Repbase, independently of our
library: every site with sequence between the flanks (mitfa, rx1) and every tbx2a locus, from insertion_sites.py's
table, is cut from its genome and annotated by RepeatMasker (-engine rmblast -no_is -nolow) against the classified
Repbase (all entries, as Carleton et al.'s CENSOR search used all of Repbase; --cichlid-only for the curated cichlid
entries). Each site's Repbase annotation is printed beside ours (from insertion_sites.py's annotation table).
With --with-version <version>, the sites are annotated with that version's whole library instead (for instance our
library combined with Repbase), to see which family wins each site (repeats_029).

Usage (from the repository root):
  python scripts/site_annotation.py --version RI_v5          # dry run
  python scripts/site_annotation.py --version RI_v5 --go

Writes, never overwriting, <scratch>/comparisons/insertion_sites/<version>_site_repbase.tsv (one row per Repbase match in
a site). Work in <scratch>/<version>/site_annotation/.
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

Site = Tuple[str, str, str, int, int]   # assembly, locus, sequence, start, end (0-based half-open)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the version whose insertion_sites.py tables to use, e.g. RI_v5")
    parser.add_argument("--cichlid-only", action="store_true", help="annotate with the curated cichlid entries only")
    parser.add_argument("--with-version", help="annotate with this version's whole library instead of Repbase "
                                               "(e.g. a library combined with Repbase: which family wins at each site)")
    parser.add_argument("--species", nargs="+", default=["Cichlidae", "Oreochromis_niloticus"])
    parser.add_argument("--release", default="20181026")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def read_sites(lines) -> List[Site]:
    """Sites with sequence: flank pairs with bases between them, and tbx2a loci."""
    out = []
    for r in csv.DictReader(lines, delimiter="\t"):
        if r.get("start") and r.get("end") and int(r["end"]) > int(r["start"]):
            out.append((r["assembly"], r["locus"], r["seq"], int(r["start"]), int(r["end"])))
    return out


def cut_sites(fasta: Path, sites: List[Site]) -> Dict[Site, str]:
    """The sequence of each site (streaming the genome once)."""
    wanted: Dict[str, List[Site]] = defaultdict(list)
    for s in sites:
        wanted[s[2]].append(s)
    out, name, parts = {}, None, []
    opener = gzip.open if str(fasta).endswith(".gz") else open

    def flush():
        if name in wanted:
            seq = "".join(parts)
            for s in wanted[name]:
                out[s] = seq[s[3]:s[4]]
    with opener(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                name, parts = line[1:].split()[0], []
            elif name in wanted:
                parts.append(line.strip())
        flush()
    return out


def parse_hits(lines) -> List[Tuple[str, int, int, str, str, float]]:
    """RepeatMasker .out -> (site name, begin, end, entry, entry label, identity), 1-based within the site."""
    out = []
    for line in lines:
        f = line.split()
        if len(f) >= 15 and f[0].isdigit():
            out.append((f[4], int(f[5]), int(f[6]), f[9], f[10], 100.0 - float(f[1])))
    return out


def main(argv=None) -> None:
    args = parse_args(argv)
    sites_file = args.scratch / "comparisons" / "insertion_sites" / f"{args.version}_insertion_sites.tsv"
    ours_file = args.scratch / "comparisons" / "insertion_sites" / f"{args.version}_insertion_sites_annotations.tsv"
    lib_dir = args.scratch / "repbase" / args.release
    repbase_fa, entries_file = lib_dir / "repbase.fa", lib_dir / "repbase_entries.tsv"
    for f in (sites_file, ours_file, repbase_fa, entries_file):
        if not f.exists():
            c.stop(f"{f} not found: run insertion_sites.py and repbase_match.py first")
    tag = f"lib_{args.with_version}" if args.with_version else ("repbase_cichlid" if args.cichlid_only else "repbase")
    target = args.scratch / "comparisons" / "insertion_sites" / f"{args.version}_site_{tag}.tsv"
    if target.exists():
        c.stop(f"{target} exists; results are never overwritten")
    work = args.scratch / args.version / {"repbase": "site_annotation", "repbase_cichlid": "site_annotation_cichlid"}.get(
        tag, f"site_annotation_{tag}")
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    with open(sites_file) as handle:
        sites = read_sites(handle)
    print(f"{args.version}: {len(sites)} sites with sequence, annotated with "
          + (f"{args.with_version}'s whole library" if args.with_version else
             "the curated cichlid entries" if args.cichlid_only else "all of Repbase") + f"\nresult: {target}")
    for s in sites:
        print(f"  {s[0]} {s[1]} {s[2]}:{s[3] + 1:,}-{s[4]:,} ({s[4] - s[3]:,} bp)")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    work.mkdir(parents=True)
    lib = c.library_files(args.with_version)[2] if args.with_version else repbase_fa
    if args.cichlid_only and not args.with_version:
        with open(entries_file) as handle:
            wanted = {r["name"] for r in csv.DictReader(handle, delimiter="\t") if r["species"] in set(args.species)}
        lib, keep = work / "curated.fa", False
        with open(repbase_fa) as src, open(lib, "x") as dst:
            for line in src:
                if line.startswith(">"):
                    keep = line[1:].split("#", 1)[0].strip() in wanted
                if keep:
                    dst.write(line)
    seqs: Dict[Site, str] = {}
    for assembly in sorted({s[0] for s in sites}):
        mine = [s for s in sites if s[0] == assembly]
        seqs.update(cut_sites(c.genome_path(assembly, c.panel_accession(assembly)), mine))
    names = {}
    with open(work / "sites.fa", "x") as handle:
        for i, s in enumerate(sites, start=1):
            names[f"site{i}"] = s
            handle.write(f">site{i}\n{seqs.get(s, '')}\n")
    pa = max(1, args.threads // c.RMBLAST_CORES_PER_JOB)
    c.run(c.container(args.sif, "RepeatMasker", "-lib", str(lib), "-engine", "rmblast", "-pa", str(pa), "-no_is", "-nolow",
                      "-dir", str(work), str(work / "sites.fa")), True, work / "repeatmasker.log")
    out_file = work / "sites.fa.out"
    hits = parse_hits(open(out_file)) if out_file.exists() else []
    with open(ours_file) as handle:
        ours = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for name, s in names.items():
        assembly, locus, seq, start, end = s
        mine = [h for h in hits if h[0] == name]
        theirs = [r for r in ours if r["assembly"] == assembly and r["locus"] == locus]
        print(f"\n{assembly} {locus} {seq}:{start + 1:,}-{end:,} ({end - start:,} bp)")
        print(f"  {'library' if args.with_version else 'Repbase'}: " + ("; ".join(f"{h[1]}-{h[2]} {h[3]} ({h[4]}, {h[5]:.0f}%)" for h in mine) or "no match"))
        print("  ours:    " + ("; ".join(f"{int(r['start']) - start + 1}-{int(r['end']) - start} {r['family']} ({r['label']})"
                                         for r in theirs) or "nothing masked"))
        for h in mine:
            rows.append([assembly, locus, seq, start, end, h[1], h[2], h[3], h[4], f"{h[5]:.1f}"])
    c.write_tsv(target, ["assembly", "locus", "seq", "site_start", "site_end", "begin_in_site", "end_in_site", "entry",
                         "entry_label", "identity"], rows)
    print(f"\ntable: {target}")


if __name__ == "__main__":
    main()
