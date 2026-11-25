#!/usr/bin/env python3
"""Download a reference genome into scratch (spec section 6, 2026-11-23).

A reference genome is a published assembly used only to compare our libraries with a published TE annotation of that
same assembly (like for like). It is never used to build or curate a library, nor for alignment, and it is not part of
the shared genome store. Each is listed in references/reference_genomes.tsv in this repository (its NCBI accession and
assembly name) and lives in <scratch>/references/<id>/.

From the accession's folder on NCBI's FTP server
(https://ftp.ncbi.nlm.nih.gov/genomes/all/<GCF|GCA>/<nnn>/<nnn>/<nnn>/<accession>_<assembly name>/) it downloads
md5checksums.txt and then, each checked against it:
  <prefix>_genomic.fna.gz        kept as genome.fna.gz
  <prefix>_genomic.gff.gz        kept as annotation.gff.gz
  <prefix>_assembly_report.txt   kept as assembly_report.txt
  <prefix>_rm.out.gz             kept as ncbi_rm.out.gz, where NCBI provides it: NCBI's own RepeatMasker run on this
                                 assembly, a second published annotation of it
It writes provenance.json (URLs, md5 sums, sizes, sequence count and length, date). Never overwrites.

Usage (from the repository root, in ga-repeats; mzebra needs internet access):
  python scripts/fetch_reference.py --id AC_fAstCal1_2          # dry run: lists the URLs
  python scripts/fetch_reference.py --id AC_fAstCal1_2 --go
"""
import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

FILES = [("_genomic.fna.gz", "genome.fna.gz", True), ("_genomic.gff.gz", "annotation.gff.gz", True),
         ("_assembly_report.txt", "assembly_report.txt", True), ("_rm.out.gz", "ncbi_rm.out.gz", False)]


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--id", required=True, help="an id in references/reference_genomes.tsv, e.g. AC_fAstCal1_2")
    parser.add_argument("--base-url", default="https://ftp.ncbi.nlm.nih.gov/genomes/all/",
                        help="NCBI's genomes/all/ folder (a local file:// copy for tests)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def folder_url(base: str, accession: str, assembly_name: str) -> str:
    """NCBI's folder for an assembly: GCF_900246225.1 + fAstCal1.2 -> GCF/900/246/225/GCF_900246225.1_fAstCal1.2/."""
    m = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", accession)
    if not m:
        raise ValueError(f"not an NCBI assembly accession: {accession}")
    folder = f"{accession}_{assembly_name.replace(' ', '_')}"
    return f"{base.rstrip('/')}/{m.group(1)}/{m.group(2)}/{m.group(3)}/{m.group(4)}/{folder}/"


def parse_md5(text: str) -> Dict[str, str]:
    """md5checksums.txt -> file name -> md5."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1].lstrip("./").lstrip("*")] = parts[0].lower()
    return out


def download(url: str, dst: Path) -> str:
    """Stream a URL into a new file; return its md5."""
    digest = hashlib.md5()
    with urllib.request.urlopen(url, timeout=120) as response, open(dst, "xb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            digest.update(chunk)
            out.write(chunk)
    return digest.hexdigest()


def main(argv=None) -> None:
    args = parse_args(argv)
    refs = c.reference_genomes()
    if args.id not in refs:
        c.stop(f"{args.id} is not in {c.REFERENCES_FILE}")
    entry = refs[args.id]
    url = folder_url(args.base_url, entry["accession"], entry["assembly_name"])
    prefix = url.rstrip("/").rsplit("/", 1)[1]
    target = c.reference_folder(args.id)
    if target.exists():
        c.stop(f"{target} exists; it is never overwritten. To fetch again: rm -r {target}")
    print(f"{args.id}: {entry['species']} {entry['assembly_name']} ({entry['accession']})\nfrom: {url}\nto:   {target}")
    for suffix, local, required in FILES:
        print(f"  {prefix}{suffix} -> {local}{'' if required else ' (if NCBI provides it)'}")
    if not args.go:
        print("dry run: nothing was downloaded. Add --go to fetch.")
        return

    target.mkdir(parents=True)
    with urllib.request.urlopen(url + "md5checksums.txt", timeout=120) as response:
        md5s = parse_md5(response.read().decode())
    record = {"id": args.id, "species": entry["species"], "accession": entry["accession"],
              "assembly_name": entry["assembly_name"], "purpose": entry["purpose"], "source_folder": url,
              "fetched": f"{datetime.now():%Y-%m-%d %H:%M:%S}", "files": {}}
    for suffix, local, required in FILES:
        name = prefix + suffix
        if name not in md5s:
            if required:
                c.stop(f"{name} is not listed in {url}md5checksums.txt")
            print(f"  {name}: not provided by NCBI")
            record["files"][local] = None
            continue
        print(f"  downloading {name} ...", flush=True)
        try:
            got = download(url + name, target / local)
        except urllib.error.URLError as error:
            c.stop(f"could not download {url + name}: {error}")
        if got != md5s[name]:
            c.stop(f"md5 mismatch for {name}: {got} downloaded, {md5s[name]} listed; remove {target} and try again")
        record["files"][local] = {"url": url + name, "md5": got, "bytes": (target / local).stat().st_size}
    lengths = c.fasta_lengths(target / "genome.fna.gz")
    record["sequences"] = len(lengths)
    record["total_bp"] = sum(lengths.values())
    with open(target / "provenance.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"\n{args.id}: {len(lengths):,} sequences, {record['total_bp']:,} bp; md5 sums match NCBI's")
    print(f"written: {target}/ (genome.fna.gz, annotation.gff.gz, assembly_report.txt"
          + (", ncbi_rm.out.gz" if record["files"].get("ncbi_rm.out.gz") else "") + ", provenance.json)")
    print(f"next: python scripts/mask_genomes.py --version RI_v5 --assemblies {args.id} --go")


if __name__ == "__main__":
    main()
