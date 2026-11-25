#!/usr/bin/env python3
"""Ensembl's soft-masked genome of a reference assembly (spec section 6, 2026-11-24).

Ensembl's soft-masked genome writes in lower case every base its own repeat pipeline masked. For fAstCal1.2 this is a
published masking of the same assembly Quah et al. 2025 used (Ensembl release 103), independent of our libraries; it
carries no classes, so it serves as a reference for which bases are repeats, not for labels.

Downloads <Species>.<assembly name>.dna_sm.toplevel.fa.gz from https://ftp.ensembl.org/pub/release-<N>/fasta/<species>/dna/
into <scratch>/references/<id>/ensembl_dna_sm.toplevel.fa.gz, checks it against Ensembl's CHECKSUMS file in that folder
(Ensembl publishes the BSD checksum of each file, as the Unix 'sum' command computes it), and writes
ensembl_provenance.json. The reference folder must exist (fetch_reference.py); nothing in it is overwritten.

Usage (from the repository root):
  python scripts/fetch_ensembl.py --id AC_fAstCal1_2 --release 103 --species astatotilapia_calliptera          # dry run
  python scripts/fetch_ensembl.py --id AC_fAstCal1_2 --release 103 --species astatotilapia_calliptera --go
"""
import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import fetch_reference as fr  # noqa: E402

LOCAL = "ensembl_dna_sm.toplevel.fa.gz"


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--id", required=True, help="a reference genome id, e.g. AC_fAstCal1_2")
    parser.add_argument("--release", type=int, required=True, help="Ensembl release, e.g. 103")
    parser.add_argument("--species", required=True, help="Ensembl's species folder, e.g. astatotilapia_calliptera")
    parser.add_argument("--base-url", default="https://ftp.ensembl.org/pub/", help="Ensembl's FTP root (file:// for tests)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def bsd_sum_python(path: Path) -> Tuple[int, int]:
    """The BSD checksum and 1 KiB block count, as 'sum -r' prints them."""
    checksum, size = 0, 0
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1 << 20)
            if not chunk:
                break
            size += len(chunk)
            for byte in chunk:
                checksum = (checksum >> 1) + ((checksum & 1) << 15)
                checksum = (checksum + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def bsd_sum(path: Path) -> Tuple[int, int]:
    """'sum -r' where available (fast), otherwise the same algorithm in Python."""
    if shutil.which("sum"):
        out = subprocess.run(["sum", "-r", str(path)], capture_output=True, text=True)
        if out.returncode == 0 and len(out.stdout.split()) >= 2:
            a, b = out.stdout.split()[:2]
            return int(a), int(b)
    return bsd_sum_python(path)


def parse_checksums(text: str) -> Dict[str, Tuple[int, int]]:
    """Ensembl's CHECKSUMS ('<sum> <blocks> <file>') -> file -> (sum, blocks)."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            out[parts[2]] = (int(parts[0]), int(parts[1]))
    return out


def main(argv=None) -> None:
    args = parse_args(argv)
    refs = c.reference_genomes()
    if args.id not in refs:
        c.stop(f"{args.id} is not in {c.REFERENCES_FILE}")
    folder = c.reference_folder(args.id)
    if not (folder / "genome.fna.gz").exists():
        c.stop(f"{folder / 'genome.fna.gz'} not found: run fetch_reference.py --id {args.id} first")
    target = folder / LOCAL
    if target.exists():
        c.stop(f"{target} exists; it is never overwritten")
    assembly = refs[args.id]["assembly_name"]
    url = f"{args.base_url.rstrip('/')}/release-{args.release}/fasta/{args.species}/dna/"
    name = f"{args.species[:1].upper()}{args.species[1:]}.{assembly}.dna_sm.toplevel.fa.gz"
    print(f"{args.id}: Ensembl release {args.release}, {name}\nfrom: {url}\nto:   {target}")
    if not args.go:
        print("dry run: nothing was downloaded. Add --go to fetch.")
        return

    with urllib.request.urlopen(url + "CHECKSUMS", timeout=120) as response:
        sums = parse_checksums(response.read().decode())
    if name not in sums:
        c.stop(f"{name} is not listed in {url}CHECKSUMS")
    print(f"  downloading {name} ...", flush=True)
    fr.download(url + name, target)
    got = bsd_sum(target)
    if got != sums[name]:
        target.unlink()
        c.stop(f"checksum mismatch for {name}: {got} downloaded, {sums[name]} listed; the download was removed")
    record = {"id": args.id, "source": "Ensembl", "release": args.release, "species": args.species, "assembly": assembly,
              "url": url + name, "bsd_sum": list(got), "bytes": target.stat().st_size,
              "fetched": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
              "use": "which bases Ensembl's repeat pipeline masked (lower case); no classes"}
    with open(folder / "ensembl_provenance.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(f"\n{name}: checksum matches Ensembl's ({got[0]} {got[1]})\nwritten: {target}")
    print(f"next: python scripts/ensembl_masking.py --version <version> --assembly {args.id} --go")


if __name__ == "__main__":
    main()
