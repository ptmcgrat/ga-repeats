"""Shared helpers for the ga-repeats pilot scripts (repeats_009).

Python 3.9 or newer (the common library's minimum). Standard library only, plus
ga_common for data locations. Facts used here were read from the installed tools'
own output (facts_for_script.txt, 2026-11):
  - TETools 2.01 image: RepeatModeler 2.0.9, RepeatMasker 4.2.6, CD-HIT 4.8.1
  - cd-hit-est inside the image: /opt/cd-hit/cd-hit-est
  - RepeatMasker -pa: each job uses 4 cores with RMBlast (RepeatMasker -help, 1.3)
  - library names: >repeatname#class/subclass (RepeatMasker -help, 1.1 -lib)
  - .out columns: RepeatMasker -help, section 3.1
  - .clstr member lines with -G 0 carry alignment coordinates (real output, 2026-11; repeats_010)
"""
from __future__ import annotations

import argparse
import csv
import functools
import gzip
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Tuple

PANEL_VERSION = "P1"
CDHIT_EST = "/opt/cd-hit/cd-hit-est"
DEFAULT_SIF = Path.home() / "containers" / "dfam-tetools-2.01.sif"
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-repeats"
RMBLAST_CORES_PER_JOB = 4


# ---------------------------------------------------------------- arguments

def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Options every script takes (spec section 9: --threads and --max-jobs)."""
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("GA_MAX_THREADS", "16")),
                        help="total threads for this script (env GA_MAX_THREADS; default 16)")
    parser.add_argument("--max-jobs", type=int,
                        default=int(os.environ.get("GA_MAX_JOBS", "1")),
                        help="jobs at once (env GA_MAX_JOBS; default 1)")
    parser.add_argument("--sif", type=Path,
                        default=Path(os.environ.get("GA_REPEATS_SIF", str(DEFAULT_SIF))),
                        help="TETools image (default ~/containers/dfam-tetools-2.01.sif)")
    parser.add_argument("--scratch", type=Path,
                        default=Path(os.environ.get("GA_REPEATS_SCRATCH", str(DEFAULT_SCRATCH))),
                        help="working space (default ~/scratch/ga-repeats)")
    parser.add_argument("--go", action="store_true",
                        help="do it; without --go the script only prints what it would do")


def stop(message: str) -> "NoReturn":  # noqa: F821
    print(f"STOP: {message}", file=sys.stderr)
    raise SystemExit(1)


# ---------------------------------------------------------------- locations

@functools.lru_cache(maxsize=1)
def _paths():
    try:
        from ga_common.paths import Paths
    except ImportError:
        stop("ga_common is not importable; run `conda activate ga-repeats` first (spec section 2)")
    return Paths.load()


def genome_path(assembly_id: str, accession_version: str) -> Path:
    return Path(str(_paths().local("genome_fasta", assembly_id=assembly_id,
                                   accession_version=accession_version)))


def version_folder(version: str) -> Path:
    """The version folder is the folder that holds the version's MANIFEST.yaml.

    library/ and masked/<assembly_id>/ inside it follow the layout in spec section 3;
    ga_common names no separate keys for them that we know of.
    """
    return Path(str(_paths().local("repeat_manifest", version=version))).parent


def config_dir() -> Path:
    """Same rule as ga_common (coordinator spec section 8)."""
    env = os.environ.get("GA_CONFIG_DIR")
    return Path(env) if env else Path.home() / "GenomeAlignments" / "ga-coordinator" / "config"


def panel_accession(assembly_id: str, panel: str = PANEL_VERSION) -> str:
    panel_file = config_dir() / f"genome_panel_{panel}.tsv"
    if not panel_file.exists():
        stop(f"panel file not found: {panel_file}")
    with open(panel_file, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("id") == assembly_id:
                return row["accession_genbank"]
    stop(f"{assembly_id} is not in {panel_file}")


def library_files(version: str) -> Tuple[str, Path, Path]:
    """(build assembly id, original families file, prefixed families file) of a version.

    A merged version (merge_libraries.py) has no single build assembly: its library is
    <version>-families.prefixed.fa, its names already carry their sources, and its
    library/ holds merge_record.json. For it, the build id is the version and the
    "original" file is the merged library itself.
    """
    lib = version_folder(version) / "library"
    prefixed = sorted(lib.glob("*-families.prefixed.fa"))
    if len(prefixed) != 1:
        stop(f"expected one *-families.prefixed.fa in {lib}, found {len(prefixed)}")
    build_id = prefixed[0].name[: -len("-families.prefixed.fa")]
    if build_id == version and (lib / "merge_record.json").exists():
        return build_id, prefixed[0], prefixed[0]
    original = lib / f"{build_id}-families.fa"
    if not original.exists():
        stop(f"missing {original}")
    return build_id, original, prefixed[0]


# ---------------------------------------------------------------- sequences

def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta_lengths(path: Path) -> Dict[str, int]:
    """Sequence name (first word of the header) -> length, in file order."""
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    size = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = size
                words = line[1:].split(None, 1)
                name = words[0] if words else ""
                if name in lengths:
                    raise ValueError(f"duplicate sequence name {name} in {path}")
                size = 0
            else:
                size += len(line.rstrip("\r\n"))
    if name is not None:
        lengths[name] = size
    return lengths


def decompress(src: Path, dst: Path) -> None:
    """Write a plain copy of a (BGZF-)gzipped file; refuses to overwrite."""
    with gzip.open(src, "rb") as fin, open(dst, "xb") as fout:
        shutil.copyfileobj(fin, fout, 1024 * 1024)


class Family(NamedTuple):
    name: str          # full first word, e.g. MZ_MZ4f_rnd-1_family-23#LTR/Gypsy
    family: str        # part before '#'
    label: str         # part after '#', e.g. LTR/Gypsy ('Unknown' if missing)
    top_class: str     # part of the label before '/', e.g. LTR
    length: int


def split_name(name: str) -> Tuple[str, str]:
    family, _, label = name.partition("#")
    return family, (label or "Unknown")


def top_class(label: str) -> str:
    return label.split("/", 1)[0] or "Unknown"


def read_families(path: Path) -> List[Family]:
    out = []
    for name, length in fasta_lengths(path).items():
        family, label = split_name(name)
        out.append(Family(name, family, label, top_class(label), length))
    return out


def source_of(name: str, prefixes: Sequence[str]) -> str:
    """Which build assembly a prefixed family name comes from (longest prefix wins)."""
    for prefix in sorted(prefixes, key=len, reverse=True):
        if name.startswith(prefix + "_"):
            return prefix
    raise ValueError(f"{name} carries none of the prefixes {list(prefixes)}")


# ---------------------------------------------------------------- cd-hit

class ClusterMember(NamedTuple):
    cluster: int
    name: str
    length: int
    representative: bool
    identity: Optional[float]  # percent identity to the representative; None for it


# Member lines: "at +/96.27%", or with local alignment (-G 0, CD-HIT 4.8.1) "at 1:161:8923:9083/+/96.27%"
# (query start:end, representative start:end, strand, identity).
_CLSTR_LINE = re.compile(r"^\d+\t(\d+)nt, >(.+?)\.\.\. "
                         r"(?:(\*)|at (?:\d+:\d+:\d+:\d+/)?(?:[+-]/)?([\d.]+)%)\s*$")


def parse_clstr(path: Path) -> List[ClusterMember]:
    members: List[ClusterMember] = []
    cluster = -1
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">Cluster "):
                cluster = int(line.split()[1])
                continue
            if not line.strip():
                continue
            match = _CLSTR_LINE.match(line)
            if not match:
                raise ValueError(f"unrecognised .clstr line in {path}: {line!r}")
            length, name, star, ident = match.groups()
            members.append(ClusterMember(cluster, name, int(length), star == "*",
                                         None if star else float(ident)))
    return members


_CLSTR_ALIGN = re.compile(r"^\d+\t\d+nt, >(.+?)\.\.\. at (\d+):(\d+):(\d+):(\d+)/([+-])/[\d.]+%\s*$")


def clstr_alignments(path: Path) -> Dict[str, Tuple[int, int, int, int, str]]:
    """Member name -> (start, end on the member, start, end on the representative, strand),
    from a cd-hit-est .clstr written with local alignment (-G 0). Representatives are absent."""
    out: Dict[str, Tuple[int, int, int, int, str]] = {}
    with open(path) as handle:
        for line in handle:
            match = _CLSTR_ALIGN.match(line.rstrip("\n"))
            if match:
                name, qs, qe, rs, re_, strand = match.groups()
                out[name] = (int(qs), int(qe), int(rs), int(re_), strand)
    return out


def word_size(identity: float) -> int:
    """cd-hit-est word length for an identity threshold (CD-HIT user's guide)."""
    if identity >= 0.95:
        return 10
    if identity >= 0.90:
        return 8
    if identity >= 0.88:
        return 7
    if identity >= 0.85:
        return 6
    if identity >= 0.80:
        return 5
    if identity >= 0.75:
        return 4
    raise ValueError("cd-hit-est cannot cluster nucleotide sequences below 0.75 identity")


# ---------------------------------------------------------------- RepeatMasker .out

class Hit(NamedTuple):
    score: int
    divergence: float
    deletion: float
    insertion: float
    query: str
    begin: int          # 1-based, inclusive
    end: int
    query_left: int
    strand: str         # '+' or 'C'
    repeat: str
    label: str          # class/family column
    repeat_id: str


def parse_out(path: Path) -> Iterator[Hit]:
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue
            yield Hit(int(fields[0]), float(fields[1]), float(fields[2]), float(fields[3]),
                      fields[4], int(fields[5]), int(fields[6]),
                      int(fields[7].strip("()")), fields[8], fields[9], fields[10], fields[14])


def hit_to_bed(hit: Hit) -> List[str]:
    """Pilot BED (BED6+4): the contract's exact fields are still open (spec section 7)."""
    strand = "-" if hit.strand == "C" else "+"
    return [hit.query, str(hit.begin - 1), str(hit.end), hit.repeat, "0", strand,
            hit.label, f"{hit.divergence:.1f}", str(hit.score), hit.repeat_id]


def check_bed(bed_rows: Iterable[Sequence[str]], lengths: Dict[str, int]) -> List[str]:
    """Problems that would stop a masked genome being regenerated (spec section 5)."""
    problems = []
    for row in bed_rows:
        name, start, end = row[0], int(row[1]), int(row[2])
        if name not in lengths:
            problems.append(f"{name} is not a sequence of the genome")
        elif not (0 <= start < end <= lengths[name]):
            problems.append(f"{name}:{start}-{end} lies outside 0-{lengths[name]}")
        if len(problems) >= 20:
            problems.append("... (further problems not listed)")
            break
    return problems


def read_bed(path: Path) -> List[Tuple[str, int, int, str, str]]:
    """(sequence, start, end, family, top class) from a pilot BED written by mask_genomes.py."""
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f[3], top_class(f[6])))
    return rows


def merge_by_sequence(intervals: Iterable[Tuple[str, int, int]]) -> Dict[str, List[Tuple[int, int]]]:
    """Half-open intervals merged per sequence: overlapping or touching spans become one."""
    by_sequence: Dict[str, List[Tuple[int, int]]] = {}
    for name, start, end in intervals:
        by_sequence.setdefault(name, []).append((start, end))
    merged: Dict[str, List[Tuple[int, int]]] = {}
    for name, spans in by_sequence.items():
        spans.sort()
        out = [list(spans[0])]
        for start, end in spans[1:]:
            if start > out[-1][1]:
                out.append([start, end])
            else:
                out[-1][1] = max(out[-1][1], end)
        merged[name] = [(a, b) for a, b in out]
    return merged


def merged_length(merged: Dict[str, List[Tuple[int, int]]]) -> int:
    return sum(end - start for spans in merged.values() for start, end in spans)


def subtract_length(a: Dict[str, List[Tuple[int, int]]], b: Dict[str, List[Tuple[int, int]]]) -> int:
    """Base pairs covered by a and not by b (both merged per sequence)."""
    total = 0
    for name, a_spans in a.items():
        b_spans = b.get(name, [])
        j = 0
        for start, end in a_spans:
            while j < len(b_spans) and b_spans[j][1] <= start:
                j += 1
            current = start
            k = j
            while k < len(b_spans) and b_spans[k][0] < end:
                if b_spans[k][0] > current:
                    total += b_spans[k][0] - current
                current = max(current, b_spans[k][1])
                if current >= end:
                    break
                k += 1
            if current < end:
                total += end - current
    return total


def intersect_length(a: Dict[str, List[Tuple[int, int]]], b: Dict[str, List[Tuple[int, int]]]) -> int:
    """Base pairs covered by both a and b (both merged per sequence)."""
    total = 0
    for name, a_spans in a.items():
        b_spans = b.get(name)
        if not b_spans:
            continue
        i = j = 0
        while i < len(a_spans) and j < len(b_spans):
            start = max(a_spans[i][0], b_spans[j][0])
            end = min(a_spans[i][1], b_spans[j][1])
            if end > start:
                total += end - start
            if a_spans[i][1] < b_spans[j][1]:
                i += 1
            else:
                j += 1
    return total


def union_length(intervals: Iterable[Tuple[str, int, int]]) -> int:
    """Base pairs covered by half-open intervals (sequence, start, end)."""
    by_sequence: Dict[str, List[Tuple[int, int]]] = {}
    for name, start, end in intervals:
        by_sequence.setdefault(name, []).append((start, end))
    total = 0
    for spans in by_sequence.values():
        spans.sort()
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start > current_end:
                total += current_end - current_start
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        total += current_end - current_start
    return total


# ---------------------------------------------------------------- running and writing

def container(sif: Path, *args: str, nice: bool = True) -> List[str]:
    command = ["apptainer", "exec", "--cleanenv", str(sif), *args]
    return (["nice", "-n", "10"] + command) if nice else command


def shell_text(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run(command: Sequence[str], go: bool, log: Optional[Path] = None) -> None:
    text = shell_text(command) + (f" > {shlex.quote(str(log))} 2>&1" if log else "")
    print(f"$ {text}", flush=True)
    if not go:
        return
    if log:
        with open(log, "a") as handle:
            result = subprocess.run(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT)
    else:
        result = subprocess.run(list(map(str, command)))
    if result.returncode != 0:
        stop(f"command failed with exit code {result.returncode}: {text}")


def copy_new(src: Path, dst: Path) -> None:
    """Copy a file to a place where nothing exists yet; never overwrites."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(src, "rb") as fin, open(dst, "xb") as fout:
        shutil.copyfileobj(fin, fout, 1024 * 1024)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with open(path, "x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def median(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2
