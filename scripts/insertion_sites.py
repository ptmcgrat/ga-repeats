#!/usr/bin/env python3
"""The known insertions anchored on published sequence (question 2, briefing C2.4).

Carleton et al. 2020 (Molecular Ecology 29:4956) give, in Supplementary Table S4, M. zebra UMD2a sequence at three
transposon indels (positions in their Table S3), and describe them (Supplementary Figs S3-S4):

  mitfa  LG5, intron 1: the UMD2a sequence on either side of the site where T. intermedius carries a 1,408 bp insertion
         (a 523 bp Rex1 match, LINE, plus 50-70 bp fragments of L1, Mariner and Copia). Without it intron 1 is 678 bp.
  rx1    LG23, about 2.5 kb upstream: the UMD2a sequence on either side of the 413 bp deletion that M. zebra and other
         short- and medium-palette species carry. The ancestral state keeps the 413 bp; long-palette species add an
         insertion of 268-421 bp inside it (384 bp of the 421 bp version match Rex1-5_AFC, LINE/Rex-Babar).
  tbx2a  LG10, about 13 kb upstream: the 1,332 bp of UMD2a that A. baenschi lacks. M. zebra (some individuals) carries a
         1,081 bp insertion there, 1,073 bp of it hAT-8_AFC (DNA/hAT).

For each assembly masked by the version, blastn (in the container) aligns these probes to the expected chromosome
(named from sequence_report.jsonl) and reports:
  mitfa, rx1  the two flanks' alignments (identity, coverage), the bases between them and the state that length implies,
              and the version's annotations there (one label per base: the highest-scoring);
  tbx2a       the locus where the probe's unique part aligns (parts aligning at several places on the chromosome are
              repeats, such as other hAT-8 copies), how much of the probe aligns there, and the annotations over it.
This checks that the probes were transcribed right (the M. zebra assemblies should align end to end), which species
carry each insertion, and whether our library finds and classifies it.

Usage (from the repository root):
  python scripts/insertion_sites.py --version RI_v5 --assemblies DT_004f MZ_MZ4f YH_YH2f MZ_NMK2024a ...          # dry run
  python scripts/insertion_sites.py --version RI_v5 --assemblies DT_004f MZ_MZ4f YH_YH2f MZ_NMK2024a ... --go

Writes, never overwriting, <scratch>/comparisons/insertion_sites/<version>_insertion_sites.tsv (one row per assembly and
locus) and _annotations.tsv (one row per annotation in a site). Work in <scratch>/<version>/insertion_sites/.
"""
import argparse
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402
import label_agreement as la  # noqa: E402

# Carleton et al. 2020, Supplementary Table S4 (M. zebra UMD2a), extracted from the supplement's PDF.
MITFA_LEFT = (
    "GAGTAGATGCGATTTTTACACTAAAAGAGTCAAAATATCTTTGCTGCGCTTGAAAACAGGACTTAAGGTTGAGGATTAGGGTTGAGTTCTAAATACTTCA"
    "TTATGTGTACCCAGTTCACAAGTCTTTGGTTAAATGAAGACAGGGAAGCTGTTTTACACACAGTCGTATTTAAATTGCTGCACGAGTTGGACAGGTGCTG"
    "AAATGTTATATAGTGCTCTTTTAGAAAAAAAGAGCTTGATGAAAAAAAATGTAAAAGAAATGTTGTTTGGACATTTCAACATCTTATCACCGTA"
)
MITFA_RIGHT = (
    "ATGTTGCAAAGCAAGCACGTTAGATTATTTGGTTTTTAAAAGGCTGTTTTCAGGGGAACTTTTAATGTAACAGAGTAGTATAGTGTGATTAATTTGGATT"
    "TGTAAGTTGCACTTTGAAGCTGGAAACATTTGAACATTTAAGTTCTTGCATTTTTAATTTGTGCTTTTACCACACCGGGACAGATGAATTATTATAATGA"
    "AAACATTCATTTTTA"
)
RX1_LEFT = (
    "CTGCATCCTCGGGCAAAATGCTGAACCCTACGTTTCCCCTGATGCATTTGAGTGAGTGTGTGAATGTGATTGTATGAGTGGGTGAATGAAAAGATGCTAT"
)
RX1_RIGHT = (
    "TTTACTGCGGGGCATGTGCAGAAAAAATACAGTCCTTTGTAGCACAGAAGTGTAACATAAAAACAAAAAGACTGTGACGAGGAATTTGAAATGATCACGG"
    "GACA"
)
TBX2A_REGION = (
    "AAAAAACTATTTATGCATTAATTACTTCATTTCACAGTGTCTAGACTGTTTGTTGAACTTGCCGCAACAGACAATTTGAAATTACAGCCAGCTTGGGAGT"
    "AAGTGATGTGACCCAAGATTGGAGCAGGATTTAGAGGTATGGTGAGTAGGGCTGGGTATCGTCACTGATTTCTAGAATCGATTCAATTCCGATTCACAAG"
    "GTCCCGAATCGATTCGATCCACGATTCGATTTAATTCGATTCGATTCGATTTAAATCTGGGAAATTTTGACAGTCAGAGATATAATTCAGATCAGGACAT"
    "TTACATATTTTTGTATCTATAAAAAGGAAGCTGACACTCTCAAGACTTTATCCAAGGTGTGAGCATCACAGCAGATGCCTTTGTGTCAAAGTAGCTGAAG"
    "ATAAAACACAGAAAAACATGTAGGTGATTTTCCTGGCCTGGGATTTTATAAAAATATTCTGCAGTACATCAAAAACGAAAGAAAACCATTAATTAATCAA"
    "CATTACCTCTGACGTTACAGCGGTTTTATTAGAGACATGGCTAAGCATTTTGCATTTTGAATAATTTTAAAAAGTTTCAATTTGTTCAGTATTGAACAGC"
    "AGAAATGAGGTTTTCTTTTCGGAAGAATGTAAAAGGGGAAAAAACAGCGGCCGACAGCGCTGTAAGCAACAGTAGACTTGTGTAACAAGCAAGCGAATAA"
    "TGAAGAAAGCGAAACTGTTCGAGAGAGAGAGAGAAAGGGAGAGGGAGAGAGAGAGAGCTGCGCATTGCAATGGAAGCAAAACAGTAAAAGAGTGAATTCA"
    "TGACGATGTTTATGTGAAGCGTTTGGATCTTCTTTTGCTGCTGATTCGGTCGATATTGTTTGGAGAGAGATCAAACTAACAGCTTTAGAATCAGCGCATA"
    "AAGGGCGTGAACACAAAGCGCGGACCCGGATCAGCGAGCTGTCGGCTTTCAGCCCCGACTGTGAGAAAGGCGACATCTAACTGATTCTGATCCGCGGGTT"
    "GCGCTGTGTGTTTAAGTCTATGTGCAGAATCCCTGTACCTGCTCTCATTTACTGTTTGGTGGTTCTTGAAATTTTGTGAGATGTCACCAGAGATTCTGGC"
    "ATTTCGGACAAAATAAATTTATATTAAAAAATCGATTCAGGATTTTAATGAATCGATTTTATGTTATCCAAGCCAGAATCGATTTTAATCGATGAATCGA"
    "TTATAAAAACCCACCCCTAATGGTGAGGTTAAGATTCTTAACTGTCTCTAAGTCCTTTGATTGTTTACCTGTAAAATCTTTTTCACCTCACACACATTAC"
    "GAACAGTTTGAAGGCAGATTGACCTGCCCAGG"
)

FLANKED = {   # locus: chromosome number, flanks, (lowest, highest, state) for the bases between the flanks
    "mitfa": ("5", MITFA_LEFT, MITFA_RIGHT, [(0, 100, "no insertion"), (1200, 1700, "insertion (1,408 bp published)")]),
    "rx1": ("23", RX1_LEFT, RX1_RIGHT, [(0, 100, "deletion (as M. zebra)"), (330, 500, "ancestral (413 bp)"),
                                        (600, 950, "insertion (413 + 268-421 bp)")]),
}
TBX2A_CHROMOSOME = "10"
MIN_IDENTITY, MIN_FLANK_COVERAGE, MAX_GAP = 90.0, 0.80, 20000


class Hit(NamedTuple):
    query: str
    subject: str
    identity: float
    qstart: int
    qend: int
    sstart: int
    send: int
    bitscore: float
    qlen: int

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def lo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def hi(self) -> int:
        return max(self.sstart, self.send)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="a version that has masked the assemblies, e.g. RI_v5")
    parser.add_argument("--assemblies", nargs="+", required=True, help="assemblies masked by it")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def parse_hits(lines) -> List[Hit]:
    """blastn -outfmt '6 qseqid sseqid pident qstart qend sstart send bitscore qlen'."""
    out = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9:
            out.append(Hit(f[0], f[1], float(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), float(f[7]), int(f[8])))
    return out


def site_between(left: Hit, right: Hit) -> Optional[Tuple[int, int]]:
    """0-based half-open stretch between a left and a right flank hit (same sequence and strand, in order), allowing
    for unaligned flank ends; None if they are not a pair."""
    if left.subject != right.subject or left.strand != right.strand:
        return None
    tail_left, head_right = left.qlen - left.qend, right.qstart - 1
    if left.strand == "+":
        if not left.hi < right.lo + 50 or right.lo - left.hi > MAX_GAP:
            return None
        start, end = left.hi + tail_left, right.lo - 1 - head_right
    else:
        if not right.hi < left.lo + 50 or left.lo - right.hi > MAX_GAP:
            return None
        start, end = right.hi + head_right, left.lo - 1 - tail_left
    return (start, max(start, end))


def best_site(hits: List[Hit], left_name: str, right_name: str):
    """The best pair of flank hits (most bitscore) and the stretch between them."""
    good = [h for h in hits if h.identity >= MIN_IDENTITY and (h.qend - h.qstart + 1) >= MIN_FLANK_COVERAGE * h.qlen]
    best = None
    for left in (h for h in good if h.query == left_name):
        for right in (h for h in good if h.query == right_name):
            site = site_between(left, right)
            if site is not None and (best is None or left.bitscore + right.bitscore > best[0]):
                best = (left.bitscore + right.bitscore, left, right, site)
    return best


def state_of(gap: int, states) -> str:
    for low, high, name in states:
        if low <= gap <= high:
            return name
    return f"other ({gap:,} bp)"


def tbx2a_locus(hits: List[Hit], window: int = 5000):
    """Cluster the region probe's hits by place; the locus is the cluster covering most of the probe's unique part
    (probe positions aligned at one place only). Returns (seq, lo, hi, probe bases aligned there, unique bases there,
    probe positions aligned at several places) or None."""
    hits = sorted((h for h in hits if h.identity >= 80 and h.qend - h.qstart + 1 >= 50), key=lambda h: (h.subject, h.lo))
    clusters: List[List[Hit]] = []
    for h in hits:
        if clusters and clusters[-1][-1].subject == h.subject and h.lo - max(x.hi for x in clusters[-1]) <= window:
            clusters[-1].append(h)
        else:
            clusters.append([h])
    if not clusters:
        return None
    covered = [set() for _ in clusters]
    for i, cl in enumerate(clusters):
        for h in cl:
            covered[i].update(range(h.qstart, h.qend + 1))
    count = Counter(p for s in covered for p in s)
    unique = {p for p, n in count.items() if n == 1}
    best = max(range(len(clusters)), key=lambda i: (len(covered[i] & unique), len(covered[i])))
    if not covered[best] & unique:
        return None
    cl = clusters[best]
    return (cl[0].subject, min(h.lo for h in cl) - 1, max(h.hi for h in cl), len(covered[best]),
            len(covered[best] & unique), len({p for p, n in count.items() if n > 1}))


def chromosome_sequences(fasta: Path, report: Path, wanted: Dict[str, Path]) -> Dict[str, str]:
    """Write each wanted chromosome (number -> output FASTA) out of the genome; returns number -> sequence name."""
    names = {}
    with open(report) as handle:
        for line in handle:
            if line.strip():
                d = json.loads(line)
                number = re.sub(r"^(LG|chr)0*", "", str(d.get("chrName") or ""), flags=re.I)
                for key in ("genbankAccession", "refseqAccession"):
                    if d.get(key) and number in wanted and d.get("role", "assembled-molecule") == "assembled-molecule":
                        names[d[key]] = number
    found, out = {}, None
    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if out:
                    out.close()
                    out = None
                seq = line[1:].split()[0]
                if seq in names:
                    found[names[seq]] = seq
                    out = open(wanted[names[seq]], "w")
                    out.write(f">{seq}\n")
            elif out:
                out.write(line)
    if out:
        out.close()
    return found


def annotations_in(rows: List[la.Row], seq: str, start: int, end: int):
    return sorted((max(a, start), min(b, end), f, l) for s, a, b, f, l in rows if s == seq and b > start and a < end)


def describe(hits) -> Tuple[int, str, int, int]:
    masked = sum(b - a for a, b, _, _ in hits)
    by_class = Counter()
    for a, b, _f, label in hits:
        by_class[la.group(label)] += b - a
    rex = sum(b - a for a, b, _f, label in hits if "Rex-Babar" in label)
    hat = sum(b - a for a, b, _f, label in hits if "hAT" in label)
    return masked, "; ".join(f"{k} {v}" for k, v in by_class.most_common()), rex, hat


def main(argv=None) -> None:
    args = parse_args(argv)
    plan = []
    for assembly in args.assemblies:
        genome = c.genome_path(assembly, c.panel_accession(assembly))
        report = genome.parent / "sequence_report.jsonl"
        bed = c.version_folder(args.version) / "masked" / assembly / f"{assembly}.bed"
        for f, what in ((genome, "genome"), (report, "sequence report"), (bed, f"{args.version}'s masking")):
            if not f.exists():
                c.stop(f"{assembly}: no {what} ({f})")
        plan.append((assembly, genome, report, bed))
    work = args.scratch / args.version / "insertion_sites"
    out_dir = args.scratch / "comparisons" / "insertion_sites"
    targets = [out_dir / f"{args.version}_insertion_sites.tsv", out_dir / f"{args.version}_insertion_sites_annotations.tsv"]
    for t in targets:
        if t.exists():
            c.stop(f"{t} exists; results are never overwritten")
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    print(f"{args.version}: mitfa (LG5), rx1 (LG23) and tbx2a (LG10) sites in " + ", ".join(p[0] for p in plan))
    c.run(c.container(args.sif, "blastn", "-version"), False)
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    work.mkdir(parents=True)
    c.run(c.container(args.sif, "blastn", "-version"), True, work / "blastn_version.txt")
    probes = work / "probes.fa"
    with open(probes, "x") as handle:
        for name, seq in (("mitfa_left", MITFA_LEFT), ("mitfa_right", MITFA_RIGHT), ("rx1_left", RX1_LEFT),
                          ("rx1_right", RX1_RIGHT), ("tbx2a_region", TBX2A_REGION)):
            handle.write(f">{name}\n{seq}\n")
    site_rows, annotation_rows = [], []
    chromosomes = {"5": "mitfa", "23": "rx1", TBX2A_CHROMOSOME: "tbx2a"}
    for assembly, genome, report, bed in plan:
        folder = work / assembly
        folder.mkdir()
        print(f"{assembly}: extracting chromosomes 5, 10 and 23 ...", flush=True)
        found = chromosome_sequences(genome, report, {n: folder / f"chr{n}.fa" for n in chromosomes})
        rows = la.one_label_per_base(la.read_labelled_bed(bed))
        for number, locus in chromosomes.items():
            if number not in found:
                site_rows.append([assembly, locus, "", "", "", "", "", "", "", "chromosome not found", "", "", "", ""])
                print(f"  {locus}: no chromosome {number} named in the sequence report")
                continue
            out = folder / f"{locus}_hits.tsv"
            c.run(c.container(args.sif, "blastn", "-query", str(probes), "-subject", str(folder / f"chr{number}.fa"),
                              "-evalue", "1e-10", "-outfmt", "6 qseqid sseqid pident qstart qend sstart send bitscore qlen",
                              "-out", str(out)), True, folder / f"{locus}_blastn.log")
            hits = parse_hits(open(out))
            if locus in FLANKED:
                _n, left_seq, right_seq, states = FLANKED[locus]
                best = best_site([h for h in hits if h.query.startswith(locus)], f"{locus}_left", f"{locus}_right")
                if best is None:
                    site_rows.append([assembly, locus, found[number], "", "", "", "", "", "", "flanks not found as a pair", "", "", "", ""])
                    print(f"  {locus}: the two flanks were not found as a pair on {found[number]}")
                    continue
                _score, left, right, (start, end) = best
                inside = annotations_in(rows, left.subject, start, end)
                masked, classes, rex, hat = describe(inside)
                state = state_of(end - start, states)
                flank_text = (f"left {left.identity:.1f}% over {left.qend - left.qstart + 1}/{left.qlen}, "
                              f"right {right.identity:.1f}% over {right.qend - right.qstart + 1}/{right.qlen}")
                site_rows.append([assembly, locus, left.subject, left.strand, start, end, end - start, flank_text, "", state,
                                  masked, classes, rex, hat])
                where = (f"{left.subject}:{start + 1:,}-{end:,}" if end > start else f"{left.subject}: flanks adjacent at {start:,}")
                print(f"  {locus}: {where} ({left.strand}) {end - start:,} bp between the flanks "
                      f"[{state}]; masked {masked:,} bp ({classes or 'nothing'}); Rex-Babar {rex:,} bp; flanks {flank_text}")
            else:
                loc = tbx2a_locus([h for h in hits if h.query == "tbx2a_region"])
                if loc is None:
                    site_rows.append([assembly, locus, found[number], "", "", "", "", "", "", "probe's unique part not found", "", "", "", ""])
                    print(f"  tbx2a: the probe's unique part was not found on {found[number]}")
                    continue
                seq, start, end, aligned, unique, repeated = loc
                inside = annotations_in(rows, seq, start, end)
                masked, classes, rex, hat = describe(inside)
                state = (f"{aligned:,} of {len(TBX2A_REGION):,} probe bases aligned here ({unique:,} unique to this place; "
                         f"{repeated:,} probe bases align at several places)")
                site_rows.append([assembly, locus, seq, "", start, end, end - start, "", aligned, state, masked, classes, rex, hat])
                print(f"  tbx2a: {seq}:{start + 1:,}-{end:,}: {state}; masked {masked:,} bp ({classes or 'nothing'}); hAT {hat:,} bp")
            for a, b, family, label in inside:
                annotation_rows.append([assembly, locus, left.subject if locus in FLANKED else seq, a, b, family, label,
                                        la.group(label), b - a])
    out_dir.mkdir(parents=True, exist_ok=True)
    c.write_tsv(targets[0], ["assembly", "locus", "seq", "strand", "start", "end", "length", "flanks", "probe_aligned",
                             "state", "masked_bp", "by_class", "rex_babar_bp", "hat_bp"], site_rows)
    c.write_tsv(targets[1], ["assembly", "locus", "seq", "start", "end", "family", "label", "class", "bp"], annotation_rows)
    print(f"tables: {targets[0]}\n        {targets[1]}")


if __name__ == "__main__":
    main()
