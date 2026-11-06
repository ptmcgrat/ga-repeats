#!/usr/bin/env python3
"""Compare two pilot libraries (question 1 pilot, spec section 4).

Answers two questions for libraries A and B (for example RI_v1 and RI_v2):

  1. Which families are shared, and which are found in only one library?
     The two prefixed libraries are pooled and clustered with cd-hit-est at several
     identity thresholds (local identity, -G 0, with a minimum alignment coverage of
     the shorter sequence, -aS). A cluster holding both prefixes is a shared family.
  2. How similar are the matching families?
     RepeatMasker compares each library's consensus sequences against the other
     library; for every family it reports the best match (highest score), its
     divergence and how much of the family the match covers.

Results go to <scratch>/comparisons/<A>_vs_<B>/ (their permanent home is still open,
spec section 7). Nothing is written without --go, and nothing is ever overwritten.

Example:
  python scripts/compare_libraries.py --a RI_v1 --b RI_v2            # dry run
  python scripts/compare_libraries.py --a RI_v1 --b RI_v2 --go
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--a", required=True, help="first version, e.g. RI_v1")
    parser.add_argument("--b", required=True, help="second version, e.g. RI_v2")
    parser.add_argument("--thresholds", default="0.80,0.85,0.90,0.95",
                        help="cd-hit-est identity thresholds (default 0.80,0.85,0.90,0.95)")
    parser.add_argument("--coverage", type=float, default=0.80,
                        help="minimum alignment coverage of the shorter sequence (default 0.80)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def cluster_rows(members: Sequence[c.ClusterMember], families: Dict[str, c.Family],
                 source_a: str, source_b: str) -> List[list]:
    rows = []
    for m in members:
        fam = families[m.name]
        rows.append([m.cluster, m.name, c.source_of(m.name, [source_a, source_b]), fam.label,
                     fam.top_class, m.length, "yes" if m.representative else "no",
                     "" if m.identity is None else f"{m.identity:.2f}"])
    return rows


def shared_unique(members: Sequence[c.ClusterMember], families: Dict[str, c.Family],
                  source_a: str, source_b: str):
    """Per class (of the cluster representative): clusters and families shared or found in one library."""
    by_cluster: Dict[int, List[c.ClusterMember]] = defaultdict(list)
    for m in members:
        by_cluster[m.cluster].append(m)
    counts: Dict[str, Counter] = defaultdict(Counter)
    for cluster_members in by_cluster.values():
        rep = next(m for m in cluster_members if m.representative)
        cls = families[rep.name].top_class
        sources = {c.source_of(m.name, [source_a, source_b]) for m in cluster_members}
        kind = "shared" if len(sources) == 2 else ("only_a" if source_a in sources else "only_b")
        classes = {families[m.name].top_class for m in cluster_members} - {"Unknown"}
        for key in (cls, "ALL"):
            counts[key][f"clusters_{kind}"] += 1
            if len(classes) > 1:
                counts[key]["clusters_with_conflicting_classes"] += 1
            for m in cluster_members:
                src = "a" if c.source_of(m.name, [source_a, source_b]) == source_a else "b"
                counts[key][f"families_{src}_{'shared' if kind == 'shared' else 'only'}"] += 1
    columns = ["clusters_shared", "clusters_only_a", "clusters_only_b", "families_a_shared",
               "families_a_only", "families_b_shared", "families_b_only",
               "clusters_with_conflicting_classes"]
    order = sorted(k for k in counts if k != "ALL") + ["ALL"]
    return [[cls] + [counts[cls][col] for col in columns] for cls in order], columns


def best_matches(out_file: Path, query_names: Dict[str, str], families: Dict[str, c.Family],
                 coverage: float) -> List[list]:
    """One row per query family: its best hit in the other library (highest score)."""
    best: Dict[str, c.Hit] = {}
    if out_file.exists():
        for hit in c.parse_out(out_file):
            if hit.query not in best or hit.score > best[hit.query].score:
                best[hit.query] = hit
    rows = []
    for short, name in query_names.items():
        fam = families[name]
        hit = best.get(short)
        if hit is None:
            rows.append([name, fam.label, fam.top_class, fam.length, "", "", "", "", "", "no"])
            continue
        covered = (hit.end - hit.begin + 1) / fam.length
        rows.append([name, fam.label, fam.top_class, fam.length, hit.repeat, hit.label, hit.score,
                     f"{hit.divergence:.1f}", f"{covered:.3f}", "yes" if covered >= coverage else "no"])
    return rows


def similarity_by_class(rows: Sequence[list], direction: str) -> List[list]:
    grouped: Dict[str, List[list]] = defaultdict(list)
    for row in rows:
        grouped[row[2]].append(row)
        grouped["ALL"].append(row)
    out = []
    for cls in sorted(k for k in grouped if k != "ALL") + ["ALL"]:
        group = grouped[cls]
        matched = [r for r in group if r[4]]
        well = [r for r in matched if r[9] == "yes"]
        med = c.median([float(r[7]) for r in matched])
        out.append([direction, cls, len(group), len(matched), len(well),
                    "" if med is None else f"{med:.1f}"])
    return out


def main(argv=None) -> None:
    args = parse_args(argv)
    thresholds = [float(t) for t in args.thresholds.split(",")]
    for t in thresholds:
        c.word_size(t)  # raises for thresholds cd-hit-est cannot use

    id_a, _, lib_a = c.library_files(args.a)
    id_b, _, lib_b = c.library_files(args.b)
    fams_a = c.read_families(lib_a)
    fams_b = c.read_families(lib_b)
    families = {f.name: f for f in fams_a + fams_b}
    if len(families) != len(fams_a) + len(fams_b):
        c.stop("the two libraries share family names; were both prefixed?")
    for f in fams_a:
        c.source_of(f.name, [id_a])
    for f in fams_b:
        c.source_of(f.name, [id_b])

    out = args.scratch / "comparisons" / f"{args.a}_vs_{args.b}"
    pa = max(1, args.threads // c.RMBLAST_CORES_PER_JOB)
    print(f"A = {args.a} ({id_a}, {len(fams_a)} families): {lib_a}")
    print(f"B = {args.b} ({id_b}, {len(fams_b)} families): {lib_b}")
    print(f"results: {out}")
    print(f"cd-hit-est: thresholds {thresholds}, coverage of the shorter sequence {args.coverage}, "
          f"{args.threads} threads; RepeatMasker: -pa {pa} (about {pa * c.RMBLAST_CORES_PER_JOB} cores)")
    if out.exists():
        if not (out / "summary.txt").exists():
            c.stop(f"{out} holds an incomplete earlier run (no summary.txt); clear it with: rm -r {out}")
        c.stop(f"{out} already exists; results are never overwritten")
    if args.go:
        out.mkdir(parents=True)

    commands: List[str] = []

    def run(command, log=None):
        commands.append(c.shell_text(command))
        c.run(command, args.go, log)

    # families table and pooled library
    pooled = out / "pooled.fa"
    if args.go:
        c.write_tsv(out / "families.tsv", ["family", "source", "label", "class", "length"],
                    [[f.name, c.source_of(f.name, [id_a, id_b]), f.label, f.top_class, f.length]
                     for f in fams_a + fams_b])
        with open(pooled, "x") as handle:
            for lib in (lib_a, lib_b):
                text = Path(lib).read_text()
                handle.write(text if text.endswith("\n") else text + "\n")

    # question 1: clustering at each threshold
    summary_rows = []
    for t in thresholds:
        prefix = out / f"cdhit_c{t:.2f}"
        run(c.container(args.sif, c.CDHIT_EST, "-i", str(pooled), "-o", str(prefix) + ".fa",
                        "-c", f"{t:.2f}", "-G", "0", "-aS", f"{args.coverage:.2f}",
                        "-n", str(c.word_size(t)), "-g", "1", "-r", "1", "-d", "0", "-M", "0",
                        "-T", str(args.threads)),
            log=out / f"cdhit_c{t:.2f}.log")
        if not args.go:
            continue
        members = c.parse_clstr(Path(str(prefix) + ".fa.clstr"))
        if len(members) != len(families):
            c.stop(f"cd-hit-est clustered {len(members)} sequences, expected {len(families)}")
        c.write_tsv(out / f"clusters_c{t:.2f}.tsv",
                    ["cluster", "family", "source", "label", "class", "length", "representative",
                     "identity_to_representative"],
                    cluster_rows(members, families, id_a, id_b))
        rows, columns = shared_unique(members, families, id_a, id_b)
        summary_rows += [[f"{t:.2f}"] + r for r in rows]
    if args.go:
        c.write_tsv(out / "shared_unique_by_class.tsv", ["threshold", "class"] + columns, summary_rows)

    # question 2: best match of each family in the other library
    similarity_rows = []
    for query_fams, query_lib, lib, tag in ((fams_b, lib_b, lib_a, f"{args.b}_in_{args.a}"),
                                            (fams_a, lib_a, lib_b, f"{args.a}_in_{args.b}")):
        workdir = out / f"repeatmasker_{tag}"
        query = workdir / "query.fa"
        names = {f"q{i:06d}": f.name for i, f in enumerate(query_fams, start=1)}
        if args.go:
            workdir.mkdir()
            sequences = {}
            current = None
            with open(query_lib) as handle:
                for line in handle:
                    if line.startswith(">"):
                        current = line[1:].split(None, 1)[0]
                        sequences[current] = []
                    else:
                        sequences[current].append(line)
            with open(query, "x") as handle:
                for short, name in names.items():
                    handle.write(f">{short}\n{''.join(sequences[name])}")
            c.write_tsv(workdir / "query_names.tsv", ["query", "family"], names.items())
        run(c.container(args.sif, "RepeatMasker", "-lib", str(lib), "-engine", "rmblast",
                        "-pa", str(pa), "-nolow", "-no_is", "-dir", str(workdir), str(query)),
            log=workdir / "repeatmasker.log")
        if not args.go:
            continue
        rows = best_matches(workdir / "query.fa.out", names, families, args.coverage)
        c.write_tsv(out / f"best_match_{tag}.tsv",
                    ["family", "label", "class", "length", "best_match", "match_label", "score",
                     "divergence", "coverage", f"coverage_at_least_{args.coverage:.2f}"], rows)
        similarity_rows += similarity_by_class(rows, tag)
    if not args.go:
        print("dry run: nothing was written. Add --go to run.")
        return
    c.write_tsv(out / "similarity_by_class.tsv",
                ["direction", "class", "families", "with_a_match", f"match_covers_{args.coverage:.2f}",
                 "median_divergence_of_matches"], similarity_rows)

    # summary and record
    lines = [f"{args.a} ({id_a}) against {args.b} ({id_b}), {datetime.now():%Y-%m-%d %H:%M}", ""]
    for label, fams in ((args.a, fams_a), (args.b, fams_b)):
        unknown = sum(1 for f in fams if f.top_class == "Unknown")
        lines.append(f"{label}: {len(fams)} families, {unknown} Unknown ({100 * unknown / len(fams):.1f}%)")
    lines += ["", "Clusters at each threshold (ALL classes): shared / only A / only B"]
    for row in summary_rows:
        if row[1] == "ALL":
            lines.append(f"  {row[0]}: {row[2]} / {row[3]} / {row[4]}")
    lines += ["", "Best matches (ALL classes): families, with a match, match covering the threshold, "
              "median divergence"]
    for row in similarity_rows:
        if row[1] == "ALL":
            lines.append(f"  {row[0]}: {row[2]}, {row[3]}, {row[4]}, {row[5]}")
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    (out / "commands.txt").write_text("\n".join(commands) + "\n")
    (out / "meta.json").write_text(json.dumps({
        "a": args.a, "b": args.b, "id_a": id_a, "id_b": id_b,
        "library_a": str(lib_a), "library_b": str(lib_b),
        "library_a_sha256": c.sha256(lib_a), "library_b_sha256": c.sha256(lib_b),
        "thresholds": thresholds, "coverage": args.coverage, "image": str(args.sif),
        "date": f"{datetime.now():%Y-%m-%d %H:%M}"}, indent=2) + "\n")
    print((out / "summary.txt").read_text())


if __name__ == "__main__":
    main()
