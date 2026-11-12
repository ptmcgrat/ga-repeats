#!/usr/bin/env python3
"""Merge several versions' libraries into one new version (spec section 4: de-duplication).

1. pool the sources' prefixed libraries (each family name already carries its build assembly);
2. cluster with cd-hit-est as compare_libraries.py does: local identity (-G 0), an identity
   threshold (-c) and a minimum alignment coverage of the shorter sequence (-aS);
3. keep cd-hit-est's representative of each cluster (its longest member);
4. give each cluster a label by the rule below and write it into the representative's header;
5. move the results into the new version's library/ folder:
     <version>-families.prefixed.fa   the merged library (what mask_genomes.py uses)
     cluster_table.tsv                every family: source, cluster, representative, identity and
                                      coverage to the representative, its own label, the cluster's
                                      label and the rule that decided it
     merge_record.json                sources and their library sha256, thresholds, commands, counts

Label rule (provisional, for the merge rehearsal; recorded in merge_record.json):
  a. no classified member: Unknown;
  b. all classified members of one class: that class (Unknown members inherit it);
  c. classified members disagree: the class given by the most source libraries;
  d. a tie: the representative's class if it is among the tied classes, else the class of the
     longest classified member among them.
  The full label (class/superfamily) is the representative's own when its class is chosen, else
  the most frequent full label among members of that class (ties: the longest member's).
  Clusters decided by c or d are flagged as conflicts in the cluster table.

Work happens in <scratch>/<version>/merge/; the version folder must not exist yet. Never overwrites.

Example:
  python scripts/merge_libraries.py --version RI_v4 --sources RI_v1 RI_v2 RI_v3            # dry run
  python scripts/merge_libraries.py --version RI_v4 --sources RI_v1 RI_v2 RI_v3 --go
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as c  # noqa: E402

RULE_TEXT = ("a: no classified member -> Unknown; b: classified members agree -> that class (Unknown members inherit); "
             "c: disagree -> class given by most source libraries; d: tie -> representative's class if tied, else the "
             "longest classified member's among the tied; full label: the representative's if its class is chosen, else "
             "the most frequent among members of that class (ties: the longest member's)")


class Member(NamedTuple):
    name: str        # full first word, e.g. MZ_MZ4f_rnd-1_family-23#LTR/Gypsy
    family: str      # before '#'
    label: str       # after '#'
    top_class: str
    length: int
    source: str      # source version, e.g. RI_v1


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="the new, merged version, e.g. RI_v4")
    parser.add_argument("--sources", nargs="+", required=True, help="versions to merge, e.g. RI_v1 RI_v2 RI_v3")
    parser.add_argument("--identity", type=float, default=0.80, help="cd-hit-est identity threshold (default 0.80)")
    parser.add_argument("--coverage", type=float, default=0.80,
                        help="minimum alignment coverage of the shorter sequence, cd-hit-est -aS (default 0.80)")
    c.add_common_args(parser)
    return parser.parse_args(argv)


def cluster_label(members: Sequence[Member], representative: Member) -> Tuple[str, str]:
    """(full label, rule letter) for one cluster."""
    classified = [m for m in members if m.top_class != "Unknown"]
    if not classified:
        return representative.label, "a"
    libraries: Dict[str, set] = defaultdict(set)
    for m in classified:
        libraries[m.top_class].add(m.source)
    if len(libraries) == 1:
        chosen, rule = next(iter(libraries)), "b"
    else:
        best = max(len(v) for v in libraries.values())
        tied = sorted(k for k, v in libraries.items() if len(v) == best)
        if len(tied) == 1:
            chosen, rule = tied[0], "c"
        elif representative.top_class in tied:
            chosen, rule = representative.top_class, "d"
        else:
            longest = max((m for m in classified if m.top_class in tied), key=lambda m: (m.length, m.name))
            chosen, rule = longest.top_class, "d"
    if representative.top_class == chosen:
        return representative.label, rule
    of_class = [m for m in classified if m.top_class == chosen]
    counts = Counter(m.label for m in of_class)
    top = max(counts.values())
    candidates = [m for m in of_class if counts[m.label] == top]
    return max(candidates, key=lambda m: (m.length, m.name)).label, rule


def coverage_of(member: Member, alignment: Optional[Tuple[int, int, int, int, str]]) -> Optional[float]:
    if alignment is None:
        return None
    qs, qe = alignment[0], alignment[1]
    return 100 * (abs(qe - qs) + 1) / member.length


def build_outputs(clusters: Dict[int, List[c.ClusterMember]], members: Dict[str, Member],
                  alignments: Dict[str, tuple]) -> Tuple[List[Tuple[str, str]], List[list], Counter]:
    """(header and representative name per cluster, cluster table rows, rule counts)."""
    headers, table, rules = [], [], Counter()
    for cluster in sorted(clusters):
        cms = clusters[cluster]
        reps = [m for m in cms if m.representative]
        if len(reps) != 1:
            c.stop(f"cluster {cluster} has {len(reps)} representatives")
        rep = members[reps[0].name]
        group = [members[m.name] for m in cms]
        label, rule = cluster_label(group, rep)
        rules[rule] += 1
        headers.append((f"{rep.family}#{label}", rep.name))
        for cm in cms:
            m = members[cm.name]
            cov = coverage_of(m, alignments.get(cm.name))
            table.append([cluster, m.name, m.source, m.length, "yes" if cm.representative else "no",
                          "" if cm.identity is None else f"{cm.identity:.2f}", "" if cov is None else f"{cov:.1f}",
                          m.label, label, rule, "yes" if rule in ("c", "d") else "no"])
    return headers, table, rules


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main(argv=None) -> None:
    args = parse_args(argv)
    if len(set(args.sources)) < 2 or len(set(args.sources)) != len(args.sources):
        c.stop("give two or more different source versions")
    if args.version in args.sources:
        c.stop("the merged version must be a new version")
    c.word_size(args.identity)  # raises for thresholds cd-hit-est cannot use

    members: Dict[str, Member] = {}
    sources = []
    for version in args.sources:
        build_id, _, prefixed = c.library_files(version)
        fams = c.read_families(prefixed)
        for f in fams:
            if f.name in members:
                c.stop(f"family name {f.name} occurs twice; were all libraries prefixed?")
            members[f.name] = Member(f.name, f.family, f.label, f.top_class, f.length, version)
        sources.append({"version": version, "build_assembly": build_id, "library": str(prefixed),
                        "library_sha256": None, "families": len(fams)})
        print(f"{version} ({build_id}): {len(fams)} families, {prefixed}")

    target = c.version_folder(args.version)
    work = args.scratch / args.version / "merge"
    print(f"new version {args.version}: {target}")
    print(f"work folder: {work}")
    print(f"cd-hit-est: identity {args.identity:.2f}, coverage of the shorter sequence {args.coverage:.2f}, local alignment")
    if target.exists():
        c.stop(f"{target} already exists; a version is never changed")
    if work.exists():
        c.stop(f"{work} exists from an earlier attempt; clear it with: rm -r {work}")
    if not args.go:
        print("dry run: nothing was run or written. Add --go to run.")
        return

    started = datetime.now()
    work.mkdir(parents=True)
    for s in sources:
        s["library_sha256"] = c.sha256(Path(s["library"]))
    pooled = work / "pooled.fa"
    with open(pooled, "x") as handle:
        for s in sources:
            text = Path(s["library"]).read_text()
            handle.write(text if text.endswith("\n") else text + "\n")
    prefix = work / f"cdhit_c{args.identity:.2f}"
    command = c.container(args.sif, c.CDHIT_EST, "-i", str(pooled), "-o", str(prefix) + ".fa",
                          "-c", f"{args.identity:.2f}", "-G", "0", "-aS", f"{args.coverage:.2f}",
                          "-n", str(c.word_size(args.identity)), "-g", "1", "-r", "1", "-d", "0", "-M", "0",
                          "-T", str(args.threads))
    c.run(command, True, work / "cdhit.log")

    clstr = Path(str(prefix) + ".fa.clstr")
    parsed = c.parse_clstr(clstr)
    if len(parsed) != len(members):
        c.stop(f"cd-hit-est clustered {len(parsed)} sequences, expected {len(members)}")
    clusters: Dict[int, List[c.ClusterMember]] = defaultdict(list)
    for m in parsed:
        clusters[m.cluster].append(m)
    headers, table, rules = build_outputs(clusters, members, c.clstr_alignments(clstr))

    representatives = read_fasta(Path(str(prefix) + ".fa"))
    merged = work / f"{args.version}-families.prefixed.fa"
    with open(merged, "x") as handle:
        for header, rep_name in headers:
            if rep_name not in representatives:
                c.stop(f"representative {rep_name} is missing from cd-hit-est's output")
            seq = representatives[rep_name]
            handle.write(f">{header}\n" + "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)))
    c.write_tsv(work / "cluster_table.tsv",
                ["cluster", "family", "source", "length", "representative", "identity_to_representative",
                 "coverage_of_family", "family_label", "cluster_label", "rule", "conflict"], table)

    cdhit_version = next((line.strip() for line in open(work / "cdhit.log") if "CD-HIT version" in line), None)
    record = {
        "version": args.version, "kind": "merged library", "sources": sources,
        "cd_hit_est": {"identity": args.identity, "coverage_of_shorter": args.coverage, "alignment": "local (-G 0)",
                       "version_line": cdhit_version, "command": c.shell_text(command)},
        "label_rule": RULE_TEXT,
        "counts": {"families_in": len(members), "clusters": len(clusters), "by_rule": dict(sorted(rules.items())),
                   "conflicts": rules["c"] + rules["d"],
                   "label_changed_on_representative": sum(1 for h, r in headers if h.split("#", 1)[1] != members[r].label)},
        "started": f"{started:%Y-%m-%d %H:%M:%S}", "finished": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
        "ga_repeats_version": (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip(),
        "image": str(args.sif),
    }
    with open(work / "merge_record.json", "x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")

    library = target / "library"
    for name in (merged.name, "cluster_table.tsv", "merge_record.json"):
        c.copy_new(work / name, library / name)
    counts = record["counts"]
    print(f"\n{len(members)} families -> {len(clusters)} clusters; rules: "
          + ", ".join(f"{k} {v}" for k, v in counts["by_rule"].items())
          + f"; conflicts {counts['conflicts']}; representative label changed {counts['label_changed_on_representative']}")
    print(f"written: {library}/ ({merged.name}, cluster_table.tsv, merge_record.json)")
    print(f"next: python scripts/mask_genomes.py --version {args.version} --assemblies DT_004f")


if __name__ == "__main__":
    main()
