"""Tests for the pilot scripts (repeats_009). Run from the repository root:

    python -m unittest discover -s tests

They need no container, no genome and no ga_common: they test the parsing and the
arithmetic on small fixtures. tests/data/RI_v1-rmod.log and builddatabase.log are the
real RI_v1 logs; example.out is the example in RepeatMasker 4.2.6's own help; cdhit_481_local.clstr holds
two real CD-HIT 4.8.1 member lines.
"""
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True  # keep __pycache__ out of the repository
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

import _common as c  # noqa: E402
import compare_libraries as cl  # noqa: E402
import summarize_masking as sm  # noqa: E402
import write_manifest as wm  # noqa: E402
import compare_masking as cm  # noqa: E402
import plot_pilot as pp  # noqa: E402
import masking_gain as mg  # noqa: E402
import label_agreement as la  # noqa: E402
import merge_libraries as ml  # noqa: E402
import evaluate_merge as em  # noqa: E402
import gene_overlap as go  # noqa: E402
import repbase_match as rbm  # noqa: E402
import gene_filter as gf  # noqa: E402
import compare_concatenation as cc  # noqa: E402
import known_insertions as ki  # noqa: E402
import own_genome_loss as ogl  # noqa: E402
import chromosome_masking as chm  # noqa: E402
import divergence_landscape as dl  # noqa: E402
import family_recovery as fr  # noqa: E402

DATA = HERE / "data"


class NamesAndClasses(unittest.TestCase):
    def test_split_name(self):
        self.assertEqual(c.split_name("MZ_MZ4f_rnd-1_family-23#LTR/Gypsy"),
                         ("MZ_MZ4f_rnd-1_family-23", "LTR/Gypsy"))
        self.assertEqual(c.split_name("plain"), ("plain", "Unknown"))
        self.assertEqual(c.top_class("LTR/Gypsy"), "LTR")
        self.assertEqual(c.top_class("Unknown"), "Unknown")
        self.assertEqual(c.top_class("PLE/Chlamys"), "PLE")

    def test_source_of_longest_prefix_wins(self):
        self.assertEqual(c.source_of("MZ_NMetZ_rnd-1_family-1#DNA", ["MZ_MZ4f", "MZ_NMetZ"]), "MZ_NMetZ")
        self.assertEqual(c.source_of("MZ_MZ4f_x#DNA", ["MZ", "MZ_MZ4f"]), "MZ_MZ4f")
        with self.assertRaises(ValueError):
            c.source_of("YH_YH2f_x#DNA", ["MZ_MZ4f"])

    def test_word_size(self):
        self.assertEqual([c.word_size(t) for t in (0.80, 0.85, 0.90, 0.95)], [5, 6, 8, 10])
        with self.assertRaises(ValueError):
            c.word_size(0.70)


class Fasta(unittest.TestCase):
    def test_lengths_and_families(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "lib.fa"
            path.write_text(">A_rnd-1_family-1#LTR/Gypsy [ Type=LTR ]\nACGT\nAC\n>A_rnd-1_family-2\nTTTT\n")
            self.assertEqual(c.fasta_lengths(path), {"A_rnd-1_family-1#LTR/Gypsy": 6, "A_rnd-1_family-2": 4})
            fams = c.read_families(path)
            self.assertEqual(fams[0].top_class, "LTR")
            self.assertEqual(fams[1].label, "Unknown")


class Clusters(unittest.TestCase):
    def test_parse_clstr(self):
        members = c.parse_clstr(DATA / "example.clstr")
        self.assertEqual(len(members), 7)
        self.assertTrue(members[0].representative)
        self.assertIsNone(members[0].identity)
        self.assertAlmostEqual(members[1].identity, 96.50)
        self.assertAlmostEqual(members[4].identity, 88.10)
        self.assertAlmostEqual(members[5].identity, 81.00)
        self.assertEqual(members[2].cluster, 1)

    def test_parse_clstr_local_alignment_lines(self):
        # The two member lines are real CD-HIT 4.8.1 output (-G 0) from the RI_v1 comparisons;
        # the representative lines around them are made up for the test.
        members = c.parse_clstr(DATA / "cdhit_481_local.clstr")
        self.assertEqual(len(members), 4)
        self.assertEqual(members[0].name, "MZ_MZ4f_rnd-1_family-373#DNA/CMC-EnSpm")
        self.assertEqual((members[0].length, members[0].representative), (167, False))
        self.assertAlmostEqual(members[0].identity, 96.27)
        self.assertTrue(members[1].representative)
        self.assertEqual(members[2].name, "YH_YH2f_rnd-2_family-534#DNA")
        self.assertAlmostEqual(members[2].identity, 87.14)

    def test_shared_unique(self):
        members = c.parse_clstr(DATA / "example.clstr")
        families = {}
        for m in members:
            fam, label = c.split_name(m.name)
            families[m.name] = c.Family(m.name, fam, label, c.top_class(label), m.length)
        rows, columns = cl.shared_unique(members, families, "MZ_MZ4f", "MZ_NMetZ")
        table = {r[0]: dict(zip(columns, r[1:])) for r in rows}
        # clusters 0 and 2 shared; cluster 1 only A; cluster 3 only B
        self.assertEqual(table["ALL"]["clusters_shared"], 2)
        self.assertEqual(table["ALL"]["clusters_only_a"], 1)
        self.assertEqual(table["ALL"]["clusters_only_b"], 1)
        self.assertEqual(table["ALL"]["families_a_shared"], 2)
        self.assertEqual(table["ALL"]["families_b_shared"], 3)
        self.assertEqual(table["ALL"]["families_a_only"], 1)
        self.assertEqual(table["ALL"]["families_b_only"], 1)
        # cluster 2 mixes DNA and LINE (Unknown ignored): one conflicting cluster
        self.assertEqual(table["ALL"]["clusters_with_conflicting_classes"], 1)
        # the class of a cluster is its representative's
        self.assertEqual(table["DNA"]["clusters_shared"], 1)
        self.assertEqual(table["DNA"]["clusters_only_b"], 1)


class RepeatMaskerOut(unittest.TestCase):
    def test_parse_out(self):
        hits = list(c.parse_out(DATA / "example.out"))
        self.assertEqual(len(hits), 7)
        first = hits[0]
        self.assertEqual((first.score, first.divergence, first.query, first.begin, first.end),
                         (207, 18.4, "Human", 75437, 75485))
        self.assertEqual((first.query_left, first.strand, first.repeat, first.label, first.repeat_id),
                         (924515, "C", "L2c", "LINE/L2", "101"))

    def test_bed_and_check(self):
        hits = list(c.parse_out(DATA / "example.out"))
        bed = [c.hit_to_bed(h) for h in hits]
        self.assertEqual(bed[0], ["Human", "75436", "75485", "L2c", "0", "-", "LINE/L2", "18.4", "207", "101"])
        self.assertEqual(bed[1][5], "+")
        self.assertEqual(c.check_bed(bed, {"Human": 1_000_000}), [])
        self.assertTrue(c.check_bed(bed, {"Human": 76000}))       # intervals past the end
        self.assertTrue(c.check_bed(bed, {"Chimp": 1_000_000}))   # unknown sequence name

    def test_best_matches(self):
        families = {"B_f1#LINE/L2": c.Family("B_f1#LINE/L2", "B_f1", "LINE/L2", "LINE", 100),
                    "B_f2#DNA": c.Family("B_f2#DNA", "B_f2", "DNA", "DNA", 50)}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "query.fa.out"
            out.write_text(
                "  SW  perc\nscore div.\n\n"
                "  500 10.0 0.0 0.0  q000001  1  90 (10) + A_x  LINE/L2  1 90 (0) 1\n"
                "  300 20.0 0.0 0.0  q000001  1  50 (50) + A_y  LINE/L2  1 50 (0) 2\n")
            rows = cl.best_matches(out, {"q000001": "B_f1#LINE/L2", "q000002": "B_f2#DNA"}, families, 0.8)
        self.assertEqual(rows[0][4], "A_x")          # highest score wins
        self.assertEqual(rows[0][8], "0.900")        # 90 of 100 bases
        self.assertEqual(rows[0][9], "yes")
        self.assertEqual(rows[1][4], "")             # no match
        summary = {r[1]: r for r in cl.similarity_by_class(rows, "B_in_A")}
        self.assertEqual(summary["ALL"][2:], [2, 1, 1, "10.0"])


class Intervals(unittest.TestCase):
    def test_union_length(self):
        self.assertEqual(c.union_length([("a", 0, 10), ("a", 5, 15), ("a", 20, 30), ("b", 0, 5)]), 30)
        self.assertEqual(c.union_length([("a", 0, 10), ("a", 10, 20)]), 20)
        self.assertEqual(c.union_length([("a", 0, 100), ("a", 10, 20)]), 100)

    def test_masking_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "clusters.tsv"
            table.write_text("cluster\tfamily\tsource\n"
                             "0\tA_f1#LTR\tA\n0\tB_f1#LTR\tB\n1\tA_f2#DNA\tA\n2\tB_f2#DNA\tB\n")
            unique = sm.only_in_one(table)
        self.assertEqual(unique["A"], {"A_f2"})
        self.assertEqual(unique["B"], {"B_f2"})
        bed = [("chr1", 0, 100, "A_f1", "LTR"), ("chr1", 50, 120, "A_f2", "DNA"),
               ("chr1", 500, 600, "A_f2", "DNA")]
        stats = sm.summarize(bed, unique["A"])
        self.assertEqual(stats["ALL"], (220, 170))
        self.assertEqual(stats["LTR"], (100, 0))
        self.assertEqual(stats["DNA"], (170, 170))


class Logs(unittest.TestCase):
    def test_builddatabase(self):
        text = (DATA / "builddatabase.log").read_text()
        self.assertEqual(wm.parse_builddatabase(text), (61, 972958384))

    def test_rmod_log(self):
        info = wm.parse_rmod_log((DATA / "RI_v1-rmod.log").read_text())
        self.assertEqual(info["repeatmodeler"], "2.0.9")
        self.assertEqual(info["threads"], "16")
        self.assertEqual(info["random_seed"], "1793333684")
        self.assertEqual((info["database_sequences"], info["database_bases"]), ("61", "972958384"))
        self.assertEqual(info["rounds"], 5)
        self.assertEqual(info["final_coverage"], "408053349 of 972958384 bp (41.94%)")
        self.assertEqual(info["round1_sample_bp"], "40808826")
        self.assertEqual((info["repeatscout_recon_families"], info["ltr_pipeline_families"],
                          info["ltr_redundant_removed"]), ("1787", "848", "239"))
        self.assertEqual(info["program_time"], "36:14:33")
        self.assertEqual(info["dfam"], "4.0")
        self.assertEqual(info["dfam_reference_sizes"], ["18011", "30659"])
        self.assertTrue(str(info["ltr_structural_analysis"]).startswith("Enabled ( GenomeTools 1.6.4"))
        self.assertIn("RepeatMasker 4.2.6", str(info["dependencies"]))

    def test_run_dates(self):
        info = wm.parse_rmod_log((DATA / "RI_v1-rmod.log").read_text())
        start = wm.run_start(info["rm_directory"])
        self.assertEqual(start, datetime(2026, 10, 30, 0, 14, 44))
        self.assertEqual(wm.add_duration(start, "36:14:33"), datetime(2026, 10, 31, 12, 29, 17))

    def test_yaml(self):
        text = "\n".join(wm.to_yaml({"version": "RI_v1", "n": 3, "ok": True, "none": None,
                                     "list": [{"a": "x", "b": ["p", "q"]}], "s": 'say "hi"'}))
        self.assertIn('version: "RI_v1"', text)
        self.assertIn("  - a: \"x\"", text)
        try:
            import yaml  # optional: only checks the text is valid YAML when PyYAML is present
        except ImportError:
            return
        self.assertEqual(yaml.safe_load(text)["list"][0]["b"], ["p", "q"])
        self.assertEqual(yaml.safe_load(text)["s"], 'say "hi"')


class MaskingOverlap(unittest.TestCase):
    """A masks s:0-10 (LTR) and s:20-30 (DNA); B masks s:5-25 (DNA); the genome is 100 bp."""

    def setUp(self):
        a = [("s", 0, 10, "A_f1", "LTR"), ("s", 20, 30, "A_f2", "DNA")]
        b = [("s", 5, 25, "B_f1", "DNA")]
        self.rows = cm.pair_rows("G", "A", "B", {"A": cm.masked_sets(a), "B": cm.masked_sets(b)}, 100)

    def test_subtract_length(self):
        m = c.merge_by_sequence
        self.assertEqual(c.subtract_length(m([("s", 0, 10)]), m([("s", 2, 4), ("s", 6, 8)])), 6)
        self.assertEqual(c.subtract_length(m([("s", 0, 10), ("s", 20, 30)]), m([("s", 5, 25)])), 10)
        self.assertEqual(c.subtract_length(m([("s", 0, 10)]), {}), 10)
        self.assertEqual(c.subtract_length(m([("s", 0, 10)]), m([("s", 0, 3), ("s", 2, 6), ("s", 9, 40)])), 3)

    def test_pair_rows(self):
        table = {(r[2], r[3]): (r[4], r[5]) for r in self.rows}
        self.assertEqual(table[("A", "ALL")], (20, 10))   # A masks 20; 10 of them B does not
        self.assertEqual(table[("B", "ALL")], (20, 10))
        self.assertEqual(table[("A", "LTR")], (10, 5))
        self.assertEqual(table[("A", "DNA")], (10, 5))
        self.assertEqual(table[("B", "DNA")], (20, 10))

    def test_overlap_bars(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "G.tsv"
            c.write_tsv(table, ["assembly", "pair", "library", "class", "masked_bp", "masked_bp_not_masked_by_other",
                                "genome_bp", "percent_masked", "percent_not_masked_by_other"], self.rows)
            bars = pp.overlap_bars(Path(tmp))
        both, only_a, only_b, a, b = bars["G: A vs B"]
        self.assertEqual((round(both, 3), round(only_a, 3), round(only_b, 3), a, b), (10.0, 10.0, 10.0, "A", "B"))


class Plots(unittest.TestCase):
    def test_class_composition_counts_each_masking_once(self):
        rows = [{"library_version": "RI_v1", "library_built_from": "X", "assembly": "G", "class": cls,
                 "percent_masked": pct} for cls, pct in (("DNA", "10.0"), ("Satellite", "1.0"), ("ALL", "12.0"))]
        comp = {"masking": rows}
        out = pp.class_composition([comp, comp])
        group = out["RI_v1 (X) on G"]
        self.assertEqual((group["DNA"], group["Other"], group["_all"]), (10.0, 1.0, 12.0))

    def test_figures_are_drawn(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not installed here")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "RI_v1_vs_RI_v2"
            folder.mkdir()
            (folder / "meta.json").write_text('{"a": "RI_v1", "b": "RI_v2", "id_a": "X", "id_b": "Y"}')
            c.write_tsv(folder / "shared_unique_by_class.tsv",
                        ["threshold", "class", "clusters_shared", "clusters_only_a", "clusters_only_b"],
                        [["0.80", cls, 5, 3, 2] for cls in ("LTR", "DNA", "LINE", "SINE", "Unknown", "ALL")])
            c.write_tsv(folder / "similarity_by_class.tsv",
                        ["direction", "class", "families", "with_a_match", "match_covers_0.80",
                         "median_divergence_of_matches"],
                        [["RI_v2_in_RI_v1", "ALL", 10, 9, 6, "2.0"], ["RI_v1_in_RI_v2", "ALL", 10, 8, 5, "2.5"]])
            overlap = Path(tmp) / "overlap"
            overlap.mkdir()
            c.write_tsv(overlap / "G.tsv", ["assembly", "pair", "library", "class", "masked_bp",
                                            "masked_bp_not_masked_by_other", "genome_bp", "percent_masked",
                                            "percent_not_masked_by_other"], MaskingOverlapRows)
            out = Path(tmp) / "plots"
            pp.main(["--comparisons", str(folder), "--overlap-dir", str(overlap), "--out", str(out), "--go"])
            made = sorted(p.name for p in out.iterdir())
        self.assertIn("fig1_shared_clusters_by_threshold.png", made)
        self.assertIn("fig5_masking_overlap.pdf", made)
        self.assertNotIn("fig4_class_composition.png", made)   # no masking summary in this folder


class LabelAgreement(unittest.TestCase):
    """A: s 0-10 LTR (f1), s 20-30 DNA (f2). B: s 5-15 DNA (g1), s 20-30 DNA (g2), s 40-50 Unknown (g3)."""

    def setUp(self):
        self.a = [("s", 0, 10, "f1", "LTR/Gypsy"), ("s", 20, 30, "f2", "DNA/hAT")]
        self.b = [("s", 5, 15, "g1", "DNA/TcMar"), ("s", 20, 30, "g2", "DNA/hAT"), ("s", 40, 50, "g3", "Unknown")]
        self.m = la.class_matrix(self.a, self.b)

    def test_class_matrix_cells(self):
        self.assertEqual(self.m[("LTR", "DNA")], 5)        # 5-10
        self.assertEqual(self.m[("LTR", "unmasked")], 5)   # 0-5
        self.assertEqual(self.m[("DNA", "DNA")], 10)       # 20-30
        self.assertEqual(self.m[("DNA", "unmasked")], 0)
        self.assertEqual(self.m[("unmasked", "DNA")], 5)   # 10-15
        self.assertEqual(self.m[("unmasked", "Unknown")], 10)

    def test_agreement_counts_only_bases_masked_by_both(self):
        self.assertEqual(la.agreement(self.m), (15, 10))

    def test_discordant_family_pairs(self):
        pairs = la.discordant_families(self.a, self.b)
        self.assertEqual(dict(pairs), {("f1", "LTR/Gypsy", "g1", "DNA/TcMar"): 5})

    def test_small_classes_are_grouped_as_other(self):
        self.assertEqual([la.group(x) for x in ("Simple_repeat", "SINE?", "RC/Helitron", "LINE/L2", "Unknown")],
                         ["Other", "Other", "Other", "LINE", "Unknown"])


def _m(name, label, length, source):
    family = name
    return ml.Member(f"{family}#{label}", family, label, label.split("/")[0], length, source)


class MergeLabelRule(unittest.TestCase):
    def test_a_no_classified_member_stays_unknown(self):
        rep = _m("A1", "Unknown", 900, "RI_v1")
        self.assertEqual(ml.cluster_label([rep, _m("B1", "Unknown", 500, "RI_v3")], rep), ("Unknown", "a"))

    def test_b_unknown_members_inherit_the_agreed_class(self):
        rep = _m("A1", "Unknown", 900, "RI_v1")
        group = [rep, _m("B1", "DNA/hAT-Ac", 500, "RI_v2"), _m("C1", "DNA/hAT-Ac", 400, "RI_v3"), _m("C2", "DNA/TcMar", 300, "RI_v3")]
        self.assertEqual(ml.cluster_label(group, rep), ("DNA/hAT-Ac", "b"))

    def test_b_keeps_the_representatives_full_label_when_its_class_is_chosen(self):
        rep = _m("A1", "LTR/Pao", 900, "RI_v1")
        group = [rep, _m("B1", "LTR/Gypsy", 500, "RI_v2"), _m("C1", "LTR/Gypsy", 400, "RI_v3")]
        self.assertEqual(ml.cluster_label(group, rep), ("LTR/Pao", "b"))

    def test_c_the_class_given_by_most_libraries_wins(self):
        rep = _m("A1", "DNA/hAT-Ac", 900, "RI_v1")
        group = [rep, _m("B1", "LTR/Gypsy", 500, "RI_v2"), _m("C1", "LTR/Gypsy", 400, "RI_v3")]
        self.assertEqual(ml.cluster_label(group, rep), ("LTR/Gypsy", "c"))

    def test_libraries_are_counted_not_families(self):
        rep = _m("A1", "DNA/hAT-Ac", 900, "RI_v1")
        group = [rep, _m("C1", "LTR/Gypsy", 500, "RI_v3"), _m("C2", "LTR/Gypsy", 400, "RI_v3")]
        self.assertEqual(ml.cluster_label(group, rep), ("DNA/hAT-Ac", "d"))   # one library each: a tie

    def test_d_tie_without_the_representative_goes_to_the_longest_member(self):
        rep = _m("A1", "Unknown", 900, "RI_v1")
        group = [rep, _m("B1", "LTR/Gypsy", 500, "RI_v2"), _m("C1", "DNA/hAT-Ac", 800, "RI_v3")]
        self.assertEqual(ml.cluster_label(group, rep), ("DNA/hAT-Ac", "d"))


class MergeOutputs(unittest.TestCase):
    def test_headers_table_and_coverage(self):
        members = {m.name: m for m in (_m("A1", "Unknown", 900, "RI_v1"), _m("B1", "LTR/Gypsy", 450, "RI_v2"),
                                       _m("C1", "DNA/hAT-Ac", 300, "RI_v3"))}
        clusters = {0: [c.ClusterMember(0, "A1#Unknown", 900, True, None), c.ClusterMember(0, "B1#LTR/Gypsy", 450, False, 91.5)],
                    1: [c.ClusterMember(1, "C1#DNA/hAT-Ac", 300, True, None)]}
        alignments = {"B1#LTR/Gypsy": (450, 1, 10, 459, "-")}
        headers, table, rules = ml.build_outputs(clusters, members, alignments)
        self.assertEqual(headers, [("A1#LTR/Gypsy", "A1#Unknown"), ("C1#DNA/hAT-Ac", "C1#DNA/hAT-Ac")])
        self.assertEqual(rules, {"b": 2})
        b1 = next(r for r in table if r[1] == "B1#LTR/Gypsy")
        self.assertEqual(b1[5:], ["91.50", "100.0", "50.0", "LTR/Gypsy", "LTR/Gypsy", "b", "no"])   # 450 of 900 bp of A1


class MergeCommand(unittest.TestCase):
    def _command(self, *extra):
        import io, contextlib, tempfile
        captured = {}
        saved = (c.run, c.library_files, c.read_families, c.version_folder)
        tmp = Path(tempfile.mkdtemp())
        def fake_run(command, go, log=None):
            captured["command"] = list(command)
            raise SystemExit(0)
        c.run = fake_run
        c.library_files = lambda v: (v, tmp / f"{v}.fa", tmp / f"{v}.fa")
        c.read_families = lambda path: [c.Family(f"{path.stem}_f1#DNA", f"{path.stem}_f1", "DNA", "DNA", 100)]
        c.version_folder = lambda v: tmp / "versions" / v
        for v in ("RI_v1", "RI_v2"):
            (tmp / f"{v}.fa").write_text(f">{v}_f1#DNA\nACGT\n")
        try:
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                ml.main(["--version", "RI_v9", "--sources", "RI_v1", "RI_v2", "--scratch", str(tmp / "scratch"), "--go", *extra])
        finally:
            c.run, c.library_files, c.read_families, c.version_folder = saved
        return captured["command"]

    def test_coverage_of_the_longer_sequence_is_passed_as_aL(self):
        command = self._command("--coverage-longer", "0.8")
        self.assertIn("-aL", command)
        self.assertEqual(command[command.index("-aL") + 1], "0.80")

    def test_without_the_option_no_aL_is_passed(self):
        self.assertNotIn("-aL", self._command())


class CopyNumberFilter(unittest.TestCase):
    """Consensus fam1 is 1,000 bp. RepeatMasker .out lines: + strand 'begin end (left)', C strand '(left) end begin'."""
    OUT = [
        "   SW   perc perc perc  query      position in query           matching  repeat     position in repeat",
        "  900  5.0  0.0  0.0  s1  1001  1900  (100)  +  fam1  LTR/Gypsy  1  900  (100)  1",        # 90% of it: full length
        "  500 10.0  0.0  0.0  s1  5001  5500  (100)  C  fam1  LTR/Gypsy  (0)  1000  501  2",      # fragments of one copy ...
        "  400 10.0  0.0  0.0  s1  5601  5950  (100)  C  fam1  LTR/Gypsy  (500)  500  151  2",     # ... 501-1000 + 151-500 = 85%
        "  300 25.0  0.0  0.0  s1  9001  9900  (100)  +  fam1  LTR/Gypsy  1  900  (100)  3",        # 75% identity: not full length
        "  100  2.0  0.0  0.0  s1  12001 12500 (100)  +  fam1  LTR/Gypsy  1  500  (500)  4",        # covers 50%: not full length
        "  100  2.0  0.0  0.0  s1  13001 13060 (100)  +  fam2  DNA        1  60   (0)    5",        # fam2 is 60 bp: under 80 bp
    ]

    def test_parsing_both_strands(self):
        hits = list(ml.consensus_hits(self.OUT))
        self.assertEqual(hits[0], ("fam1", "1", 900, 5.0, 1, 900))
        self.assertEqual(hits[1], ("fam1", "2", 500, 10.0, 501, 1000))

    def test_full_length_copies_join_fragments_and_apply_80_80_80(self):
        counts = ml.full_length_copies(ml.consensus_hits(self.OUT), {"fam1": 1000, "fam2": 60})
        self.assertEqual(counts["fam1"], 2)   # copies 1 and 2
        self.assertEqual(counts["fam2"], 0)


class ConcatenationTest(unittest.TestCase):
    def _tmp(self):
        import tempfile
        return Path(tempfile.mkdtemp())

    def test_fragmentation_counts_interspersed_annotations_and_elements(self):
        bed = self._tmp() / "g.bed"
        bed.write_text("#chrom\tstart\tend\tfamily\tscore\tstrand\tclass_family\tdivergence\tsw_score\trm_id\n"
                       "chr1\t0\t100\tfamA\t0\t+\tLTR/Gypsy\t5\t300\t1\n"
                       "chr1\t150\t250\tfamA\t0\t+\tLTR/Gypsy\t5\t280\t1\n"      # second fragment of element 1
                       "chr1\t300\t330\t(CA)n\t0\t+\tSimple_repeat\t2\t40\t2\n"  # left out
                       "chr2\t0\t50\tfamB\t0\t+\tDNA/hAT\t9\t200\t1\n")         # element 1 of another sequence
        self.assertEqual(cc.fragmentation(bed), {"annotations": 3, "elements": 2, "median_length": 100})

    def test_cluster_members_keeps_clusters_of_two_or_more(self):
        table = self._tmp() / "cluster_table.tsv"
        table.write_text("cluster\tfamily\tsource\n0\tA_rnd-1_family-1#LTR/Gypsy\tRI_v1\n"
                         "0\tB_rnd-1_family-7#LTR/Gypsy\tRI_v2\n1\tC_rnd-2_family-4#DNA\tRI_v3\n")
        self.assertEqual(cc.cluster_members(table), {"A_rnd-1_family-1": "0", "B_rnd-1_family-7": "0"})

    def test_family_split(self):
        member_of = {"f1": "0", "f2": "0", "f3": "1", "f4": "1", "f5": "1"}
        rows = [("s", 0, 90, "f1", "LTR"), ("s", 100, 110, "f2", "LTR"), ("s", 200, 240, "f3", "DNA"),
                ("s", 300, 335, "f4", "DNA"), ("s", 400, 425, "f5", "DNA"), ("s", 500, 600, "other", "LINE")]
        split = cc.family_split(rows, member_of)
        self.assertEqual(dict(split["0"]), {"f1": 90, "f2": 10})
        s = cc.split_summary(split)
        self.assertEqual((s["clusters"], s["bp"], s["clusters_no_member_half"], s["bp_no_member_half"], s["clusters_two_members_10pct"]),
                         (2, 200, 1, 100, 2))
        self.assertAlmostEqual(s["top_member_share"], (90 + 40) / 200)

    def test_run_minutes_from_the_run_record(self):
        scratch = self._tmp()
        folder = scratch / "RI_v8" / "masking" / "DT_004f"
        folder.mkdir(parents=True)
        (folder / "run_record.json").write_text('{"started": "2026-11-17 10:05", "finished": "2026-11-17 14:35"}')
        self.assertEqual(cc.run_minutes(scratch, "RI_v8", "DT_004f"), 270)
        self.assertIsNone(cc.run_minutes(scratch, "RI_v5", "DT_004f"))


class KnownInsertionsTest(unittest.TestCase):
    GFF = [
        "##gff-version 3",
        "chr1\tRefSeq\tgene\t101\t600\t.\t+\t.\tID=gene-rx1;Name=rx1;gene=rx1",
        "chr1\tRefSeq\tmRNA\t101\t600\t.\t+\t.\tID=rna-1;Parent=gene-rx1",
        "chr1\tRefSeq\texon\t101\t200\t.\t+\t.\tParent=rna-1",
        "chr1\tRefSeq\texon\t301\t400\t.\t+\t.\tParent=rna-1",
        "chr1\tRefSeq\tCDS\t151\t200\t.\t+\t0\tParent=rna-1",
        "chr1\tRefSeq\tCDS\t301\t350\t.\t+\t1\tParent=rna-1",
        "chr2\tRefSeq\tgene\t500\t1100\t.\t-\t.\tID=gene-mitfa;Name=mitfa;gene=mitfa",
        "chr2\tRefSeq\tmRNA\t500\t1100\t.\t-\t.\tID=rna-2;Parent=gene-mitfa",
        "chr2\tRefSeq\texon\t1000\t1100\t.\t-\t.\tParent=rna-2",
        "chr2\tRefSeq\texon\t800\t900\t.\t-\t.\tParent=rna-2",
        "chr2\tRefSeq\texon\t500\t600\t.\t-\t.\tParent=rna-2",
        "chr2\tRefSeq\tCDS\t800\t1050\t.\t-\t0\tParent=rna-2",
        "chr2\tRefSeq\tgene\t5\t50\t.\t+\t.\tID=gene-other;Name=other",
    ]

    def test_gene_models_and_regions(self):
        genes = {g["name"]: g for g in ki.read_gene_models(self.GFF, ki.GENES)}
        self.assertEqual(set(genes), {"rx1", "mitfa"})
        self.assertEqual(ki.intron1(genes["rx1"]["transcripts"]["rna-1"]["exons"], "+"), (200, 300))
        self.assertEqual(ki.intron1(genes["mitfa"]["transcripts"]["rna-2"]["exons"], "-"), (900, 999))   # between the two highest exons
        self.assertEqual(ki.regions(genes["rx1"], 100), [(50, 150)])                                   # 100 bp 5' of the CDS start (151)
        self.assertEqual(ki.upstream(genes["mitfa"]["transcripts"]["rna-2"]["cds"], "-", 100), (1050, 1150))
        self.assertEqual(ki.regions(genes["mitfa"], 100), [(900, 999)])

    def test_annotations_are_clipped_to_the_region(self):
        rows = [("chr2", 850, 950, "fA", "LINE/Rex-Babar"), ("chr2", 980, 1200, "fB", "DNA/hAT"), ("chr1", 900, 990, "fC", "LTR")]
        self.assertEqual(ki.annotations_in(rows, "chr2", 900, 999),
                         [(900, 950, "fA", "LINE/Rex-Babar"), (980, 999, "fB", "DNA/hAT")])


class OwnGenomeLossTest(unittest.TestCase):
    def test_lost_bases_per_family(self):
        own = [("s", 0, 100, "fA", "LTR/Gypsy"), ("s", 200, 300, "fB", "DNA/hAT"), ("t", 0, 10, "fC", "Unknown")]
        merged = [("s", 50, 250, "x", "LTR/Gypsy")]
        lost = ogl.lost_by_family(own, merged)
        self.assertEqual(dict(lost), {("fA", "LTR/Gypsy"): 50, ("fB", "DNA/hAT"): 50, ("fC", "Unknown"): 10})

    def test_cluster_info(self):
        import tempfile
        table = Path(tempfile.mkdtemp()) / "cluster_table.tsv"
        table.write_text("cluster\tfamily\trepresentative\tcluster_label\n"
                         "0\tA_f1#LTR/Gypsy\tyes\tLTR/Gypsy\n0\tB_f7#LTR/Gypsy\tno\tLTR/Gypsy\n1\tC_f4#DNA\tyes\tDNA\n")
        info = ogl.cluster_info(table)
        self.assertEqual(info["B_f7"], {"cluster": "0", "role": "member", "size": 2, "representative": "A_f1", "cluster_label": "LTR/Gypsy"})
        self.assertEqual(info["C_f4"]["role"], "representative")


class ChromosomeMaskingTest(unittest.TestCase):
    REPORT = [
        '{"chrName": "3", "role": "assembled-molecule", "genbankAccession": "CM3.1", "refseqAccession": "NC_3.1", "length": 1000}',
        '{"chrName": "13", "role": "assembled-molecule", "genbankAccession": "CM13.1", "refseqAccession": "NC_13.1", "length": 500}',
        '{"chrName": "Un", "role": "unplaced-scaffold", "genbankAccession": "JA1.1", "length": 100}',
    ]

    def test_names_and_lg3(self):
        genbank, refseq = chm.chromosome_names(self.REPORT)
        self.assertEqual(genbank, {"CM3.1": ("3", 1000), "CM13.1": ("13", 500), "JA1.1": ("unplaced", 100)})
        self.assertEqual(refseq, {"NC_3.1": ("3", 1000), "NC_13.1": ("13", 500)})
        self.assertTrue(all(chm.LG3.match(x) for x in ("3", "LG3", "chr03")))
        self.assertFalse(any(chm.LG3.match(x) for x in ("13", "LG13", "23")))

    def test_per_chromosome(self):
        genbank, _ = chm.chromosome_names(self.REPORT)
        raw = [("CM3.1", 0, 100, "f", "LTR/Gypsy", 300), ("CM3.1", 50, 150, "g", "DNA/hAT", 200), ("JA1.1", 0, 10, "h", "Unknown", 50)]
        stats = chm.per_chromosome(genbank, raw, [("CM3.1", 0, 100, "f", "LTR/Gypsy"), ("CM3.1", 100, 150, "g", "DNA/hAT"),
                                                 ("JA1.1", 0, 10, "h", "Unknown")])
        self.assertEqual((stats["3"]["length"], stats["3"]["masked"], stats["3"]["LTR"], stats["3"]["DNA"]), (1000, 150, 100, 50))
        self.assertEqual((stats["unplaced"]["length"], stats["unplaced"]["masked"]), (100, 10))
        self.assertEqual(stats["13"]["masked"], 0)


class DivergenceLandscapeTest(unittest.TestCase):
    def test_bins_and_median(self):
        lines = ["#chrom\tstart\tend\tfamily\tscore\tstrand\tclass_family\tdivergence\tsw_score\trm_id",
                 "s\t0\t100\tf\t0\t+\tLTR/Gypsy\t2.5\t300\t1",
                 "s\t100\t130\tg\t0\t+\tLTR/Gypsy\t12.0\t200\t2",
                 "s\t200\t210\t(CA)n\t0\t+\tSimple_repeat\t1.0\t20\t3",
                 "s\t300\t400\th\t0\t+\tDNA/hAT\t70.0\t100\t4"]
        land = dl.landscape(lines, 1.0, 50.0)
        self.assertEqual(dict(land["LTR"]), {2.0: 100, 12.0: 30})
        self.assertEqual(dict(land["DNA"]), {50.0: 100})     # everything from 50% on goes to the last bin
        self.assertNotIn("Other", land)                       # simple repeats left out
        self.assertEqual(dl.weighted_median(land["LTR"]), 2.0)


class FamilyRecoveryTest(unittest.TestCase):
    def test_parsing_both_strands(self):
        lines = ["  500  2.0  0.0  0.0  q000001  1  980  (0)  +  CICH1  DNA/hAT  1  980  (20)  1",
                 "  400  3.0  0.0  0.0  q000002  1  500  (0)  C  CICH1  DNA/hAT  (500)  500  1  2"]
        self.assertEqual(fr.parse_hits(lines), [("q000001", "CICH1", 98.0, 1, 980), ("q000002", "CICH1", 97.0, 1, 500)])

    def test_flynn_categories(self):
        perfect = [("q1", "F", 98.0, 1, 980)]
        good = [("q1", "F", 99.0, 1, 500), ("q2", "F", 99.0, 480, 1000)]
        present = [("q1", "F", 85.0, 1, 900)]
        absent = [("q1", "F", 99.0, 1, 600)]
        self.assertEqual(fr.classify(perfect, 1000)["category"], "perfect")
        self.assertEqual(fr.classify(good, 1000)["category"], "good")
        self.assertEqual(fr.classify(present, 1000)["category"], "present")
        self.assertEqual(fr.classify(absent, 1000)["category"], "not found")
        self.assertEqual(fr.classify([], 1000)["category"], "not found")
        best = fr.classify(good, 1000)
        self.assertEqual((best["sequences_80"], round(best["coverage_95"], 3)), (2, 1.0))


class MergedVersionLookup(unittest.TestCase):
    def test_a_merged_version_uses_its_own_library(self):
        import tempfile
        root = Path(tempfile.mkdtemp())
        lib = root / "RI_v4" / "library"
        lib.mkdir(parents=True)
        (lib / "RI_v4-families.prefixed.fa").write_text(">A1#LTR/Gypsy\nACGT\n")
        (lib / "merge_record.json").write_text("{}")
        saved = c.version_folder
        c.version_folder = lambda version: root / version
        try:
            build_id, original, prefixed = c.library_files("RI_v4")
        finally:
            c.version_folder = saved
        self.assertEqual((build_id, original, prefixed), ("RI_v4", lib / "RI_v4-families.prefixed.fa", lib / "RI_v4-families.prefixed.fa"))


class MergeEffect(unittest.TestCase):
    def test_lost_and_gained_against_the_union_of_sources(self):
        merged = [("s", 0, 40, "m1", "LTR/Gypsy", 500), ("s", 90, 100, "m2", "DNA/hAT", 100)]
        sources = [[("s", 0, 50, "a1", "LTR/Gypsy", 400)], [("s", 30, 60, "b1", "DNA/hAT", 300)]]
        e = em.merge_effect(merged, sources)
        self.assertEqual((e["merged"], e["union"], e["both"], e["lost"], e["gained"]), (50, 60, 40, 20, 10))
        self.assertEqual(e["lost_by_class"].get("LTR"), 10)    # 40-50 is LTR in the sources (score 400 beats 300)
        self.assertEqual(e["lost_by_class"].get("DNA"), 10)    # 50-60
        self.assertEqual(e["gained_by_class"].get("DNA"), 10)  # 90-100


GFF = [
    "##gff-version 3\n",
    "s1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene-g1;gene=g1\n",
    "s1\tGnomon\tCDS\t11\t100\t.\t+\t0\tID=cds-1;Parent=rna-1;gene=g1;product=myosin%2C heavy chain\n",
    "s1\tGnomon\tCDS\t51\t150\t.\t+\t0\tID=cds-2;Parent=rna-2;gene=g1;product=myosin%2C heavy chain\n",
    "s1\tGnomon\tCDS\t301\t400\t.\t+\t0\tID=cds-3;Parent=rna-3;gene=t1;product=transposase-like protein\n",
    "s1\tGnomon\texon\t1\t500\t.\t+\t.\tID=exon-1;Parent=rna-1;gene=g1\n",
]


class GeneOverlap(unittest.TestCase):
    def setUp(self):
        self.genes = go.read_cds(GFF)

    def test_cds_merged_across_transcripts_and_products_decoded(self):
        self.assertEqual(self.genes["g1"]["spans"], {"s1": [(10, 150)]})
        self.assertEqual(self.genes["g1"]["product"], "myosin, heavy chain")
        self.assertFalse(self.genes["g1"]["te_like"])
        self.assertTrue(self.genes["t1"]["te_like"])

    def test_masked_cds_by_class_family_and_gene(self):
        rows = [("s1", 0, 40, "fA", "Unknown"), ("s1", 140, 160, "fB", "Simple_repeat"), ("s1", 300, 400, "fC", "DNA/hAT")]
        classes, families, per_gene = go.summarise(self.genes, rows)
        total, by_class = classes["genes without TE products"]
        self.assertEqual((total, by_class), (140, {"Unknown": 30, "Other": 10}))
        self.assertEqual(classes["genes with TE products"], (100, {"DNA": 100}))
        self.assertEqual(families[("fA", "Unknown")], {"bp": 30, "genes": ["g1"], "named": ["g1"], "te_genes": 0})
        self.assertEqual(families[("fC", "DNA/hAT")]["te_genes"], 1)
        self.assertEqual(per_gene["g1"], ("myosin, heavy chain", False, 140, 40, ("fA", "Unknown")))

    def test_te_product_words(self):
        for product, expected in (("reverse transcriptase", True), ("Pol polyprotein", True), ("piggyBac transposable element-derived protein 4", True),
                                  ("hatching enzyme", False), ("polymerase delta", False), ("integrin beta-1", False)):
            self.assertEqual(bool(go.TE_PRODUCT.search(product)), expected, product)


SEQS_EMBL = """CC   header
XX
ID   ACROBAT1     repeatmasker; DNA;  ???;  12 BP.
CC   consensus
XX
SQ   Sequence 12 BP;
     gggcccggtt tt                                                       12
//
ID   ALIEN1_CA     repeatmasker; DNA;  ???;  8 BP.
SQ   Sequence 8 BP;
     acgtacgt                                                            8
//
""".splitlines(True)
META_EMBL = """ID   ACROBAT1     repeatmasker; DNA;  ???;  12 BP.
KW   DNA/PIF.
CC   RepeatMasker Annotations:
CC        Type: DNA
CC        SubType: PIF
CC        Species: Danio
CC        SearchStages: 
//
ID   ALIEN1_CA     repeatmasker; DNA;  ???;  8 BP.
CC   RepeatMasker Annotations:
CC        Type: Satellite
CC        SubType: 
CC        Species: 
//
""".splitlines(True)


class RepbaseLibrary(unittest.TestCase):
    def test_the_two_halves_are_joined_by_name(self):
        entries = rbm.build_repbase(SEQS_EMBL, META_EMBL)
        self.assertEqual(entries, [("ACROBAT1", "DNA", "PIF", "Danio", "GGGCCCGGTTTT"), ("ALIEN1_CA", "Satellite", "", "", "ACGTACGT")])
        self.assertEqual([rbm.label_of(e[1], e[2]) for e in entries], ["DNA/PIF", "Satellite"])

    def test_a_sequence_without_metadata_stops(self):
        with self.assertRaises(SystemExit):
            rbm.build_repbase(SEQS_EMBL, META_EMBL[:8])

    def test_support_needs_a_transposon_class(self):
        self.assertIn("DNA", rbm.TE_CLASSES)
        self.assertNotIn("Satellite", rbm.TE_CLASSES)


class GeneFilter(unittest.TestCase):
    GENOME = {"s": "AAATGGCTTTGCCCAAAAGGGTAGTTT"}   # ATG GCT TTG | CCC (intron AAAA) GGG TAG

    def test_plus_strand_two_exons(self):
        t = {"sequence": "s", "strand": "+", "segments": [(2, 14, 0), (18, 24, 0)]}
        self.assertEqual(gf.protein_of(t, self.GENOME), "MALPG")

    def test_minus_strand_and_phase(self):
        fwd = "ATGGCTTAA"                      # M A *
        genome = {"s": "CC" + gf.revcomp(fwd) + "GG"}
        t = {"sequence": "s", "strand": "-", "segments": [(2, 11, 0)]}
        self.assertEqual(gf.protein_of(t, genome), "MA")
        t2 = {"sequence": "s", "strand": "-", "segments": [(2, 11, 1)]}   # skip one base: TGG CTT AA
        self.assertEqual(gf.protein_of(t2, genome), "WL")

    def test_transposon_products_and_short_proteins_are_left_out(self):
        genome = {"s": "ATG" + "GCT" * 40 + "TAA"}
        tx = {"t1": {"sequence": "s", "strand": "+", "gene": "g1", "product": "myosin", "segments": [(0, 126, 0)]},
              "t2": {"sequence": "s", "strand": "+", "gene": "g2", "product": "transposase", "segments": [(0, 126, 0)]},
              "t3": {"sequence": "s", "strand": "+", "gene": "g3", "product": "tiny", "segments": [(0, 30, 0)]},
              "t4": {"sequence": "s", "strand": "+", "gene": "LOC1", "product": "uncharacterized protein LOC1", "segments": [(0, 126, 0)]}}
        proteins, skipped = gf.host_proteins(tx, genome)
        self.assertEqual(list(proteins), ["g1"])
        self.assertEqual(skipped, {"transposon product": 1, "internal stop or under 30 aa": 1, "uncharacterized or unnamed": 1})

    def test_evidence_a_needs_every_genome(self):
        genes = {"G1": [{"te_product": "no", "percent_masked": "80", "top_family": "f1", "product": "NACHT, LRR and PYD domains-containing protein 3"}] * 5,
                 "G2": [{"te_product": "no", "percent_masked": "80", "top_family": "f1", "product": "NLR family CARD domain-containing protein 3"}] * 5}
        fams = {"G1": [{"family": "f2", "other_genes_hit": "40", "named_genes_hit": "25"}, {"family": "f3", "other_genes_hit": "60", "named_genes_hit": "21"}],
                "G2": [{"family": "f2", "other_genes_hit": "40", "named_genes_hit": "19"}, {"family": "f3", "other_genes_hit": "60", "named_genes_hit": "22"}]}
        a = gf.evidence_a(genes, fams, 5, 20)
        self.assertTrue(a["f1"]["flag"])
        self.assertEqual(a["f1"]["by"], "gene families")
        self.assertFalse(a["f2"]["flag"])          # 19 named genes on G2; uncharacterized models do not count
        self.assertEqual(a["f2"]["genes_hit"], [25, 19])
        self.assertEqual((a["f3"]["flag"], a["f3"]["by"]), (True, "named genes"))

    def test_evidence_a_needs_the_named_gene_count(self):
        with self.assertRaises(SystemExit):
            gf.evidence_a({"G1": []}, {"G1": [{"family": "f2", "other_genes_hit": "25"}]}, 5, 20)

    def test_strong_and_weak_evidence_b(self):
        self.assertEqual(gf.b_strength(("h", 1e-30, 300.0, 0.45), "early endosome antigen 1", 0.30), "strong")
        self.assertEqual(gf.b_strength(("h", 1e-30, 300.0, 0.10), "protein NLRC3", 0.30), "strong")
        self.assertEqual(gf.b_strength(("h", 1e-30, 300.0, 0.10), "trace amine-associated receptor 13c", 0.30), "strong")
        self.assertEqual(gf.b_strength(("h", 1e-30, 300.0, 0.10), "early endosome antigen 1", 0.30), "weak")

    def test_a_report_from_an_earlier_rule_is_stale(self):
        import tempfile
        d = Path(tempfile.mkdtemp())
        (d / "old.tsv").write_text("family\tevidence_A\tevidence_B\n")
        (d / "new.tsv").write_text("family\tevidence_A\tevidence_B\tevidence_B_strength\ttag_category\n")
        self.assertTrue(gf.report_is_stale(d / "old.tsv"))
        self.assertFalse(gf.report_is_stale(d / "new.tsv"))
        self.assertFalse(gf.report_is_stale(d / "absent.tsv"))

    def test_evidence_b_compares_with_transposon_proteins(self):
        host = gf.best_hits(["q1 h1 1e-30 300", "q1 h2 1e-20 200", "q2 h1 1e-30 300", "q3 h3 1e-8 120"])
        te = gf.best_hits(["q2 T1 1e-40 200"])
        self.assertEqual(host["q1"], ("h1", 1e-30, 300.0, None))
        self.assertEqual(gf.evidence_b(host, te, 1e-10, 2.0), {"q1"})   # q2: 300 < 2 x 200; q3: evalue too weak

    def test_host_match_coverage_on_both_strands(self):
        hits = gf.best_hits(["q1 h1 1e-30 300 101 400 1000", "q2 h2 1e-20 200 500 201 1000"])
        self.assertAlmostEqual(hits["q1"][3], 0.30)
        self.assertAlmostEqual(hits["q2"][3], 0.30)

    def test_files_of_the_first_method_are_found(self):
        import tempfile
        work = Path(tempfile.mkdtemp())
        for name in ("host_proteins.fa", "blastx_host.tsv", "host_db.pdb", "host_named.fa", "filter_report.tsv", "blastx_te.tsv"):
            (work / name).write_text("")
        self.assertEqual(sorted(p.name for p in gf.stale_files(work)), ["blastx_host.tsv", "host_db.pdb", "host_proteins.fa"])

    def test_decision_rule(self):
        self.assertEqual(gf.decide(True, "strong", True), ("tag", "A and B"))
        self.assertEqual(gf.decide(True, "", True), ("keep", "A only; Repbase transposon support"))
        self.assertEqual(gf.decide(True, "weak", False), ("tag", "A only; no Repbase support"))
        self.assertEqual(gf.decide(False, "strong", False), ("tag", "B only; no Repbase support"))
        self.assertEqual(gf.decide(False, "weak", False), ("keep", "weak B only (small host match, not a gene family): not tagged"))
        self.assertEqual(gf.decide(False, "", False), ("keep", ""))

    def test_tag_categories(self):
        cases = {"protein NLRC3": "NLR", "tripartite motif-containing protein 16": "TRIM", "olfactory receptor 10Z1": "OR",
                 "trace amine-associated receptor 13c": "TAAR", "protocadherin alpha-2": "PCDH", "zinc finger protein 668": "ZNF",
                 "histone H1.0": "HIST", "histone-lysine N-methyltransferase SETD2": "OTHER", "tubulin alpha-1A chain": "TUB",
                 "tubulin-specific chaperone A": "OTHER", "homeodomain-interacting protein kinase 1": "OTHER"}
        for product, expected in cases.items():
            self.assertEqual(gf.category_of(product), expected, product)
        self.assertFalse(gf.HOST_FAMILIES.search("histone-lysine N-methyltransferase SETD2"))

    def test_tag_category_prefers_a_known_family(self):
        self.assertEqual(gf.tag_category("protein NLRC3", "strong", {"flag": True, "category": "TRIM"}), "NLR")
        self.assertEqual(gf.tag_category("myosin-7", "strong", {"flag": True, "category": "TRIM"}), "TRIM")
        self.assertEqual(gf.tag_category("myosin-7", "strong", {"flag": False}), "OTHER")
        self.assertEqual(gf.tag_category("", "", {"flag": True, "category": "MULTI"}), "MULTI")
        self.assertEqual(gf.tagged_name("MZ_MZ4f_rnd-1_family-163", "NLR", "LTR/Gypsy"), "MZ_MZ4f_rnd-1_family-163_hostgene-NLR#LTR/Gypsy")

    def test_a_filtered_version_uses_its_own_library(self):
        import tempfile
        root = Path(tempfile.mkdtemp())
        lib = root / "RI_v6" / "library"
        lib.mkdir(parents=True)
        (lib / "RI_v6-families.prefixed.fa").write_text(">A1#LTR/Gypsy\nACGT\n")
        (lib / "filter_record.json").write_text("{}")
        saved = c.version_folder
        c.version_folder = lambda version: root / version
        try:
            self.assertEqual(c.library_files("RI_v6")[0], "RI_v6")
        finally:
            c.version_folder = saved


class OneLabelPerBase(unittest.TestCase):
    def test_simple_repeat_inside_a_transposon_keeps_the_transposon_label(self):
        rows = [("s", 0, 100, "f1", "LTR/Gypsy", 500), ("s", 40, 60, "(CA)n", "Simple_repeat", 30)]
        self.assertEqual(la.one_label_per_base(rows), [("s", 0, 100, "f1", "LTR/Gypsy")])

    def test_overlapping_hits_split_at_the_higher_score(self):
        rows = [("s", 0, 60, "fa", "DNA/hAT", 300), ("s", 50, 100, "fb", "LTR/Gypsy", 400)]
        self.assertEqual(la.one_label_per_base(rows), [("s", 0, 50, "fa", "DNA/hAT"), ("s", 50, 100, "fb", "LTR/Gypsy")])

    def test_equal_scores_go_to_the_longer_annotation_and_gaps_stay_unmasked(self):
        rows = [("s", 0, 30, "short", "DNA", 200), ("s", 10, 80, "long", "LINE/L2", 200), ("s", 90, 95, "g", "Unknown", 50)]
        self.assertEqual(la.one_label_per_base(rows),
                         [("s", 0, 10, "short", "DNA"), ("s", 10, 80, "long", "LINE/L2"), ("s", 90, 95, "g", "Unknown")])

    def test_no_base_is_counted_twice_in_the_matrix(self):
        a = la.one_label_per_base([("s", 0, 100, "f1", "LTR/Gypsy", 500), ("s", 40, 60, "(CA)n", "Simple_repeat", 30)])
        b = la.one_label_per_base([("s", 0, 100, "g1", "DNA/hAT", 400)])
        m = la.class_matrix(a, b)
        self.assertEqual(la.agreement(m), (100, 0))
        self.assertNotIn(("Other", "DNA"), m)


class MaskingGain(unittest.TestCase):
    """On a 100-bp genome: v1 masks 0-30 (DNA), v2 20-50 (LTR), v3 40-60 (Unknown); together 0-60."""

    def setUp(self):
        sets = {"v1": cm.masked_sets([("s", 0, 30, "f1", "DNA")]),
                "v2": cm.masked_sets([("s", 20, 50, "f2", "LTR")]),
                "v3": cm.masked_sets([("s", 40, 60, "f3", "Unknown")])}
        self.rows = mg.gain_rows("G", sets, 100)

    def steps(self, order):
        return [(r[5], r[6]) for r in self.rows if r[1] == order and r[4] == "ALL"]

    def test_each_step_adds_only_new_bases(self):
        self.assertEqual(self.steps("v1 > v2 > v3"), [(30, 30), (20, 50), (10, 60)])
        self.assertEqual(self.steps("v3 > v2 > v1"), [(20, 20), (20, 40), (20, 60)])

    def test_every_order_ends_at_the_same_total(self):
        finals = {r[6] for r in self.rows if r[4] == "ALL" and r[2] == 3}
        self.assertEqual(finals, {60})
        self.assertEqual(len({r[1] for r in self.rows}), 6)

    def test_step_means_add_up_to_the_total(self):
        means = mg.step_means(self.rows)
        self.assertAlmostEqual(means[1], 80 / 3, places=2)   # first library: 30, 30 or 20 bp
        self.assertAlmostEqual(sum(means.values()), 60.0, places=2)

    def test_added_bases_carry_the_added_librarys_class(self):
        row = next(r for r in self.rows if r[1] == "v1 > v2 > v3" and r[2] == 2 and r[4] == "LTR")
        self.assertEqual(row[5], 20)


MaskingOverlapRows = [["G", "A vs B", "A", "ALL", 20, 10, 100, "20.000", "10.000"],
                      ["G", "A vs B", "B", "ALL", 20, 10, 100, "20.000", "10.000"]]


if __name__ == "__main__":
    unittest.main()
