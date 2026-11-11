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
