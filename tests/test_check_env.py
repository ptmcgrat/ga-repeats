"""Tests for scripts/check_env.py. They use stand-in data only: no apptainer, minimap2, image or
ga_common is needed, so they run the same on the laptop and on mzebra."""

import contextlib
import datetime
import hashlib
import io
import sys
import tempfile
import types
import unittest
from unittest import mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_env as ce  # noqa: E402

CPUINFO = (
    "processor\t: 0\nflags\t\t: fpu vme sse4_2 avx avx2 fma\n\n"
    "processor\t: 1\nflags\t\t: fpu vme sse4_2 avx avx2 fma\n\n"
    "processor\t: 2\nflags\t\t: fpu vme sse4_2 avx\n"
)
CONFIG = ('<cactus_workflow_config>\n'
          '  <graphmap zipAlleles="1" GAFOverlapFilterExact="1" minigraphSortInput="mash"/>\n'
          '</cactus_workflow_config>\n')
PANGENOME_HELP = ("usage: cactus-pangenome jobStore seqFile --outDir --outName --reference "
                  "[--refContigs] [--otherContig] [--gfa] [--clip CLIP] [--filter] [--maxCores] "
                  "[--mgCores] [--mapCores] [--consCores] [--indexCores] [--maxMemory]\n")
GRAPHMAP_HELP = "usage: cactus-graphmap ... [--delFilter DELFILTER]\n"
SPLIT_HELP = "usage: cactus-graphmap-split ... [--permissiveContigFilter [PERMISSIVECONTIGFILTER]]\n"


def stand_in_runner(config=CONFIG):
    """Answers the commands check_env.py runs, with canned stand-in output."""
    def runner(cmd, timeout=600):
        line = " ".join(str(part) for part in cmd)
        if cmd == ["apptainer", "--version"]:
            return 0, "apptainer version 1.5.3 (stand-in)\n", ""
        if cmd[1:3] == ["env", "export"]:
            return 0, "name: {}\ndependencies: [stand-in]\n".format(cmd[-1]), ""
        if cmd == ["minimap2", "--version"]:
            return 0, "2.31-stand-in\n", ""
        if cmd == ["minimap2"]:
            return 1, "", "Usage: minimap2 [options] <target.fa> [query.fa] (stand-in)\n"
        if "head -n 1" in line:
            return 0, "#!/usr/bin/env python3\n", ""
        if "import os, cactus" in line:
            return 0, "/opt/stand-in/cactus\n", ""
        if line.endswith("cat /opt/stand-in/cactus/" + ce.CONFIG_NAME):
            return (0, config, "") if config else (1, "", "no such file")
        if "find / -xdev" in line:
            return 0, "", ""
        if line.endswith("cactus-pangenome --help"):
            return 0, PANGENOME_HELP, ""
        if line.endswith("cactus-graphmap --help"):
            return 0, GRAPHMAP_HELP, ""
        if line.endswith("cactus-graphmap-split --help"):
            return 0, SPLIT_HELP, ""
        return 127, "", "not found (stand-in)"
    return runner


class Helpers(unittest.TestCase):
    def test_count_avx2(self):
        self.assertEqual(ce.count_avx2(CPUINFO), 2)

    def test_count_avx2_needs_the_whole_word(self):
        self.assertEqual(ce.count_avx2("flags\t: avx512f avx2x\n"), 0)

    def test_option_present(self):
        text = "  --clip CLIP\n  --clipNonMinigraph\n  --refContigs REFCONTIGS"
        self.assertTrue(ce.option_present("--clip", text))
        self.assertTrue(ce.option_present("--refContigs", text))
        self.assertFalse(ce.option_present("--filter", text))
        self.assertFalse(ce.option_present("--clip", "only --clipNonMinigraph here"))

    def test_attributes_present(self):
        found = ce.attributes_present(["zipAlleles", "missingAttribute"], CONFIG)
        self.assertEqual(len(found["zipAlleles"]), 1)
        self.assertEqual(found["missingAttribute"], [])

    def test_interpreter_from_shebang(self):
        self.assertEqual(ce.interpreter_from_shebang("#!/usr/bin/env python3"), "python3")
        self.assertEqual(ce.interpreter_from_shebang("#!/opt/venv/bin/python3"), "/opt/venv/bin/python3")
        self.assertIsNone(ce.interpreter_from_shebang("import os"))

    def test_sha256_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "blob"
            data = b"stand-in image " * 100000
            path.write_bytes(data)
            self.assertEqual(ce.sha256_file(path, chunk=4096), hashlib.sha256(data).hexdigest())


class Runs(unittest.TestCase):
    """End to end, with a stand-in image file, cpuinfo, tools and ga_common."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.sif = self.tmp / "containers" / "cactus_stand-in.sif"
        self.sif.parent.mkdir()
        self.sif.write_bytes(b"not a real image")
        self.cpuinfo = self.tmp / "cpuinfo"
        self.cpuinfo.write_text(CPUINFO)
        self.scratch = self.tmp / "scratch" / "ga-alignment"
        self.now = datetime.datetime(2026, 1, 2, 3, 4, 5, tzinfo=datetime.timezone.utc)
        # a stand-in ga_common, removed again after each test
        package = types.ModuleType("ga_common")
        package.__file__ = "/stand-in/ga_common/__init__.py"
        paths = types.ModuleType("ga_common.paths")

        class Paths:
            @classmethod
            def load(cls):
                return cls()
        paths.Paths = Paths
        package.paths = paths
        saved = {name: sys.modules.get(name) for name in ("ga_common", "ga_common.paths")}
        sys.modules["ga_common"] = package
        sys.modules["ga_common.paths"] = paths

        def restore():
            for name, module in saved.items():
                if module is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = module
        self.addCleanup(restore)

    def argv(self, *extra):
        return ["--sif", str(self.sif), "--scratch", str(self.scratch),
                "--cpuinfo", str(self.cpuinfo)] + list(extra)

    def run_main(self, argv, runner=None):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = ce.main(argv, runner=runner or stand_in_runner(), now=self.now)
        return code, out.getvalue(), err.getvalue()

    def report(self):
        folder = self.scratch / "env" / "20260102-030405"
        rows = [line.split("\t") for line in (folder / "report.tsv").read_text().splitlines()[1:]]
        return folder, {row[0]: (row[1], row[2]) for row in rows}

    def test_dry_run_writes_nothing(self):
        code, out, _ = self.run_main(self.argv())
        self.assertEqual(code, 0)
        self.assertIn("dry run", out)
        self.assertFalse(self.scratch.exists())

    def test_go_writes_the_report(self):
        code, out, _ = self.run_main(self.argv("--go"))
        self.assertEqual(code, 0, out)
        folder, rows = self.report()
        for name in ("report.md", "minimap2_usage.txt", "cactus-pangenome_help.txt",
                     "cactus-graphmap_help.txt", "cactus_config.xml"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertEqual(rows["avx2"][0], "ok")
        self.assertEqual(rows["ga_common"][0], "ok")
        self.assertEqual(rows["cactus image sha256"][1], hashlib.sha256(b"not a real image").hexdigest())
        self.assertEqual(rows["option --delFilter"], ("ok", "in cactus-graphmap"))
        self.assertEqual(rows["option --permissiveContigFilter"], ("ok", "in cactus-graphmap-split"))
        # a stage whose help cannot be read is flagged, but does not fail the run
        self.assertEqual(rows["cactus-graphmap-join --help"][0], "warn")
        self.assertEqual(rows["config zipAlleles"][0], "ok")

    def test_conda_environment_is_exported(self):
        stand_in = {"CONDA_EXE": "/stand-in/conda", "CONDA_DEFAULT_ENV": "ga-alignment"}
        with mock.patch.dict("os.environ", stand_in):
            code, _, _ = self.run_main(self.argv("--go"))
        self.assertEqual(code, 0)
        folder, rows = self.report()
        self.assertEqual(rows["conda env export"], ("info", "saved to conda_env_export.yml"))
        self.assertIn("name: ga-alignment", (folder / "conda_env_export.yml").read_text())

    def test_never_overwrites(self):
        self.run_main(self.argv("--go"))
        code, _, err = self.run_main(self.argv("--go"))
        self.assertEqual(code, 2)
        self.assertIn("never overwrites", err)

    def test_no_avx2_fails(self):
        self.cpuinfo.write_text("processor\t: 0\nflags\t\t: fpu sse4_2\n")
        code, _, _ = self.run_main(self.argv("--go"))
        self.assertEqual(code, 1)
        self.assertEqual(self.report()[1]["avx2"][0], "fail")

    def test_missing_config_attribute_is_flagged_not_failed(self):
        config = CONFIG.replace(' zipAlleles="1"', "")
        code, _, _ = self.run_main(self.argv("--go"), runner=stand_in_runner(config))
        self.assertEqual(code, 0)
        self.assertEqual(self.report()[1]["config zipAlleles"][0], "warn")

    def test_missing_pangenome_help_fails(self):
        def runner(cmd, timeout=600):
            if " ".join(cmd).endswith("cactus-pangenome --help"):
                return 1, "", "stand-in failure"
            return stand_in_runner()(cmd, timeout)
        code, _, _ = self.run_main(self.argv("--go"), runner=runner)
        self.assertEqual(code, 1)
        self.assertEqual(self.report()[1]["cactus-pangenome --help"][0], "fail")

    def test_missing_image_fails(self):
        self.sif.unlink()
        code, _, _ = self.run_main(self.argv("--go"))
        self.assertEqual(code, 1)
        self.assertEqual(self.report()[1]["cactus image"][0], "fail")


if __name__ == "__main__":
    unittest.main()
