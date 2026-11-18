#!/usr/bin/env python3
"""Record the evidence that the alignment project's environment works (repository ga-alignment).

Checks, in order:
  host        host name, user, CPU threads, load average
  avx2        CPU threads whose flags include avx2 (the Cactus release images need it)
  disk        free space in the home, scratch and containers folders
  python      the interpreter running this script; the conda environment, exported to a file
  ga_common   the coordinator's library imports, and Paths.load() reads its configuration
  apptainer   version, from the tool's own output
  minimap2    version, its usage text, and its manual page if the environment has one
  cactus      the image file's size and sha256; the help texts of cactus-pangenome and of the
              stages cactus-graphmap, -graphmap-split and -graphmap-join; whether the options the
              plan uses appear in them; whether the configuration attributes the plan relies on
              appear in the image's default config

A dry run unless --go. With --go it writes one new folder, <scratch>/env/<UTC timestamp>/, with
report.tsv, report.md and the captured texts. It never overwrites: if that folder exists it stops.
It writes nothing else, installs nothing, and runs nothing but the tools' version and help commands.
Exit status: 0 if no check failed, 1 if one did, 2 if it refused to overwrite.
"""

import argparse
import datetime
import getpass
import hashlib
import os
import re
import shutil
import socket
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SIF = Path.home() / "containers" / "cactus_v3.3.0.sif"
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"

# Options the plan uses (proposal, section 4), looked up in the installed version's help texts.
PLANNED_OPTIONS = [
    "--reference", "--refContigs", "--otherContig", "--gfa", "--clip", "--filter",
    "--maxCores", "--mgCores", "--mapCores", "--consCores", "--indexCores", "--maxMemory",
    "--permissiveContigFilter", "--delFilter",
]
# Configuration attributes described in the uploaded Minigraph-Cactus documentation.
CONFIG_ATTRIBUTES = ["zipAlleles", "GAFOverlapFilterExact", "minigraphSortInput"]
CONFIG_NAME = "cactus_progressive_config.xml"
# The all-in-one command first (required), then the stages whose options the plan uses.
HELP_TOOLS = ["cactus-pangenome", "cactus-graphmap", "cactus-graphmap-split", "cactus-graphmap-join"]


# ---------------------------------------------------------------- helpers (pure, tested)

def count_avx2(cpuinfo_text):
    """Number of 'flags' lines (one per CPU thread) that list avx2 as a whole word."""
    return sum(1 for line in cpuinfo_text.splitlines()
               if line.startswith("flags") and re.search(r"(?<!\S)avx2(?!\S)", line))


def option_present(option, text):
    """True if the option appears as a whole option name (so --clip does not match --clipXyz)."""
    return re.search(r"(?<![\w-])" + re.escape(option) + r"(?![\w-])", text) is not None


def attributes_present(attributes, config_text):
    """For each attribute, the config lines that mention it (an empty list if none)."""
    found = {}
    for attribute in attributes:
        pattern = re.compile(r"\b" + re.escape(attribute) + r"\b")
        found[attribute] = [line.strip() for line in config_text.splitlines() if pattern.search(line)]
    return found


def interpreter_from_shebang(first_line):
    """The interpreter named on a script's first line, or None if it is not a #! line."""
    if not first_line.startswith("#!"):
        return None
    words = first_line[2:].strip().split()
    if not words:
        return None
    if os.path.basename(words[0]) == "env" and len(words) > 1:
        return words[-1]
    return words[0]


def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def first_line(text):
    for line in text.splitlines():
        if line.strip():
            return line.strip()
    return ""


def nearest_existing(path):
    path = Path(path).expanduser()
    while not path.exists() and path != path.parent:
        path = path.parent
    return path


def run(cmd, timeout=600):
    """Run a command; return (returncode, stdout, stderr). A missing program gives 127."""
    try:
        done = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return done.returncode, done.stdout, done.stderr
    except FileNotFoundError as error:
        return 127, "", str(error)
    except subprocess.TimeoutExpired:
        return 124, "", "timed out after {} s".format(timeout)


# ---------------------------------------------------------------- the checks

def check_ga_common():
    try:
        import ga_common  # noqa: F401
    except Exception as error:  # report any import problem, do not crash
        return "fail", "import failed: {}".format(error)
    where = getattr(ga_common, "__file__", "?")
    version = "?"
    try:
        from importlib import metadata
        for name in ("ga_common", "ga-common"):
            try:
                version = metadata.version(name)
                break
            except metadata.PackageNotFoundError:
                continue
    except Exception:
        pass
    try:
        from ga_common.paths import Paths
        Paths.load()
    except Exception as error:
        return "fail", "{} (version {}); Paths.load() failed: {}".format(where, version, error)
    return "ok", "{} (version {}); Paths.load() read the configuration".format(where, version)


def find_cactus_config(sif, runner):
    """Locate the image's default configuration file and return (path, text, problem)."""
    base = ["apptainer", "exec", "--cleanenv", str(sif)]
    python = "python3"
    rc, out, _ = runner(base + ["sh", "-c", 'head -n 1 "$(command -v cactus-pangenome)"'])
    if rc == 0:
        python = interpreter_from_shebang(first_line(out)) or python
    rc, out, _ = runner(base + [python, "-c",
                                "import os, cactus; print(os.path.dirname(cactus.__file__))"])
    candidates = []
    if rc == 0 and out.strip():
        candidates.append(out.strip().splitlines()[-1] + "/" + CONFIG_NAME)
    for candidate in candidates:
        rc, text, _ = runner(base + ["cat", candidate])
        if rc == 0 and text:
            return candidate, text, None
    rc, out, _ = runner(base + ["sh", "-c",
                                "find / -xdev -name {} 2>/dev/null | head -n 1".format(CONFIG_NAME)],
                        timeout=900)
    if rc == 0 and out.strip():
        candidate = out.strip().splitlines()[0]
        rc, text, _ = runner(base + ["cat", candidate])
        if rc == 0 and text:
            return candidate, text, None
    return None, "", "could not locate {} in the image".format(CONFIG_NAME)


def collect(args, runner, outdir):
    """Run every check. Returns the report rows; writes captured texts into outdir."""
    rows = []

    def add(check, status, value):
        rows.append((check, status, " ".join(str(value).split())))

    def save(name, text):
        (outdir / name).write_text(text)
        return name

    # host
    add("host", "info", "{} as {}; {} CPU threads".format(socket.gethostname(), getpass.getuser(),
                                                          os.cpu_count()))
    try:
        add("load average", "info", "{:.1f} {:.1f} {:.1f} (1, 5, 15 min)".format(*os.getloadavg()))
    except OSError:
        add("load average", "info", "not available")

    # avx2
    try:
        threads = count_avx2(Path(args.cpuinfo).read_text())
        add("avx2", "ok" if threads > 0 else "fail",
            "{} CPU threads list avx2".format(threads) if threads > 0 else
            "no CPU thread lists avx2: the Cactus release images will not run (legacy binaries would)")
    except OSError as error:
        add("avx2", "fail", "could not read {}: {}".format(args.cpuinfo, error))

    # disk
    for label, path in (("home", Path.home()), ("scratch", args.scratch),
                        ("containers", Path(args.sif).expanduser().parent)):
        existing = nearest_existing(path)
        usage = shutil.disk_usage(str(existing))
        add("disk " + label, "info", "{:.0f} GB free of {:.0f} GB at {}".format(
            usage.free / 1e9, usage.total / 1e9, existing))

    # python and environment
    add("python", "info", "{} at {}; conda env {}".format(
        sys.version.split()[0], sys.executable, os.environ.get("CONDA_DEFAULT_ENV", "none")))
    conda_exe, env_name = os.environ.get("CONDA_EXE"), os.environ.get("CONDA_DEFAULT_ENV")
    if conda_exe and env_name:
        rc, out, err = runner([conda_exe, "env", "export", "-n", env_name])
        add("conda env export", "info" if rc == 0 else "warn",
            "saved to " + save("conda_env_export.yml", out) if rc == 0 else first_line(err) or rc)
    else:
        add("conda env export", "info", "not run inside a conda environment")
    add("throttling", "info", "--threads {} --max-jobs {} (recorded; this script is single-threaded)".format(
        args.threads, args.max_jobs))

    # ga_common
    status, value = check_ga_common()
    add("ga_common", status, value)

    # apptainer
    rc, out, err = runner(["apptainer", "--version"])
    add("apptainer", "ok" if rc == 0 else "fail", first_line(out) if rc == 0 else first_line(err) or rc)

    # minimap2
    rc, out, err = runner(["minimap2", "--version"])
    add("minimap2", "ok" if rc == 0 else "fail", first_line(out) if rc == 0 else first_line(err) or rc)
    rc, out, err = runner(["minimap2"])  # without arguments it prints its usage, including the presets
    usage = out + err
    add("minimap2 usage", "info" if usage.strip() else "warn",
        "saved to " + save("minimap2_usage.txt", usage) if usage.strip() else "no usage text")
    manpage = Path(os.environ.get("CONDA_PREFIX", "/nonexistent")) / "share" / "man" / "man1" / "minimap2.1"
    if manpage.is_file():
        add("minimap2 manual page", "info", "saved to " + save("minimap2_manpage.1", manpage.read_text()))
    else:
        add("minimap2 manual page", "info", "not in this environment")

    # cactus image
    sif = Path(args.sif).expanduser()
    if not sif.is_file():
        add("cactus image", "fail", "{} not found".format(sif))
        return rows
    add("cactus image", "ok", "{} ({:.2f} GB)".format(sif, sif.stat().st_size / 1e9))
    add("cactus image sha256", "info", sha256_file(sif))

    helps = {}
    for tool in HELP_TOOLS:
        rc, out, err = runner(["apptainer", "exec", "--cleanenv", str(sif), tool, "--help"])
        text = out + err
        ok = rc == 0 and bool(text.strip())
        if ok:
            helps[tool] = text
        add(tool + " --help", "ok" if ok else ("fail" if tool == "cactus-pangenome" else "warn"),
            "saved to " + save(tool + "_help.txt", text) if ok else first_line(err) or rc)
    for option in PLANNED_OPTIONS:
        where = [tool for tool, text in helps.items() if option_present(option, text)]
        add("option " + option, "ok" if where else "warn",
            "in " + ", ".join(where) if where else "in no help text: do not use it as written")

    path, text, problem = find_cactus_config(sif, runner)
    if problem:
        add("cactus config", "warn", problem)
    else:
        add("cactus config", "ok", "{} (saved to {})".format(path, save("cactus_config.xml", text)))
        for attribute, lines in attributes_present(CONFIG_ATTRIBUTES, text).items():
            add("config " + attribute, "ok" if lines else "warn",
                lines[0][:160] if lines else
                "absent from this image's default config: decide the Cactus version (proposal, section 7)")
    return rows


# ---------------------------------------------------------------- reports

def write_reports(outdir, rows, args, version):
    with open(outdir / "report.tsv", "w") as handle:
        handle.write("check\tstatus\tvalue\n")
        for check, status, value in rows:
            handle.write("{}\t{}\t{}\n".format(check, status, value))
    lines = ["# Environment report ({})".format(version), "",
             "- written: {} UTC".format(outdir.name),
             "- image: {}".format(Path(args.sif).expanduser()),
             "", "| check | status | value |", "|---|---|---|"]
    lines += ["| {} | {} | {} |".format(c, s, v.replace("|", "\\|")) for c, s, v in rows]
    problems = [(c, s, v) for c, s, v in rows if s in ("fail", "warn")]
    lines += ["", "## Needs attention" if problems else "## Nothing needs attention", ""]
    lines += ["- **{}** ({}): {}".format(c, s, v) for c, s, v in problems]
    (outdir / "report.md").write_text("\n".join(lines) + "\n")


def read_version():
    try:
        return (REPO / "VERSION").read_text().strip()
    except OSError:
        return "unknown version"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--sif", type=Path, default=DEFAULT_SIF,
                        help="the Cactus image (default: %(default)s)")
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH,
                        help="the project's scratch folder (default: %(default)s)")
    parser.add_argument("--cpuinfo", type=Path, default=Path("/proc/cpuinfo"), help=argparse.SUPPRESS)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("GA_MAX_THREADS") or 0),
                        help="the coordinator's throttle (recorded; 0 = no limit)")
    parser.add_argument("--max-jobs", type=int, default=int(os.environ.get("GA_MAX_JOBS") or 0),
                        help="the coordinator's throttle (recorded; 0 = no limit)")
    parser.add_argument("--go", action="store_true", help="run the checks and write the report")
    return parser.parse_args(argv)


def main(argv=None, runner=run, now=None):
    args = parse_args(argv)
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y%m%d-%H%M%S")
    outdir = Path(args.scratch).expanduser() / "env" / stamp
    version = read_version()
    if not args.go:
        print("check_env.py ({}): dry run.".format(version))
        print("With --go it would check: host and load, avx2, disk, python, ga_common, apptainer,")
        print("minimap2 (version, usage, manual page), and the Cactus image")
        print("  {}".format(Path(args.sif).expanduser()))
        print("(size, sha256, help texts, planned options, config attributes),")
        print("and write a new folder:\n  {}".format(outdir))
        print("Nothing has been run or written. Add --go to run the checks.")
        return 0
    if outdir.exists():
        print("STOP: {} already exists, and this script never overwrites.".format(outdir), file=sys.stderr)
        print("To redo it, remove it by hand:  rm -r {}".format(outdir), file=sys.stderr)
        return 2
    outdir.mkdir(parents=True)
    rows = collect(args, runner, outdir)
    write_reports(outdir, rows, args, version)
    failed = [c for c, s, _ in rows if s == "fail"]
    warned = [c for c, s, _ in rows if s == "warn"]
    print("Report written to {}".format(outdir))
    print("{} checks: {} failed, {} need attention.".format(len(rows), len(failed), len(warned)))
    for check in failed + warned:
        print("  - " + check)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
