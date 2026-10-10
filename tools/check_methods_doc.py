#!/usr/bin/env python3
"""Keep docs/methods.md in step with the analysis.

  python3 tools/check_methods_doc.py coverage
      every analysis script (scripts/**, R/lib/*.R) is named in docs/methods.md
  python3 tools/check_methods_doc.py values
      every value in the methods document's "Values that come from the data" table equals the
      value in outputs/extended/tables/manuscript_numbers.csv (skipped if that file is absent)
  python3 tools/check_methods_doc.py staged <commit-message-file>
      (git commit-msg hook) if the commit changes analysis files, docs/methods.md must be in the
      same commit, or the message must contain a line "Methods-unchanged: <reason>"
  python3 tools/check_methods_doc.py range <base> <head>
      (CI) the same rule for every commit in base..head

Analysis files: scripts/**, R/lib/** (except make_paths.R), resources/**, env/renv.lock.
Standard library only. Run from the repository root. Exit status 1 on failure.
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "methods.md"
NUMBERS = ROOT / "outputs" / "extended" / "tables" / "manuscript_numbers.csv"
TRAILER = re.compile(r"^Methods-unchanged:\s*\S", re.MULTILINE)
VALUES_HEADING = "## Values that come from the data"


def is_analysis_file(path):
    p = path.replace("\\", "/")
    if p == "R/lib/make_paths.R":
        return False
    return p.startswith(("scripts/", "R/lib/", "resources/")) or p == "env/renv.lock"


def analysis_scripts():
    files = [p for p in (ROOT / "scripts").rglob("*") if p.suffix in (".R", ".sh") and p.is_file()]
    files += [p for p in (ROOT / "R" / "lib").glob("*.R") if p.name != "make_paths.R"]
    return sorted(str(p.relative_to(ROOT)) for p in files)


def coverage():
    text = DOC.read_text(encoding="utf-8")
    missing = []
    for f in analysis_scripts():
        short = f.split("/", 1)[1] if f.startswith("scripts/") else f   # "01_gse204943/01_minfi_preprocess.R"
        if short not in text and f not in text:
            missing.append(f)
    if missing:
        print("docs/methods.md does not describe these scripts:\n  " + "\n  ".join(missing))
        return 1
    print(f"coverage: all {len(analysis_scripts())} analysis scripts are described in docs/methods.md")
    return 0


def values():
    if not NUMBERS.exists():
        print(f"values: {NUMBERS.relative_to(ROOT)} not found; run `make audit` first (check skipped)")
        return 0
    current = {}
    with NUMBERS.open(encoding="utf-8") as h:
        for row in csv.DictReader(h):
            current[row["statement"]] = row["value"]
    text = DOC.read_text(encoding="utf-8")
    if VALUES_HEADING not in text:
        print(f"docs/methods.md has no '{VALUES_HEADING}' section")
        return 1
    section = text.split(VALUES_HEADING, 1)[1].split("\n## ", 1)[0]
    rows = [r for r in section.splitlines() if r.startswith("|") and not re.match(r"^\|\s*-", r)][1:]
    bad, n = [], 0
    for r in rows:
        cells = [c.strip() for c in r.strip().strip("|").split("|")]
        if len(cells) < 2:
            continue
        statement, quoted = cells[0].strip("`"), cells[1].strip("`")
        n += 1
        if statement not in current:
            bad.append(f"statement not in manuscript_numbers.csv: {statement}")
        elif current[statement] != quoted:
            bad.append(f"{statement}\n      methods.md: {quoted}\n      current:    {current[statement]}")
    if bad:
        print("docs/methods.md quotes values that no longer match the current results:\n  " + "\n  ".join(bad))
        print("Update docs/methods.md (and the text that uses these values), then rerun this check.")
        return 1
    print(f"values: all {n} data-derived values in docs/methods.md match manuscript_numbers.csv")
    return 0


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout


def rule(files, message, label):
    changed = [f for f in files if is_analysis_file(f)]
    if not changed or "docs/methods.md" in files or TRAILER.search(message):
        return True
    print(f"{label}: analysis files changed without updating docs/methods.md:\n  " + "\n  ".join(changed))
    print("Update docs/methods.md in the same commit, or, if the methods really are unchanged (e.g. a\n"
          "comment, plot styling or refactor), add a line to the commit message:\n"
          "  Methods-unchanged: <reason>")
    return False


def staged(msg_file):
    files = git("diff", "--cached", "--name-only").split()
    message = Path(msg_file).read_text(encoding="utf-8")
    return 0 if rule(files, message, "commit") else 1


def commit_range(base, head):
    ok = True
    for sha in git("rev-list", "--no-merges", f"{base}..{head}").split():
        files = git("diff-tree", "--no-commit-id", "--name-only", "-r", "--root", sha).split()
        ok &= rule(files, git("log", "-1", "--format=%B", sha), f"commit {sha[:8]}")
    return 0 if ok else 1


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "coverage":
        sys.exit(coverage())
    if cmd == "values":
        sys.exit(values())
    if cmd == "staged" and len(sys.argv) == 3:
        sys.exit(staged(sys.argv[2]) or coverage())
    if cmd == "range" and len(sys.argv) == 4:
        sys.exit(commit_range(sys.argv[2], sys.argv[3]) or coverage())
    print(__doc__)
    sys.exit(2)
