#!/usr/bin/env python3
"""Checks that the live docs name only files and modules that exist.

A name in backquotes that ends in a file extension of the tree (`ae_m.f90`,
`scripts/run_test.sh`) must be a file of the tree, by path or by name, and a
name ending in _m (`logmesh_m`) a module of src/ or tests/. Left out:
- the historical records of HISTORICAL, which name things as they were at
  their dates;
- names with a placeholder in capitals (`NAME.out`, `test_NAME.f90`) or a
  wildcard;
- the files of ALLOWED, which the tree does not hold on purpose.

usage: docrefs.py [ROOT]   (ROOT the source tree; default this script's ..)
"""

import pathlib
import re
import sys

# records of plans, audits and changes, as of their dates
HISTORICAL = {
    "changes_from_4.0.1.md",
    "consistency_audit.md",
    "fortran_best_practices.md",
    "fr_merge_plan.md",
    "library_plan.md",
    "numerics_audit.md",
    "release_notes.md",
    "smoothness_plan.md",
    "test_suite_speed.md",
}
# files named on purpose that the tree does not hold
ALLOWED = {
    "CMakeUserPresets.json",  # personal presets, ignored by git
    "cov.json",  # written by scripts/coverage.py
    "make.inc",  # the upstream make build, dropped
    "vrel.f90",  # the 4.0.1 file, named in known_issues.md
}
# the directories of the tree that hold files (not the builds)
TREE = ("src", "scripts", "tests", "doc", "cmake", ".github", ".vscode")
EXTENSIONS = "f90|F90|inc|py|sh|md|toml|dat|txt|cmake|yaml|yml|json|out"
NAME = re.compile(rf"`([^`\s]+\.(?:{EXTENSIONS}))`|`(\w+_m)`")
PLACEHOLDER = re.compile(r"[A-Z]{3,}|[*<|]")


def tree_files(root):
    """The paths of the tree's files, relative to root, and their names."""
    paths = {p.name for p in root.iterdir() if p.is_file()}
    for top in TREE:
        for path in (root / top).rglob("*"):
            if path.is_file():
                paths.add(str(path.relative_to(root)))
    return paths, {pathlib.Path(p).name for p in paths}


def modules(root):
    found = set()
    for path in [*(root / "src").rglob("*.[fF]90"), *(root / "tests").rglob("*.f90")]:
        text = path.read_text(errors="replace")
        found |= {m.lower() for m in re.findall(r"(?im)^\s*module\s+(\w+)\s*$", text)}
    return found


def main():
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else __file__).resolve()
    if root.is_file():
        root = root.parent.parent
    paths, names = tree_files(root)
    mods = modules(root)
    docs = [root / "README.md"] + sorted(
        p for p in (root / "doc").glob("*.md") if p.name not in HISTORICAL
    )
    missing = []
    for doc in docs:
        for number, line in enumerate(doc.read_text().splitlines(), 1):
            for match in NAME.finditer(line):
                name = match.group(1) or match.group(2)
                if PLACEHOLDER.search(name) or name in ALLOWED:
                    continue
                if match.group(1):
                    found = name in paths or pathlib.Path(name).name in names
                else:
                    found = name.lower() in mods
                if not found:
                    missing.append(f"{doc.relative_to(root)}:{number}: {name}")
    for item in missing:
        print(item)
    if missing:
        sys.exit(
            f"docrefs.py: {len(missing)} names of files or modules that do not exist"
        )
    print(f"docrefs.py: {len(docs)} docs name only files and modules that exist")


if __name__ == "__main__":
    main()
