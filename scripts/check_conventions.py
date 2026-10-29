#!/usr/bin/env python3
"""Checks the conventions of doc/coding_standards.md that neither fortitude
nor gfortran checks, in the Fortran of src/ (the imported files excluded).

usage: check_conventions.py [FILE]... [--all] [--update-baseline]
                            [--baseline FILE] [--list]

Each finding is a rule, a file, a procedure and a detail (no line number, so
that it survives edits elsewhere in the file).  The findings already in the
baseline (scripts/conventions_baseline.txt, one per line, counted) are
allowed: the check fails only on new ones, so the conventions hold for new
code while the old code is fixed, and --update-baseline rewrites the
baseline after a fix (it only shrinks, unless --list shows why not).  With
files, only those are checked (pre-commit); with --all or none, every
file of src/ (the ctest test).

The rules (doc/coding_standards.md and doc/consistency_audit.md):
  sections  a dummy argument is declared under the label of its intent
            (!Input variables for intent(in), ...), a local under !Local
            variables; a module procedure with arguments has the labels
  header    a module procedure has a !! description and one @param per dummy
            argument, in order (a,b and a..b name several)
  glossary  the @param text of a name in doc/glossary.md starts with the
            glossary's text for it
  message   an ERROR or WARNING message names the routine that writes it,
            its host, or a public procedure of its module; its format
            starts with a blank line ('(/...') and prints numbers only as
            i0, f12.6 or es12.4; no string in it is continued across lines
            (the continuation would copy the indentation into the text).
            The messages of 4.0.1 (scripts/messages_4.0.1.txt) are exempt
  stop      a module procedure stops only with error stop 1 (rule 5)
  ierr      an ierr returned by a call is used after it (or passed up as
            the procedure's own)
  moddoc    a module's variables and constants, and a type's components,
            have a !> or !< comment, one name each
  license   the file starts with the copyright and GPL notice
  onechar   no one-letter names (rule 5), type components included
  meshstep  the mesh step is not re-derived from rr (log(rr(...)) outside
            logmesh_m: lm_step)
  use       a module's use statements serve its own declarations (the
            procedures use what they call)
  format    no d edit descriptor (an exponent printed D): e or es
  kinds     no _dp literals (1.0d0) and no real64/real128 kinds outside
            constants_m
  names     module x_m is in x_m.f90; the public procedures of a module in
            the prefix table of doc/coding_standards.md start with its prefix
"""

import argparse
import collections
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASELINE = ROOT / "scripts" / "conventions_baseline.txt"
MESSAGES_401 = [
    line
    for line in (ROOT / "scripts" / "messages_4.0.1.txt").read_text().splitlines()
    if line and not line.startswith("#")
]


def read_glossary(path):
    """{name: text} of the table after the <!-- glossary --> marker."""
    terms, inside = {}, False
    for line in path.read_text().splitlines():
        if line.startswith("<!-- glossary"):
            inside = True
            continue
        if inside:
            if not line.startswith("|"):
                if terms:
                    break
                continue
            cells = [c.strip() for c in line.strip("|").split("|")]
            name = cells[0].strip("`")
            if name and name != "name" and not set(name) <= {"-"}:
                terms[name] = cells[1].split(" / ")
    return terms


GLOSSARY = read_glossary(ROOT / "doc" / "glossary.md")


def read_prefixes(path):
    """{module: prefix} of the table after the <!-- prefixes --> marker."""
    table, inside = {}, False
    for line in path.read_text().splitlines():
        if line.startswith("<!-- prefixes"):
            inside = True
            continue
        if inside:
            if not line.startswith("|"):
                if table:
                    break
                continue
            cells = [c.strip().strip("`") for c in line.strip("|").split("|")]
            if cells[0] and cells[0] != "module" and not set(cells[0]) <= {"-"}:
                table[cells[0]] = cells[1]
    return table


PREFIXES = read_prefixes(ROOT / "doc" / "coding_standards.md")
# the numeric edit descriptors a message may use
MESSAGE_NUMBERS = {"i0", "f12.6", "es12.4"}
IMPORTED = {"exc_off_pbe_m.f", "m_libxc_list.F90"}
# the procedures that may stop (rule 5): the command line and oncvpsp_main
STOPS = {"parse_command_line", "oncvpsp_main"}

RE_UNIT = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:integer|real|logical|complex|character|type)"
    r"(?:\s*\([^)]*\))?\s+)?(?:(?:pure|elemental|recursive)\s+)*"
    r"(subroutine|function|program|module)\s+(\w+)\s*(?:\(([^)]*)\))?",
    re.I,
)
RE_END = re.compile(r"^\s*end\s*(subroutine|function|program|module)\b", re.I)
RE_DECL = re.compile(
    r"^(integer|real|logical|complex|character|type\s*\(|class\s*\(|procedure\s*\()"
    r"(.*?)::(.*)$",
    re.I,
)
LABELS = {
    "in": "Input variables",
    "inout": "Input/output variables",
    "out": "Output variables",
    None: "Local variables",
}
STRING = re.compile(r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"")


def strip_comment(line):
    """The code of a line, without its comment."""
    quote = None
    for ii, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return line[:ii]
    return line


def open_string_at_end(line):
    """Whether a line ends inside a string (a string continued with &)."""
    quote = None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return False
    return quote is not None


def statements(lines):
    """(first line number, joined code) of each statement."""
    ii = 0
    while ii < len(lines):
        code = strip_comment(lines[ii]).rstrip()
        start = ii
        if not code.strip():
            ii += 1
            continue
        parts = [code]
        while parts[-1].endswith("&") and ii + 1 < len(lines):
            parts[-1] = parts[-1][:-1]
            ii += 1
            nxt = strip_comment(lines[ii]).rstrip()
            while not nxt.strip() and ii + 1 < len(lines):
                ii += 1
                nxt = strip_comment(lines[ii]).rstrip()
            nxt = nxt.lstrip()
            parts.append(nxt[1:] if nxt.startswith("&") else nxt)
        ii += 1
        yield start + 1, "".join(parts).strip()


def split_top(text):
    """text split at the commas outside parentheses."""
    out, depth, buf = [], 0, ""
    for ch in text:
        depth += ch in "(["
        depth -= ch in ")]"
        if ch == "," and depth == 0:
            out.append(buf)
            buf = ""
        else:
            buf += ch
    out.append(buf)
    return out


class Unit:
    """A module, program or procedure, with what the rules need of it."""

    def __init__(self, kind, name, args, line, host):
        self.kind, self.name, self.line, self.host = kind, name, line, host
        self.args = [a.strip().lower() for a in args.split(",") if a.strip()]
        self.decls = []  # (line, intent, names, attrs, label)
        self.body = []  # (line, statement)
        self.header = []  # the !! lines
        self.public = set()
        self.interface = []  # the statements of its interface blocks


def parse(path):
    """The units of a file (modules, procedures) and its lines."""
    lines = path.read_text().splitlines()
    units, stack = [], []
    interface = typedef = 0
    for ln, st in statements(lines):
        low = st.lower()
        if re.match(r"(abstract\s+)?interface\b", low):
            interface += 1
            continue
        if re.match(r"end\s*interface\b", low):
            interface -= 1
            continue
        if interface:
            if stack:
                stack[-1].interface.append(st)
            continue
        if re.match(r"end\s*type\b", low):
            typedef = 0
            continue
        if re.match(r"type\s*(,[^:]*)?(::)?\s*\w+\s*$", low) and not low.startswith(
            "type("
        ):
            typedef = ln
            if stack:
                stack[-1].body.append((ln, st))
            continue
        mm = RE_UNIT.match(st)
        if mm and not low.startswith(("end", "module procedure")):
            unit = Unit(
                mm.group(1).lower(),
                mm.group(2).lower(),
                mm.group(3) or "",
                ln,
                stack[:],
            )
            rm = re.search(r"\bresult\s*\(\s*(\w+)", low)
            unit.result = rm.group(1) if rm else unit.name
            jj = ln
            while strip_comment(lines[jj - 1]).rstrip().endswith("&"):
                jj += 1
            while jj < len(lines) and not lines[jj].strip():
                jj += 1
            while jj < len(lines) and lines[jj].lstrip().startswith("!!"):
                unit.header.append(lines[jj].strip())
                jj += 1
            units.append(unit)
            stack.append(unit)
            continue
        if RE_END.match(low) or (low == "end" and stack):
            stack.pop()
            continue
        if not stack:
            continue
        cur = stack[-1]
        cur.body.append((ln, st))
        pm = re.match(r"public\s*::\s*(.*)$", low)
        if pm and cur.kind == "module":
            cur.public |= {n.strip() for n in pm.group(1).split(",")}
        dm = RE_DECL.match(st)
        if dm:
            attrs = (dm.group(1) + dm.group(2)).lower()
            it = re.search(r"intent\s*\(\s*(in\s*out|inout|in|out)\s*\)", attrs)
            intent = it.group(1).replace(" ", "") if it else None
            names = [
                m.group(1).lower()
                for m in (re.match(r"\s*(\w+)", x) for x in split_top(dm.group(3)))
                if m
            ]
            label = None
            for kk in range(ln - 2, cur.line - 1, -1):
                lm = re.match(
                    r"!\s*(Input/output|Input|Output|Local) variables",
                    lines[kk].strip(),
                )
                if lm:
                    label = lm.group(1) + " variables"
                    break
            cur.decls.append((ln, intent, names, attrs, label, bool(typedef)))
    return units, lines


def param_names(header):
    """The dummy names of the @param lines, in order."""
    names = []
    for line in header:
        pm = re.match(r"!!\s*@param\s+(\S+)", line)
        if not pm:
            continue
        spec = re.sub(r"\(.*", "", pm.group(1)).lower()
        if ".." in spec:
            names.append(("range",) + tuple(spec.split("..")))
        else:
            names.extend(n for n in spec.split(",") if n)
    return names


def param_texts(header):
    """(name, text) of each @param line naming one dummy, continuation joined."""
    out = []
    for ii, line in enumerate(header):
        pm = re.match(r"!!\s*@param\s+(\w+)\s+(.*)$", line)
        if not pm:
            continue
        text = pm.group(2).strip()
        for nxt in header[ii + 1 :]:
            if re.match(r"!!\s*@param\b", nxt) or not nxt.strip("! "):
                break
            text += " " + nxt.lstrip("!").strip()
        out.append((pm.group(1).lower(), text))
    return out


def header_matches(params, args):
    """Whether the @param names cover args in order."""
    ii = 0
    for item in params:
        if isinstance(item, tuple):
            first, last = item[1], item[2]
            if ii >= len(args) or args[ii] != first or last not in args[ii:]:
                return False
            ii = args.index(last, ii) + 1
        else:
            if ii >= len(args) or args[ii] != item:
                return False
            ii += 1
    return ii == len(args)


def check_file(path, rel):
    """The findings of one file: (rule, file, procedure, detail)."""
    found = []
    units, lines = parse(path)

    def add(rule, unit, detail, line=0):
        found.append(((rule, rel, unit.name if unit else "-", detail), line))

    if not (lines and lines[0].strip() == "!" and "Copyright" in "".join(lines[1:4])):
        if not path.name.startswith("oncvpsp") or len(lines) > 10:
            add("license", None, "no copyright and GPL notice at the top")
    modules = [u for u in units if u.kind == "module"]
    module = modules[0] if modules else None
    # names
    # a stub stands in for the module when its library is off
    stem = path.stem.lower().replace("_stub_m", "_m")
    if module and module.name != stem:
        add("names", module, f"module {module.name} in {path.name}")
    if module and module.name in PREFIXES:
        prefix = PREFIXES[module.name]
        for unit in units:
            if (
                unit.kind in ("subroutine", "function")
                and unit.host == [module]
                and unit.name in module.public
                and not unit.name.startswith(prefix)
                and unit.name != prefix.rstrip("_")
            ):
                add("names", unit, f"public {unit.name} without the prefix {prefix}")
    for unit in units:
        procedure = unit.kind in ("subroutine", "function")
        in_module = any(h.kind == "module" for h in unit.host)
        internal = any(
            h.kind in ("subroutine", "function", "program") for h in unit.host
        )
        # sections
        if procedure:
            for ln, intent, names, attrs, label, _ in unit.decls:
                if "parameter" in attrs or label is None:
                    continue
                want = LABELS[intent]
                if intent is None and unit.kind == "function" and unit.result in names:
                    want = LABELS["out"]
                elif intent is None and set(names) & set(unit.args):
                    want = LABELS["in"]  # dummy procedures
                if label != want:
                    add(
                        "sections",
                        unit,
                        f"{','.join(names)} ({intent or 'local'}) under !{label}",
                        ln,
                    )
        if (
            procedure
            and unit.args
            and any(h.kind == "module" for h in unit.host)
            and not any(h.kind != "module" for h in unit.host)
            and all(
                label is None
                for _, _, _, attrs, label, _ in unit.decls
                if "parameter" not in attrs
            )
        ):
            add("sections", unit, "no section labels", unit.line)
        if procedure:
            for ln, st in unit.body:
                if re.match(r"use\b", st.lower()):
                    for kk in range(ln - 2, unit.line - 1, -1):
                        if re.match(
                            r"!\s*(Input/output|Input|Output|Local) variables",
                            lines[kk].strip(),
                        ):
                            add(
                                "sections",
                                unit,
                                f"{st.split(',')[0]} under a section label",
                                ln,
                            )
                            break
        # header
        if procedure and in_module and not internal:
            if not any(
                line.startswith("!! ") and "@param" not in line for line in unit.header
            ):
                add("header", unit, "no !! description", unit.line)
            params = param_names(unit.header)
            if unit.args and not params:
                add("header", unit, "no @param", unit.line)
            for name, text in param_texts(unit.header):
                want = GLOSSARY.get(name)
                if want and not text.startswith(tuple(want)):
                    add("glossary", unit, f"@param {name}: {text[:40]}", unit.line)
            if unit.args and params and not header_matches(params, unit.args):
                add(
                    "header",
                    unit,
                    "@param names differ from the dummy arguments",
                    unit.line,
                )
        # one-letter names
        for ln, _, names, _, _, _ in unit.decls:
            for name in names:
                if len(name) == 1:
                    add("onechar", unit, f"name {name}", ln)
        # use at module level only for what the module's own declarations need
        if unit.kind == "module":
            spec = [st for _, st in unit.body if not re.match(r"use\b", st.lower())]
            spec_text = " ".join(spec + unit.interface).lower()
            for ln, st in unit.body:
                um = re.match(r"use\s+(\w+)\s*,\s*only\s*:(.*)$", st.lower())
                if um:
                    names = [n.split("=>")[-1].strip() for n in um.group(2).split(",")]
                    used = [
                        n
                        for n in names
                        if re.search(r"\b%s\b" % re.escape(n), spec_text)
                    ]
                    if not used:
                        add(
                            "use",
                            unit,
                            f"use {um.group(1)} at module level, not used there",
                            ln,
                        )
        # module documentation
        if unit.kind == "module":
            for ln, _, names, attrs, _, _ in unit.decls:
                text = lines[ln - 1]
                above = lines[ln - 2].strip() if ln > 1 else ""
                documented = above.startswith("!>") or "!<" in text
                if not documented:
                    add("moddoc", unit, f"{','.join(names)} undocumented", ln)
                elif len(names) > 1:
                    add("moddoc", unit, f"{','.join(names)} documented together", ln)
        # statements
        hosts = {h.name for h in unit.host if h.kind != "module"}
        allowed = {unit.name} | hosts | (module.public if module else set())
        stmts = unit.body
        for kk, (ln, st) in enumerate(stmts):
            low = st.lower()
            if procedure and in_module and unit.name not in STOPS:
                sm = re.search(r"\b(error\s+stop|stop)\b\s*(.*)$", low)
                if sm and not re.match(r"(write|print)\b", low):
                    if sm.group(1) == "stop" or sm.group(2).strip() != "1":
                        add(
                            "stop",
                            unit,
                            f"{sm.group(1)} {sm.group(2).strip()}".strip(),
                            ln,
                        )
            if re.match(r"(if\s*\(.*?\)\s*)?write\s*\(", low):
                for fmt in re.findall(r"'(\([^']*\))'", st):
                    if re.search(r"(?<=[,(p])d\d+\.\d+", fmt.lower()):
                        add("format", unit, "d edit descriptor", ln)
            if re.match(r"(if\s*\(.*?\)\s*)?write\s*\(", low):
                joined = re.sub(r"(['\"])\s*//\s*(['\"])", "", st)
                text = " ".join(a or b for a, b in STRING.findall(joined))
                mm = re.search(r"(\w+):\s*(ERROR|WARNING)\b", text)
                if mm and any(m in text for m in MESSAGES_401):
                    mm = None
                if mm:
                    prefix = mm.group(1).lower()
                    if prefix not in allowed:
                        add("message", unit, f"prefix {prefix}: {mm.group(2)}", ln)
                    fm = re.match(
                        r"(?:if\s*\(.*?\)\s*)?write\s*\([^,]*,\s*'([^']*)'", st, re.I
                    )
                    if fm and not fm.group(1).startswith("(/"):
                        add(
                            "message",
                            unit,
                            f"no blank line before {prefix}: {mm.group(2)}",
                            ln,
                        )
                    if fm:
                        for desc in re.findall(
                            r"(?<![a-z])(1p|[a-z]+\d+(?:\.\d+)?)", fm.group(1).lower()
                        ):
                            if desc != "a" and desc not in MESSAGE_NUMBERS:
                                add("message", unit, f"edit descriptor {desc}", ln)
                    jj = ln - 1
                    while jj < len(lines) and strip_comment(
                        lines[jj]
                    ).rstrip().endswith("&"):
                        if open_string_at_end(lines[jj].rstrip()[:-1]):
                            add("message", unit, "string continued across lines", ln)
                        jj += 1
            cm = re.match(r"(?:if\s*\(.*\)\s*)?call\s+(\w+)\s*\((.*)\)\s*$", st, re.I)
            if cm:
                args = [a.strip().lower() for a in split_top(cm.group(2))]
                for arg in args:
                    if re.fullmatch(r"ierr\w*", arg) and arg not in unit.args:
                        after = " ".join(s for _, s in stmts[kk + 1 :]).lower()
                        if not re.search(r"\b%s\b" % arg, after):
                            add(
                                "ierr",
                                unit,
                                f"{arg} of {cm.group(1).lower()} unused",
                                ln,
                            )
            if path.name != "logmesh_m.f90" and re.search(
                r"\blog\s*\(\s*rr\s*\([^()]*\)\s*/\s*rr\s*\(", low
            ):
                add("meshstep", unit, "log(rr(...)) re-derives the mesh step", ln)
            if path.name != "constants_m.f90":
                if re.search(r"\d_dp\b", low):
                    add("kinds", unit, "_dp literal", ln)
                if re.search(r"\breal(64|128)\b", low):
                    add("kinds", unit, "real64/real128 outside constants_m", ln)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--all", action="store_true", help="every file of src/")
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    parser.add_argument("--update-baseline", action="store_true")
    parser.add_argument(
        "--list", action="store_true", help="every finding, baseline or not"
    )
    args = parser.parse_args()

    if args.files and not args.all:
        paths = [p.resolve() for p in args.files]
    else:
        paths = sorted((ROOT / "src").glob("*.[fF]90"))
    paths = [
        p
        for p in paths
        if p.name not in IMPORTED
        and p.suffix in (".f90", ".F90")
        and p.parent == ROOT / "src"
    ]
    found = []
    for path in paths:
        found += check_file(path, path.relative_to(ROOT).as_posix())
    keys = ["\t".join(f) for f, _ in found]
    lines = {"\t".join(f): line for f, line in found}

    if args.update_baseline:
        args.baseline.write_text(
            "# scripts/check_conventions.py: the findings allowed, one per line\n"
            + "".join(k + "\n" for k in sorted(keys))
        )
        print(f"{len(keys)} findings written to {args.baseline.name}")
        return 0

    checked = {p.relative_to(ROOT).as_posix() for p in paths}
    base = collections.Counter(
        line
        for line in (
            args.baseline.read_text().splitlines() if args.baseline.exists() else []
        )
        if line and not line.startswith("#") and line.split("\t")[1] in checked
    )
    new = collections.Counter(keys) - base
    fixed = base - collections.Counter(keys)

    def show(key):
        rule, rel, proc, detail = key.split("\t")
        return f"{rel}:{lines.get(key, 0)}: {rule}: {proc}: {detail}"

    if args.list:
        for key in sorted(keys):
            print(show(key))
    for key in sorted(new.elements()):
        print("new: " + show(key))
    if fixed:
        print(
            f"{sum(fixed.values())} findings of the baseline are fixed: "
            "run with --update-baseline"
        )
    if new:
        print(
            f"{sum(new.values())} findings not in the baseline "
            "(doc/coding_standards.md; scripts/check_conventions.py says each rule)"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
