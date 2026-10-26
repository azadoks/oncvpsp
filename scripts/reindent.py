#!/usr/bin/env python3
"""Re-indent free-form Fortran in Hamann's layout.

  - unit statements (program, module, subroutine, function) and the
    statements directly inside them in column 2;
  - each block body 2 columns right of its opener; else/case/contains-like
    lines at the opener's column;
  - full-line comments in column 1;
  - continuation lines start with '&' in column 1; the text after it keeps
    its position relative to the statement's first line;
  - preprocessor lines unchanged.

The object files compile byte-identical before and after.  Run by
pre-commit on the Fortran files (see doc/coding_standards.md).

Usage: reindent.py [--check] file...
       reindent.py --stdin [name]

With --stdin the text comes from standard input and the re-indented text
goes to standard output, which is how an editor runs it as a formatter
(doc/developer_setup.md); name is the file's path, used in the warnings and
to leave the imported files alone.  Warnings go to standard error.
"""

import re
import sys

BASE = 1
STEP = 2

OPEN_RE = [
    re.compile(r"^(\w+\s*:\s*)?do\b(?!\s*=)"),
    re.compile(r"^(\w+\s*:\s*)?select\s*(case|type)\b"),
    re.compile(r"^(abstract\s+)?interface\b"),
    re.compile(r"^type\s*(,[^:]*)?::\s*\w+$"),
    re.compile(r"^type\s+\w+$"),
    re.compile(r"^(\w+\s*:\s*)?associate\s*\("),
    re.compile(r"^(\w+\s*:\s*)?block$"),
    re.compile(r"^enum\b"),
]
MID_RE = [
    re.compile(r"^else\b"),
    re.compile(r"^elseif\b"),
    re.compile(r"^elsewhere\b"),
    re.compile(r"^case\s*(\(|default\b)"),
    re.compile(r"^(class|type)\s+is\s*\("),
    re.compile(r"^class\s+default\b"),
]
CLOSE_RE = re.compile(
    r"^end\s*(do|if|select|where|forall|associate|block|interface|type|enum)\b"
)


def split_code(line, in_str):
    """Return (code, comment, in_str_at_end) for a physical line.

    in_str is the quote character of a character context continued from the
    previous line, or None.  Strings are replaced by 'S' in code.
    """
    code = []
    i = 0
    q = in_str
    n = len(line)
    if q is not None:
        # a continued character context resumes after a leading '&'
        s = line.lstrip()
        if s.startswith("&"):
            i = len(line) - len(s) + 1
    while i < n:
        ch = line[i]
        if q:
            if ch == q:
                if i + 1 < n and line[i + 1] == q:
                    i += 2
                    continue
                q = None
                code.append("S")
            elif ch == "&" and line[i + 1 :].strip() == "":
                # continued inside the string
                return "".join(code), "", q
            i += 1
            continue
        if ch in "'\"":
            q = ch
        elif ch == "!":
            return "".join(code), line[i:], None
        else:
            code.append(ch)
        i += 1
    return "".join(code), "", q


def balanced_end(s, start):
    """Index just after the parenthesis group that starts at s[start]."""
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "(":
            depth += 1
        elif s[i] == ")":
            depth -= 1
            if depth == 0:
                return i + 1
    return len(s)


def classify(stmt):
    """'open', 'mid', 'close', 'midopen' or None for one statement."""
    s = stmt.strip().lower()
    s = re.sub(r"\s+", " ", s)
    if not s:
        return None
    if CLOSE_RE.match(s) or re.match(
        r"^end(do|if|select|where|forall|interface|type)\b", s
    ):
        return "close"
    m = re.match(r"^(else\s*if)\s*\(", s)
    if m:
        return "mid"
    for r in MID_RE:
        if r.match(s):
            return "mid"
    m = re.match(r"^(\w+\s*:\s*)?if\s*\(", s)
    if m:
        j = balanced_end(s, s.index("(", m.end() - 1))
        return "open" if s[j:].strip() == "then" else None
    m = re.match(r"^(\w+\s*:\s*)?(where|forall)\s*\(", s)
    if m:
        j = balanced_end(s, s.index("(", m.end() - 1))
        return "open" if s[j:].strip() == "" else None
    for r in OPEN_RE:
        if r.match(s):
            return "open"
    return None


def reindent(text, fname):
    lines = text.split("\n")
    out = []
    depth = 0
    pp = []  # preprocessor stack: [depth at #if, depth at end of first branch]
    i = 0
    warn = []
    n = len(lines)
    while i < n:
        line = lines[i]
        st = line.strip()
        if st == "":
            out.append("")
            i += 1
            continue
        if st.startswith("#"):
            d = st[1:].strip()
            if re.match(r"if", d):
                pp.append([depth, None])
            elif re.match(r"(else|elif)", d):
                if pp[-1][1] is None:
                    pp[-1][1] = depth
                depth = pp[-1][0]
            elif re.match(r"endif", d):
                top = pp.pop()
                if top[1] is not None and top[1] != depth:
                    warn.append(
                        f"{fname}:{i + 1}: branches end at depths {top[1]} and {depth}"
                    )
                    depth = max(depth, top[1])
            out.append(line.rstrip())
            i += 1
            continue
        if st.startswith("!"):
            out.append(st)
            i += 1
            continue
        # a statement: gather its physical lines
        first = i
        phys = [i]
        code, com, q = split_code(line, None)
        codes = [code]
        strq = [None]
        while code.rstrip().endswith("&") or q is not None:
            i += 1
            while i < n and (
                lines[i].strip() == ""
                or lines[i].strip().startswith("!")
                or lines[i].strip().startswith("#")
            ):
                phys.append(i)
                codes.append(None)
                strq.append(None)
                i += 1
            if i >= n:
                break
            phys.append(i)
            strq.append(q)
            code, com, q = split_code(lines[i], q)
            codes.append(code)
        i += 1
        joined = " ".join(
            c.strip().rstrip("&").lstrip("&") for c in codes if c is not None
        )
        # statements separated by ';'
        kinds = [classify(seg) for seg in joined.split(";")]
        kinds_set = [k for k in kinds if k]
        if len(kinds_set) > 1:
            warn.append(f"{fname}:{first + 1}: several block statements on one line")
        kind = kinds_set[0] if kinds_set else None
        if kind in ("close", "mid"):
            depth -= 1
        if depth < 0:
            warn.append(f"{fname}:{first + 1}: negative depth")
            depth = 0
        ind = BASE + STEP * depth
        orig = len(lines[first]) - len(lines[first].lstrip())
        delta = ind - orig
        if lines[first].lstrip()[:1].isdigit():
            warn.append(f"{fname}:{first + 1}: statement label")
        out.append(" " * ind + lines[first].strip())
        prev_code = codes[0]
        for k, j in enumerate(phys[1:], 1):
            pl = lines[j]
            ps = pl.strip()
            if ps == "":
                out.append("")
                continue
            if ps.startswith("#"):
                out.append(pl.rstrip())
                continue
            if ps.startswith("!"):
                out.append(ps)
                continue
            if strq[k] is not None:
                # continued character context: the text is part of the string
                if ps.startswith("&"):
                    out.append(ps)
                else:
                    warn.append(
                        f"{fname}:{j + 1}: string continued without a leading &, left as is"
                    )
                    out.append(pl.rstrip())
                prev_code = codes[k]
                continue
            lead = len(pl) - len(pl.lstrip())
            if ps.startswith("&"):
                rest = ps[1:]
                ws = len(rest) - len(rest.lstrip())
                body = rest.lstrip()
                pc = prev_code.rstrip()[:-1] if prev_code is not None else ""
                if (
                    ws == 0
                    and pc
                    and not pc[-1].isspace()
                    and (pc[-1].isalnum() or pc[-1] == "_")
                    and body[:1].isalnum()
                ):
                    # a token split across the lines
                    out.append("&" + rest)
                else:
                    # position of the text after '&', moved with the statement
                    col = lead + 1 + ws + delta
                    out.append("&" + " " * max(1, col - 1) + body)
            else:
                col = lead + delta
                out.append("&" + " " * max(1, col - 1) + ps)
            prev_code = codes[k]
        if kind in ("open", "mid"):
            depth += 1
    if depth != 0:
        warn.append(f"{fname}: depth {depth} at end of file")
    if pp:
        warn.append(f"{fname}: unclosed preprocessor block")
    return "\n".join(out), warn


# Imported code, kept as written upstream: the same list as in
# .pre-commit-config.yaml and fortitude.toml.  pre-commit never passes these
# files; an editor formatting on save does, and gets them back unchanged.
IMPORTED = (
    "src/exc_off_pbe_m.f",
    "src/functionals.F90",
    "src/m_libxc_list.F90",
    "src/gg1cc.f90",
    "src/gp1cc.f90",
    "src/gpp1cc.f90",
)


def from_stdin(name):
    text = sys.stdin.buffer.read().decode("latin-1")
    if name.replace("\\", "/").endswith(IMPORTED):
        new = text
    else:
        new, warn = reindent(text, name)
        for w in warn:
            print("warning:", w, file=sys.stderr)
    sys.stdout.buffer.write(new.encode("latin-1"))
    return 0


def main():
    args = sys.argv[1:]
    if args and args[0] == "--stdin":
        return from_stdin(args[1] if len(args) > 1 else "<stdin>")
    check = "--check" in args
    files = [a for a in args if a != "--check"]
    bad = 0
    for f in files:
        text = open(f, encoding="latin-1").read()
        new, warn = reindent(text, f)
        for w in warn:
            print("warning:", w, file=sys.stderr)
        if new != text:
            bad += 1
            if check:
                print("would re-indent", f)
            else:
                open(f, "w", encoding="latin-1").write(new)
    return 1 if (check and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
