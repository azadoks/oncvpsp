"""The input as TOML text: a dictionary in the layout of the TOML input
(doc/users_guide.md), its tables and arrays of tables, written for the
library's parser (toml-f), which checks it"""

import json
import math
import os
from collections.abc import Mapping


def _bare(key):
    return key.replace("_", "").replace("-", "").isalnum()


def _key(key):
    return key if _bare(key) else json.dumps(key)


def _value(value):
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_value(vv) for vv in value) + "]"
    raise TypeError(f"no TOML value for {type(value).__name__}: {value!r}")


def _is_tables(value):
    return (
        isinstance(value, (list, tuple))
        and len(value) > 0
        and all(isinstance(vv, Mapping) for vv in value)
    )


def _table(name, table, header, lines):
    """the table under its header ([name], [[name]] or none), its values
    then its subtables and arrays of tables; None values left out"""
    table = {kk: vv for kk, vv in table.items() if vv is not None}
    nested = [
        (kk, vv)
        for kk, vv in table.items()
        if isinstance(vv, Mapping) or _is_tables(vv)
    ]
    plain = [(kk, vv) for kk, vv in table.items() if kk not in dict(nested)]
    if header and (plain or not nested or header.startswith("[[")):
        lines.append(header)
    for kk, vv in plain:
        lines.append(f"{_key(kk)} = {_value(vv)}")
    if plain:
        lines.append("")
    for kk, vv in nested:
        sub = f"{name}.{_key(kk)}" if name else _key(kk)
        if isinstance(vv, Mapping):
            _table(sub, vv, f"[{sub}]", lines)
        else:
            for item in vv:
                _table(sub, item, f"[[{sub}]]", lines)


def dumps(data):
    """the TOML text of a dictionary of tables, its None values left out
    (TOML has no null), as if not given"""
    lines = []
    _table("", data, None, lines)
    return "\n".join(lines) + "\n"


def merged(base, tables):
    """base with each table of tables merged over its own, key by key"""
    out = {kk: dict(vv) if isinstance(vv, Mapping) else vv for kk, vv in base.items()}
    for kk, vv in tables.items():
        if isinstance(vv, Mapping) and isinstance(out.get(kk), Mapping):
            out[kk].update(vv)
        else:
            out[kk] = vv
    return out


def input_text(inp=None, tables=None):
    """the TOML text, as bytes, of the input inp (a dictionary in the TOML
    layout, the TOML text, or the path of a TOML file: a path-like, or a
    string of one line ending in .toml) with the tables merged over it"""
    tables = tables or {}
    if isinstance(inp, str) and "\n" not in inp and inp.endswith(".toml"):
        inp = os.fspath(os.path.abspath(inp))
        inp = open(inp, "rb").read()
    if isinstance(inp, os.PathLike):
        with open(inp, "rb") as ff:
            inp = ff.read()
    if isinstance(inp, (bytes, str)) and not tables:
        return inp.encode() if isinstance(inp, str) else inp
    if isinstance(inp, (bytes, str)):
        import tomllib

        inp = tomllib.loads(inp.decode() if isinstance(inp, bytes) else inp)
    return dumps(merged(inp or {}, tables)).encode()
