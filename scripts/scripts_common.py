"""Helpers shared by the scripts, with only the standard library: numbers as
Fortran prints them, and the TOML input read with tomllib, changed and
written back.
"""

import re
import tomllib


def fortran_float(text):
    """The number of TEXT, with a Fortran D exponent (1.0D-02) or an E one."""
    return float(text.replace("D", "E").replace("d", "e"))


# TOML input: read with tomllib, one value changed, written back


def parse_key(key):
    """('table', 'name', index, '*' or None) of KEY, e.g. pseudopotentials.rc[0]
    or pseudopotentials.rc[*] (every element)."""
    match = re.fullmatch(r"(\w+)\.(\w+)(?:\[(\d+|\*)\])?", key)
    if not match:
        raise ValueError(
            f"KEY must be table.name, table.name[i] or table.name[*], not {key!r}"
        )
    table, name, index = match.groups()
    return table, name, index if index in (None, "*") else int(index)


def get_value(data, key):
    table, name, index = parse_key(key)
    value = data[table][name]
    return value if index in (None, "*") else value[index]


def set_value(data, key, value):
    table, name, index = parse_key(key)
    if index in (None, "*"):
        # an optional table the input leaves out, e.g. [log_mesh]
        data.setdefault(table, {})[name] = value
    else:
        data[table][name][index] = value


def apply_sets(data, sets):
    """data with each KEY=VALUE of sets applied, VALUE read as TOML."""
    for item in sets or []:
        key, _, text = item.partition("=")
        set_value(data, key.strip(), tomllib.loads(f"v = {text}")["v"])
    return data


def toml_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(value, list):
        if not value:
            return "[]"
        return "[ " + ", ".join(toml_value(v) for v in value) + ",]"
    raise TypeError(f"no TOML form for {value!r}")


def dump_toml(data):
    """TOML text for the input: its top-level keys, then its arrays of tables
    and its tables (of scalars and arrays), with floats written so that they
    read back exactly; a value None, which TOML cannot hold, is left out."""

    def pairs(table):
        return [f"{k} = {toml_value(v)}" for k, v in table.items() if v is not None]

    def tables(value):
        return isinstance(value, list) and value and isinstance(value[0], dict)

    top = {k: v for k, v in data.items() if not isinstance(v, dict) and not tables(v)}
    out = [*pairs(top), ""] if top else []
    for name, value in data.items():
        if tables(value):
            for entry in value:
                out += [f"[[{name}]]", *pairs(entry), ""]
    for name, value in data.items():
        if isinstance(value, dict):
            out += [f"[{name}]", *pairs(value), ""]
    return "\n".join(out)
