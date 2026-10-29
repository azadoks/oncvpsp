# Python and Julia

ONCVPSP runs as a library, `liboncvpsp`, with a C interface
(`include/oncvpsp.h`). The Python package `oncvpsp` (`python/oncvpsp`) and
the Julia package `ONCVPSP` (`julia/ONCVPSP`) call it: a run, or one stage
of a run, in the calling process, its results as arrays, with no files
written or parsed. The results are those of the executables bit for bit;
the text output, if kept, is theirs but for the banner.

## Building

The library is built with the executables (`ONCVPSP_WITH_C_API`, on by
default; see [Developer Setup](developer_setup.md)), as
`build/lib/liboncvpsp.so`. Each package finds it through
`ONCVPSP_LIBRARY`, else beside its modules (a wheel's), else on the
system's library path:

```sh
export ONCVPSP_LIBRARY=$PWD/build/lib/liboncvpsp.so
export PYTHONPATH=$PWD/python                     # Python
julia --project=julia/ONCVPSP                     # Julia
```

`pip install .` builds the library and the Python package together
(scikit-build-core, `pyproject.toml`), the library inside the package.

## The input

The input is the TOML input of the executables (`--toml-input`, the
[Users Guide](users_guide.md)): the path of a TOML file, its text, or a
dictionary in its layout. Keyword arguments are tables merged over it, key
by key, so that one input serves a scan. `None` (`nothing`) values are
left out, as if not given, TOML having no null:

```python
import oncvpsp

base = oncvpsp.read_toml("Si.toml")
for qcut in (6.0, 7.0, 8.0):
    pp = dict(base["pseudopotentials"], qcut=[qcut, 8.0, 8.0])
    run = oncvpsp.run(base, irel=1, pseudopotentials=pp)
```

`irel` is 0 (non-), 1 (scalar-) or 2 (fully relativistic), or `"nr"`,
`"sr"`, `"fr"`, as the executables `oncvpspnr.x`, `oncvpsp.x` and
`oncvpspr.x`.

## Runs and stages

| Python | Julia | What |
|---|---|---|
| `run(inp, irel)` | `ONCVPSP.run(inp; irel)` | the whole run |
| `all_electron(inp, irel)` | `all_electron(inp; irel)` | the all-electron atom (fields `ae.*`) |
| `run(inp, ae=ae)` | `ONCVPSP.run(inp; ae)` | a run from that atom, the same as the whole run |
| `channel(inp, ae, l, mode)` | `channel(inp, ae, l; mode)` | one channel (fields `channel.*`) |
| `has_node(inp, ae, l)` | `has_node(inp, ae, l)` | whether its first pseudo-wave has a node |
| `coupled(run, icmod, fcfact, rcfact)` | `coupled(run, icmod, fcfact, rcfact)` | the model core and the tests again with another model core |
| `confine(run, ...)` | `confine(run; ...)` | the confined atom, as `oncvpsp_confine.x` (fields `confine.*`) |

A channel's `mode` is `"channel"` (the pseudo-waves, semi-local potential
and convergence curves), `"projectors"` (also the local potential and the
channel's projectors) or `"node"`, which stops at the node test of the
first projector, the cheapest way to scan for a node. The input of a
channel or a run from an atom must be of that atom: the same nucleus,
configuration, functional, mesh and speed of light.

## Results

Each function returns a handle. `handle[name]` copies a field into a new
array, in the Fortran shape (NumPy arrays in Fortran order, Julia arrays
as they are); `keys(handle)` lists the fields it holds, `oncvpsp.FIELDS`
(`ONCVPSP.fields`) all of them, by their Fortran names: `pot.*` the
pseudopotential (`src/pseudo_m.f90`), `tests.*` its tests
(`src/test_results_m.f90`), `core.*` the model core, the rest of a run's
results (`src/results_m.f90`). `handle.dims` gives the sizes, and
`handle.write("Si.upf")` (`ONCVPSP.write(handle, "Si.upf")`) the psp8, UPF
or PSML file the executables write.

```python
run = oncvpsp.run("Si.toml")
rr, vkb = run["pot.rr"], run["pot.vkb"]   # vkb[:, iproj, l, ikap]
print(run["tests.convergence.cross"][..., 0])
```

`oncvpsp.build_info()` (`build_info()`) gives the library's build: its
git commit, the state of the source against it and the build time, the
executables' `build:` line, by which results kept can be tied to the
ONCVPSP that made them.

`text=True`, the default, keeps the text output (`handle.text`), `False`
discards it, and a path writes it to that file.

## Failures

A failure throws `OncvpspError`: its `ifail`, its kind and stage
(`kind_name`, `stage_name`), the channel (`ll`) and the message's number
(`msg`, a `msg_*` of `src/failure_m.f90`, kept from one version to the
next), the message itself from the text output. With `check=False`
(`check=false`) it is left in the handle instead: `handle.ifail` and
`handle.failure`. Nothing stops the process; a handle is freed when it is
no longer referenced, or by `close`.

## Notes

- A Python or Julia of its own toolchain (conda's) may load its own
  libquadmath, whose quad precision rounds otherwise than the compiler's:
  the results then differ from the executables' at 1e-9. Preloading the
  compiler's (`LD_PRELOAD=$(gfortran -print-file-name=libquadmath.so.0)`)
  restores them bit for bit, as the tests do.
- The library is loaded global (`RTLD_GLOBAL`), as MKL needs when it is
  linked as its separate libraries; with `mkl_rt` (CMake's
  `BLA_VENDOR=Intel10_64_dyn`) it does not.
- Handles share nothing: runs in several threads at once are independent,
  but one handle must not be used by two threads at once.
