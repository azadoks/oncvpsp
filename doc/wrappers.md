# Python and Julia

ONCVPSP runs as a library, `liboncvpsp`, with a C interface
(`include/oncvpsp.h`). The Python package `oncvpsp` (`python/oncvpsp`) and
the Julia package `ONCVPSP` (`julia/ONCVPSP`) call it: a run, or one stage
of a run, in the calling process, its results as arrays, with no files
written or parsed. The results are those of the executables bit for bit;
the text output, if kept, is theirs but for the banner.

## Building

The library is built with the executables (`ONCVPSP_WITH_C_API`, on by
default), as `build/lib/liboncvpsp.so`. Each package finds it through
`ONCVPSP_LIBRARY`, else beside its modules (a wheel's), else on the
system's library path:

```sh
export ONCVPSP_LIBRARY=$PWD/build/lib/liboncvpsp.so
export PYTHONPATH=$PWD/python                     # Python
julia --project=julia/ONCVPSP                     # Julia
```

`pip install .` builds the library and the Python package together
(scikit-build-core, `pyproject.toml`), the library inside the package.
[Building and Installing](building.md) has the build options, the
installation of each package and an example of each.

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
    # only qcut changes; the table's other keys are base's
    run = oncvpsp.run(base, irel=1, pseudopotentials={"qcut": [qcut, 8.0, 8.0]})
```

`irel` is 0 (non-), 1 (scalar-) or 2 (fully relativistic), or `"nr"`,
`"sr"`, `"fr"` (Julia: also `:nr`, `:sr`, `:fr`; Python: `oncvpsp.Irel`),
as the executables `oncvpspnr.x`, `oncvpsp.x` and `oncvpspr.x`.

## Runs and stages

| Python | Julia | Returns | What |
|---|---|---|---|
| `run(inp, irel)` | `oncvpsp(inp; irel)` | `Run` | the whole run |
| `all_electron(inp, irel)` | `all_electron(inp; irel)` | `AllElectron` | the all-electron atom (fields `ae.*`) |
| `run(inp, ae=ae)` | `oncvpsp(inp; ae)` | `Run` | a run from that atom, the same as the whole run |
| `channel(inp, ae, l, mode)` | `channel(inp, ae, l; mode)` | `Channel` | one channel (fields `channel.*`) |
| `has_node(inp, ae, l)` | `has_node(inp, ae, l)` | a Boolean | whether its first pseudo-wave has a node |
| `coupled(run, icmod, fcfact=, rcfact=)` | `coupled(run, icmod; fcfact, rcfact)` | `Run` | the model core and the tests again with another model core |
| `confine(run, ...)` | `confine(run; ...)` | `Confinement` | the confined atom, as `oncvpsp_confine.x` (fields `confine.*`) |

In Julia, `ONCVPSP.run` is `oncvpsp` by the C function's name (unexported,
as Base has a `run`), and `Channel` is `ONCVPSP.Channel` (Base has a
`Channel` too). `confine`'s `radii` and `walls` are `None` (`nothing`) for
their defaults and empty for none.

A channel's `mode` is `"channel"` (the pseudo-waves, semi-local potential
and convergence curves), `"projectors"` (also the local potential and the
channel's projectors) or `"node"`, which stops at the node test of the
first projector, the cheapest way to scan for a node. The input of a
channel or a run from an atom must be of that atom: the same nucleus,
configuration, functional, mesh and speed of light.

### Runs assembled from channels

A channel depends on its own inputs and the all-electron atom alone, so a
search that changes one channel at a time can keep channels and assemble
runs from them:

| Python | Julia | What |
|---|---|---|
| `local_potential(inp, ae)` | `local_potential(inp, ae)` | the local potential alone, a `Local` (fields `local.*`) |
| `channel(inp, ae, l, "projectors", local=loc, tests=True)` | `channel(inp, ae, l; mode=:projectors, loc, tests=true)` | a channel on that local potential, with its own tests |
| `assemble(inp, ae, channels, local=loc)` | `assemble(inp, ae, channels; loc)` | the run from those channels |

`assemble` takes one channel of each l from 0 to lmax, in any order, and
for lloc < 4 the local channel from the `Local` (or a channel of l =
lloc). Each channel, and the local potential, must be of the input and the
atom given, checked by what each was made from; a run assembled from them
is then the whole run's bit for bit (its text output has none of the
channels' lines). In Julia the keyword is `loc`, `local` being a keyword
of the language.

A channel's tests are its log derivatives (`channel.logder.*`, on the
input's window or `energies=(emin, emax, de)`) and its ghost tests
(`channel.ghosts.*`). A whole run makes them on the local potential after
switching it, beyond the largest core radius rcmax, to the ionic tail of
the pseudo-atom's own density, which depends on every channel and the
model core; a channel alone has the potential before that switch. Its log
derivatives are therefore the run's bit for bit when they read the
potential no further than rcmax (Python `ch.exact`, Julia `ONCVPSP.isexact(ch)`, the field
`channel.logder.exact`), and otherwise differ from them at round-off in
the test inputs. Its ghost tests reach beyond rcmax always: in the test
inputs their counts are the run's and their energies within 1e-6 Ha.

## Results

Each function returns a handle, a mapping (Python's `Mapping`, Julia's
`AbstractDict`) of the fields it holds to arrays. `handle[name]` copies a
field into a new array, in the Fortran shape (NumPy arrays in Fortran
order, Julia arrays as they are), and `keys(handle)` lists them; the names
are the Fortran ones: `pot.*` the pseudopotential (`src/pseudo_m.f90`),
`tests.*` its tests (`src/test_results_m.f90`), `core.*` the model core,
the rest of a run's results (`src/results_m.f90`). `field_info(name)`
(`fieldinfo(name)`) says what a field is: its axes, units and a sentence,
from the library (`oncv_field_info`); `oncvpsp.FIELDS` (`ONCVPSP.fields`)
has them all.

```python
>>> oncvpsp.field_info("tests.convergence.cross")
tests.convergence.cross ([ecut,err], level, iproj, l, ikap) [Ha]: (Ecut, err)
where the residual kinetic energy per electron crosses each level, ...
```

An axis is `r` the radial mesh, `iproj` the projector, `l`, `ikap` (1, or
2 for j=l-1/2 when fully relativistic), and so on (`include/oncvpsp.h`);
an axis of different quantities is their names in brackets, `[ecut,err]`.
The densities of the fields are 4π times the density (∫ ρ r² dr is the
charge); the projectors `vkb` are r·β, normalized.

`handle.dims` gives the sizes, and `run.write("Si.upf")` (Julia:
`write("Si.upf", run)`) the psp8, UPF or PSML file the executables write.
A run also gives the physics without the indices:

| Python | Julia | What |
|---|---|---|
| `run.mesh` | `mesh(run)` | the radial mesh (bohr) |
| `run.local_potential` | `local_potential(run)` | the local potential (Ha) |
| `run.semilocal(l)` | `semilocal(run, l)` | the unscreened semi-local potential of `l` |
| `run.projectors(l, ikap)` | `projectors(run, l; ikap)` | `(energies, vkb)`, the `nproj` projectors of `l` |
| `run.kernel(l, ikap)` | `kernel(run, l; ikap)` | the nonlocal operator Σᵢ vkbᵢ(r) eᵢ vkbᵢ(r′), (mmax, mmax) |
| `run.valence_density` | `valence_density(run)` | the pseudo valence density (e/bohr³) |
| `run.core_density`, `run.model_core` | `core_density(run)`, `model_core(run)` | the core density and its model |
| `run.convergence(l, ikap, iproj)` | `convergence(run, l; ikap, iproj)` | rows (Ecut, err) of the pseudo-wave's convergence |

A channel of mode `"projectors"` has `projectors` and `kernel` too. The
kernel, not the projectors, is the thing to compare between two
pseudopotentials: the projectors are one basis of it among many.

```python
with oncvpsp.run("Si.toml") as run:   # freed at the end of the block
    rr = run.mesh
    print(run.convergence(0))         # l=0, its first pseudo-wave
    kk = run.kernel(1)
```

`oncvpsp.build_info()` (`build_info()`) gives the library's build: its
git commit, the state of the source against it and the build time, the
executables' `build:` line, by which results kept can be tied to the
ONCVPSP that made them.

`text=True`, the default, keeps the text output (`handle.text`), `False`
discards it, and a path writes it to that file.

## Failures

A failure throws `OncvpspError`: its `ifail`, its kind and stage, the
channel (`ll`) and the message's number (`msg`, a `msg_*` of
`src/failure_m.f90`, kept from one version to the next), the message
itself from the text output. In Python the kind and stage are
`oncvpsp.FailureKind` and `oncvpsp.Stage`, and the error is of the class
of its kind: `InputError`, `ConvergenceError`, `NotFoundError`,
`ConstructionError` (a node in a pseudo-wave as its `NodeError`),
`NumericalError`, `FileIOError`, `InternalError`. In Julia they are
Symbols, `err.kind === :convergence`, `err.stage === :ae`, and
`ONCVPSP.isnode(err)` tells a node.

```python
try:
    run = oncvpsp.run(base, pseudopotentials={"qcut": [9.0, 8.0, 8.0]})
except oncvpsp.NodeError:
    ...
except oncvpsp.OncvpspError as err:
    print(err.kind, err.stage, err.ll)
```

With `check=False` (`check=false`) it is left in the handle instead:
`handle.ifail` and `handle.failure`. Nothing stops the process; a handle
is freed when it is no longer referenced, by `close`, at the end of a
`with` block (Python) or of a `do` block (Julia's `oncvpsp(inp) do run
... end`, and `all_electron`'s).

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
