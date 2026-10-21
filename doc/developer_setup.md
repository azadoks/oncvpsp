# Developer setup

What a clone needs for working on the code: the tools, the commit hooks,
the builds, the tests and an editor that keeps to the layout. The rules
themselves are in [Coding Standards](coding_standards.md).

## Branches

Development happens on `oncvpsp-opt`. `master` mirrors upstream
([oncvpsp/oncvpsp](https://github.com/oncvpsp/oncvpsp)) and takes no
commits of its own, so check the branch before starting work.

## Tools

The hooks, the linters and the editor integration are Python tools. CI
uses Python 3.12. Install them once, outside any project environment, with
[pipx](https://pipx.pypa.io) (or `pip install --user`):

```bash
pipx install pre-commit
pipx install fortitude-lint==0.9.2   # the version in .pre-commit-config.yaml
pipx install fortls
pipx install ruff
```

pre-commit keeps its own copies of the hooks' tools, at the versions pinned
in `.pre-commit-config.yaml`; the installs above are for the command line
and the editor, and should match those versions.

## Commit hooks

```bash
pre-commit install          # once per clone
pre-commit run --all-files  # what CI does on every push
```

The hooks (listed in [Coding Standards](coding_standards.md#linting)) trim
whitespace, re-indent the Fortran with `scripts/reindent.py`, and run
fortitude, shellcheck and ruff. fortitude only reports; to apply the fixes
it knows, and re-indent afterwards (fortitude's fixes can change the
layout):

```bash
fortitude check --fix src/file.f90
python3 scripts/reindent.py src/file.f90
```

`reindent.py --check file...` reports without changing anything, and
`reindent.py --stdin` filters standard input to standard output, which is
how an editor uses it (below).

## Builds

CMake is the build for day-to-day work, with the presets in
`CMakePresets.json`:

```bash
cmake --preset default            # RelWithDebInfo in build/, the reference build
cmake --build --preset default --parallel
ctest --preset default
```

The `debug` preset (in `build-debug/`) adds run-time checks and
floating-point traps, and is the one to use when something goes wrong.
The `study` preset (in `build-study/`) is the `release` build (`-O3
-ffast-math -funroll-loops`, the flags of the upstream `make.inc`) with
`-march=native`, kept in its own directory for the oncvpsp-opt pseudopotential
study. It is about 35 % faster than the `default` build (29_Cu and 83_Bi,
scalar-relativistic) and passes the same tests, as it moves the nonlocal
kernels by about 1e-8, but only runs on CPUs like the one that built it.

`ONCVPSP_WITH_FR=OFF` leaves out the fully-relativistic executables and their
tests.

`ONCVPSP_BLAS_DOT=ON` computes the overlaps of the pseudo-atom's outward
integration inside the first core radius (`outwf_overlaps`, most of their
work) by one `dgemm` of the linked BLAS instead of a compensated sum in the
code. It is faster with an optimized BLAS, but the results then depend on
the BLAS at round-off (the references pass either way), so it is off by
default.

`CMakeUserPresets.json` (ignored by git) is the place for personal presets
such as compiler paths or `Libxc_ROOT`.

CMake is the only build: the upstream make build (`make.inc`, `src/Makefile`
and `tests/data/TEST.sh`) was dropped, as ctest covers its tests. A source
file is added or removed in `src/CMakeLists.txt`.

After changing any routine's argument list, run
`scripts/check_interfaces.sh` (it needs a CMake build in `build/`), which
checks every call against the definition of the routine it calls.

## Tests

`ctest --preset default` runs everything; `-L <label>` selects a part of it
(`unit`, `sr`, `r`, `nr`, `split`, `mesh`, `toml`, `psml`, `pipe`,
`errors`, `echo`, `psp8`, `nist`) and `-R <regex>` a test by name. Each test
leaves its outputs in its own directory under `build/tests`.

### The Python of the comparisons

The tests that compare outputs (`scripts/compare_output.py`, below) need a
Python with [pyoncvpsp](https://github.com/azadoks/pyoncvpsp), which parses
the outputs and the psp8 and PSML files in them, and
[upf-tools](https://pypi.org/project/upf-tools/), which parses the UPF
files. The configuration finds it through three cache variables:

- `ONCVPSP_TEST_PYTHON`: a Python that has both, e.g. that of a pyoncvpsp
  checkout (`uv sync` there, then
  `-DONCVPSP_TEST_PYTHON=/path/to/pyoncvpsp/.venv/bin/python`, after adding
  upf-tools to it with `uv pip install upf-tools`).
- Otherwise the configuration makes a virtual environment in the build tree,
  `build/test-venv`, with `uv` if it is on `PATH` (else with `python3 -m
  venv` and pip), and installs `ONCVPSP_PYONCVPSP` and `ONCVPSP_UPF_TOOLS`
  into it. `ONCVPSP_PYONCVPSP` is a pip requirement, by default the
  `oncvpsp-tests` branch of azadoks/pyoncvpsp, or a source directory, which
  is installed editable, so that changes to it need no reconfiguration:

  ```bash
  cmake --preset default -DONCVPSP_PYONCVPSP=$HOME/src/pyoncvpsp
  ```

  `ONCVPSP_UPF_TOOLS` is a pip requirement, `upf-tools==0.2.0` by default.
  The environment is made again when either changes; remove
  `build/test-venv` to update a pyoncvpsp taken from a repository.

The summary at the end of the configuration names the Python (`Test
Python:`). Without one that imports both, CMake warns and the tests that
compare outputs are disabled (ctest lists them as `Not Run (Disabled)`);
the others still run.

### Comparing outputs

`scripts/compare_output.py A.out B.out` parses two outputs and compares them
quantity by quantity, so that the comparison does not depend on how the
lines line up: the energies, eigenvalues, test configurations, model core,
log derivatives, diagnostics, plot data, input echoes and the psp8, UPF
(and with `--psml` PSML) files. It prints the largest differences for each
quantity and each one beyond its tolerance. A number passes if it agrees
within a relative or absolute 1e-5 (relative to |a| + |b|, the test of
`fldiff.pl -easy`, which the tests used before) and within the rule for its
quantity in `TOLERANCES`, which tightens it where 1e-5 says little, to about
two units of the last printed digit where the noise of the compiler flags
(2e-10 Ha in the all-electron energies) is below it (the terms of the
all-electron energy to 1e-8 Ha, the eigenvalues of the test configurations
to 1e-7 Ha, d2Exc to 3e-8 Ha, the input echoes to 1e-12); integers are
compared as numbers too, apart from the numbers of steps of the solvers,
which may differ by one, and text, with its numbers masked, must be the
same. With `--identical`, for outputs that
must be the same (the TOML and text inputs, the split and the single
program, a pipe and a file, the input echoes run again), the numbers must
agree to 1e-10, and the lines must be the same too, which covers what
pyoncvpsp does not parse.

The projectors and their energies are never compared: they are not unique
(any change of basis among the projectors of one l, or l and j, that keeps
sum_ij |beta_i> D_ij <beta_j| leaves the pseudopotential the same). The
nonlocal kernel of each block of the psp8, UPF and PSML files is compared
instead, which must agree to 1e-6 (1e-10 with `--identical`) in the
Frobenius norm, relative to the kernel, and for a spin-orbit kernel to the
scalar-relativistic kernel of its l; the local potential is compared point
by point, to the base test.

### References

The reference outputs are in `tests/refs`, compared with
`compare_output.py`; changes that move the numbers beyond the tolerances
need new references. Before updating them, check what changed:

- `scripts/compare_output.py tests/refs/NAME.out build/tests/NAME/NAME.out`
  (the reference tests leave this report in `build/tests/NAME/NAME.diff`).
- `scripts/kernel_diff.py OLD NEW`, which needs only the standard library,
  compares the pseudopotentials of two outputs, or psp8 or UPF files, by
  what a calculation uses, the nonlocal kernel of each block and the local
  potential.
- The `nist` tests (`scripts/nist_test.py`) compare the all-electron atom
  with the NIST atomic reference data to its stated accuracy (about 1e-6
  Ha), with the data's functionals and speed of light (iexc 6 and 7,
  `c_light`). The few cases that miss it by a little (`EXPECTED_FAILURES`:
  some ScRLDA lanthanides, which other codes miss too, and Fr) are
  expected to fail, and a test fails if one of them passes. ctest runs the
  light atoms, a few heavy ones and those of the expected failures; the
  whole table, 183 atoms and cations per set, is the script run by hand with
  `--jobs`.

`scripts/ref_patch.py NEW.out tests/refs/NAME.out` then updates a reference
minimally: it keeps every line within the base test (1e-5) and takes the new
one only where a number moves beyond it or the text changes. Run the
reference test afterwards, as a quantity with a tighter rule can still
differ.

## Editors

`.editorconfig` tells any editor that honours it the whitespace rules: LF
line endings, no trailing whitespace, one final newline, spaces not tabs,
2-column Fortran indentation, 80 columns. The indentation of a whole file
is `scripts/reindent.py`'s job, and an editor can run it as a formatter
because it reads standard input and writes standard output with `--stdin`.
In vim, for example, `:%!python3 scripts/reindent.py --stdin %` re-indents
the buffer.

### VS Code

The `.vscode` directory holds workspace settings, tasks and extension
recommendations, so opening the folder is most of the setup.

1. Install the tools above and build once with the `default` preset (the
   gfortran linter reads the module files from `build/`).

2. Open the folder. VS Code offers to install the recommended extensions
   (`.vscode/extensions.json`); accept, or install them from the
   Extensions view:

    - **Modern Fortran** (`fortran-lang.linter-gfortran`): syntax,
      navigation and hover through fortls, and compile diagnostics through
      gfortran. Its own formatters are disabled by the workspace settings
      so that only reindent formats the files.
    - **Custom Local Formatters** (`jkillian.custom-local-formatters`):
      registers `scripts/reindent.py --stdin` as the formatter for
      free-form Fortran.
    - **Fortitude** (`plasmafair.fortitude`): fortitude's diagnostics in
      the editor and its fixes as quick fixes, with the settings in
      `fortitude.toml`. It uses the `fortitude` on `PATH`, hence the
      pinned install above. The extension is in preview; if it misbehaves,
      the `pre-commit: this file` task runs the same checks.
    - **Ruff** and **Python** for the Python scripts, **EditorConfig** for
      `.editorconfig`, and **CMake Tools**, which picks up the presets.

3. That is all. The workspace settings (`.vscode/settings.json`) do the
   rest:

    - Free-form Fortran (`.f90`, `.F90`, `.inc`) is formatted on save by
      reindent, which gets the file's path and leaves the imported files
      (`functionals.F90`, `gg1cc.f90` etc., the list in `reindent.py`)
      unchanged. The Format Document command runs it by hand.
    - The gfortran linter gets the include paths of the `default` build
      and the `WITH_TOML` definition, so it compiles each file as the build
      does. Without a build in `build/`, expect "can't open module file"
      diagnostics; fortls and fortitude work without one.
    - fortls reads `.fortls` at the top of the tree: the sources in `src`
      and `tests/unit`, `.inc` files included, with `WITH_TOML` defined.
      Hovering over a routine, at its definition or at a `call`, shows its
      `!!` header: the description and the `@param` arguments. Hovering
      over a constant shows its `!>` or `!<` comment (see the Header
      rules in [Coding Standards](coding_standards.md#conventions)).
    - The Python files are formatted on save by ruff.

   The tasks (Terminal > Run Task) are `pre-commit: this file`,
   `pre-commit: all files`, `fortitude --fix then reindent: this file` and
   `check interfaces`.

The settings are per-workspace, so they do not change other projects.
Personal settings for this clone go in the user settings, or in
`.vscode/settings.json` without committing them.

## Documentation

The site under `doc/` is built with [MkDocs](https://www.mkdocs.org) and
deployed by CI from `master`. To preview it:

```bash
pip install -r requirements.txt
mkdocs serve
```

New pages go in the `nav` of `mkdocs.yml`.
