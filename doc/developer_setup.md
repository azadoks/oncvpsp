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
`CMakeUserPresets.json` (ignored by git) is the place for personal presets
such as compiler paths or `Libxc_ROOT`.

The legacy make build (`make.inc` from `make.inc.example`, then `make all`
in the top directory) must keep working. Before pushing, check that it still
builds, in a scratch copy of the tree so that `src/` is not littered with
`.o` and `.mod` files. When adding or removing a source file, update both
`src/CMakeLists.txt` and `src/Makefile` (`OBJS_ALL` and the module
dependencies).

After changing any routine's argument list, run
`scripts/check_interfaces.sh` (it needs a CMake build in `build/`), which
checks every call against the definition of the routine it calls.

## Tests

`ctest --preset default` runs everything; `-L <label>` selects a part of it
(`unit`, `sr`, `r`, `nr`, `split`, `mesh`, `toml`, `psml`, `pipe`,
`errors`, `echo`, `psp8`) and `-R <regex>` a test by name. The reference
outputs are in `tests/refs` and are compared as written, so changes that
alter the numbers need new references. Each test leaves its outputs in its
own directory under `build/tests`.

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
