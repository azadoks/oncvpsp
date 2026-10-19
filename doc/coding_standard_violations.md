# Coding standard violations

Style and consistency deviations found by auditing the 82 Fortran files in
`src/` that `src/CMakeLists.txt` compiles (including the optional Libxc/TOML
files and the two main programs). The vendored `src/xmlf90-wxml/` library is
excluded. Every item has since been fixed or kept on purpose (see
[Status](#status)); the rest of the page describes the code as it was.
Line numbers refer to commit `7e9480a`, before any Known Issues fixes.

This page covers layout, naming, comments and consistency only. Suspected bugs
are in [Known Issues](known_issues.md). Language-level practice (`implicit
none`, `GOTO`, `COMMON`, explicit interfaces, kinds, pointers, error handling)
belongs in `fortran_best_practices.md`. Where a finding touches those areas it
is only cross-referenced here.

## The standard

There is an **explicit** standard, [Coding Standards](coding_standards.md)
(`doc/coding_standards.md`). It is short and mostly about structure: free-form
Fortran 95, one subroutine per file, no one-letter variable names,
input/output/local variables grouped, a short description next to the
`subroutine` statement, and every argument described in comments after it. It
ends with "There are probably other conventions I've been following". There is
no `CONTRIBUTING`, `.editorconfig`, fprettify/findent/fortitude config,
pre-commit hook or CI lint step. CI (`.github/workflows/cmake_ctest.yaml`)
only builds and runs the tests.

Everything the explicit standard leaves out is **inferred** here from the
majority of the Hamann-authored files, measured with scripts (counts below).
The inferred conventions are:

| Topic | Inferred convention | Majority |
|---|---|---|
| Header | `!` + GPLv3 block, "Copyright (c) 1989-2019 by D. R. Hamann…" | 68/82 files |
| Unit statement column | column 2 (one leading space) | 64/82 files |
| Block indent step | 2 spaces | 1025/1821 block bodies (56%) |
| Continuation | trailing `&` plus a leading `&` in column 1 | 843 continuations |
| Keyword case | lower case | all but 2 files |
| Block ends | `end do`, `end if`, `else if`, `end subroutine name` | 808/839, 612/637, 88/90, 110/119 |
| Relational operators | `==`, `/=`, `<`, … | 880 vs 91 uses |
| Real declarations | `real(dp)` with a local `integer, parameter :: dp=kind(1.0d0)` | 880/917 decls |
| Declarations | always `::` | all but 14 lines |
| Section labels | `!Input variables`, `!Input/output variables`, `!Output variables`, `!Local variables`, in that order | 44 routines exact |
| Last statement | `return` before `end subroutine` | 84/119 units |
| Output | `write(6,'(…)')` with an inline format string; no labeled `FORMAT` | 951 inline, 0 labeled |
| Line length | ≤ 80 columns in practice (no hard limit) | 167 lines over 80, 0 over 132 |
| Scalar names | doubled letters (`ii`, `jj`, `ll`, `rr`, `zz`) instead of one letter | — |

The newer module files (`input_text_m.f90`, `input_toml_m.f90`,
`functionals.F90`, most of `m_psmlout.f90`) follow a second, fprettify-like
style: column-1 unit statements, 3-space indents, trailing-only `&`,
one-per-line declarations with `!>` doc comments. Imported code
(`exc_off_pbe.f` from Perdew's group, `gg1cc.f90`/`gp1cc.f90`/`gpp1cc.f90` from
ABINIT, `functionals.F90` from Octopus, `m_libxc_list.F90` from SIESTA) keeps
its upstream style. `exc_off_pbe.f` says so itself (lines 21–22). Most of the
violations below come from these two groups. The project has since
exempted the imported code and converted the module files to the legacy
style (see [Status](#status)).

## Status

The mechanical pass (October 2026) used
[fortitude](https://fortitude.readthedocs.io) 0.9.2, which now runs as a
pre-commit hook and in CI (see
[Coding Standards](coding_standards.md#linting)). It covers every file
except the imported ones (`exc_off_pbe.f`, `functionals.F90`,
`gg1cc`/`gp1cc`/`gpp1cc`, `m_libxc_list.F90`). fprettify was tried first,
but its 0.3.7 release applies the operator replacements only while
reformatting whitespace, which would restyle every legacy file. No output
changed: the make build gives bit-identical outputs before and after.

Fixed:

- trailing whitespace, in all files (tabs remain only in `exc_off_pbe.f`);
- `enddo`, `endif` → `end do`, `end if`;
- `.eq.`-style relational operators → `==` etc.;
- `character*N` → `character(len=N)`, and `::` in every declaration;
- `CONTAINS` → `contains`;
- the specific intrinsics `dlog`, `dsqrt`, `dabs`, `dexp`, `dmin1`,
  `dmax1`, `idint` and `dfloat` → `log`, `sqrt`, `abs`, `exp`, `min`,
  `max`, `int` and `real(…,dp)`;
- `intent` on every dummy argument, with the missing section labels in
  `ae_file_m.f90` and `vkboutwf.f90`;
- `only` on every `use` statement;
- `implicit none` in `m_psmlout`, and its `lsymb` is a `parameter`;
- assumed-size dummies (`*`) in `dpnint`, `vkbphsft` and `vkboutwf`
  given explicit shapes (those in `adams.inc` are kept, so that the
  integrators can inline its functions);
- `(/ … /)` → `[ … ]`, superfluous semicolons, a `case default` in each
  `select case`, and lines over 100 columns;
- the section labels, spelled as the four canonical ones (`!Input
  variables`, `!Input/output variables`, `!Output variables`, `!Local
  variables`, and `!Output variables - printing only`), and the typos
  `vafriables`, `vaiables`, `OUTPuT` and `1989-201r`;
- the stray first lines, and `run_diag_sr_so_r.f90`'s copyright line;
- the [cosmetic Known Issues](known_issues.md#cosmetic).

A second pass (October 2026, commits `f5a2ab7` to `f7f76ca`) settled the
items under "Needs manual work or a decision" in [Fixing](#fixing). The
conventions it applied are now written down in
[Coding Standards](coding_standards.md#conventions). It fixed:

- **Indentation**, in Hamann's layout in every file, the module files
  included: unit statements in column 2, 2-space steps, comments in
  column 1, a leading `&` in column 1. `scripts/reindent.py` does this
  and runs as a pre-commit hook. The objects compiled byte-identical.
- **Line length**: at most 80 columns, enforced by fortitude.
- **Headers**: a description for every routine and a comment for every
  dummy argument; the stale argument comments corrected or removed.
- **Section labels**: canonical labels in canonical order in every
  routine, with the misfiled declarations moved.
- **Names**: no one-letter names; `m_psmlout` renamed `psmlout_m`;
  `double precision` in it replaced by `real(dp)`.
- **Dead code**: the commented-out code removed; the alternative built-in
  meshes in the main programs kept as prose. `oncvpsp_r.F90`'s `cmtskp`
  and `read_error` and the uncompiled files had already gone.
- **Statement forms**: `if(`, `write(` etc. without a space, and `return`
  before the end of every procedure.
- **Errors**: the command-line errors of the drivers in the documented
  form on unit 6, `stdout` replaced by unit 6, and `error stop 1` in
  `modcore_budget` (which exited with status 0), `ae_file_m` and
  `input_toml_m`.

The make build's outputs stayed identical, and the parse trees changed
only where the renames, the error stops and the added `return`s were.

Kept on purpose:

- the imported files, which keep their upstream form;
- copyright and license text: the years in `linout.f90`, `ldiracf*.f90`
  and `exc_off_pbe.f`, and the missing headers of `input_text_m.f90` and
  `input_toml_m.f90`, whose authors would have to add them;
- the text of messages from 4.0.1, including the `test_data:` prefix in
  `check_data`, the `_r` routines' non-`_r` prefixes, and the
  list-directed `lmax` and `calling psmlout` prints, which the reference
  outputs contain (see [Error messages](coding_standards.md#error-messages));
- `derivs` in `exc_libxc.f90`, an external routine that `functionals.F90`
  calls too, so it can't become internal, and a file of its own would need
  a change to every user's `make.inc` (`OBJS_LIBXC`);
- unused locals (e.g. `maxset` in `lschpsbar`, `cfloc` in `check_data`),
  which fortitude's default rules don't flag.

The pass also found a bug, now in
[Known Issues](known_issues.md#p2-open).

## Summary

| Category | Lines / items | Files |
|---|---|---|
| More than one top-level unit per file | 5 files (17 units) | 5 |
| File name ≠ unit name | 3 | 3 |
| `.F90` without preprocessor directives | 1 | 1 |
| Tracked sources not compiled | 3 | 3 |
| License header missing / truncated / non-standard | 11 | 11 |
| Stray first line before the header | 5 | 5 |
| Copyright year inconsistent or mistyped | 7 | 7 |
| Routine without a description comment | 20 routines | 4 |
| Arguments not described in comments (rule 9) | 233 of 1385 args, 63 routines | 38 |
| Argument comments naming non-arguments (stale) | 23 | 16 |
| Section labels missing entirely | 27 routines | 7 |
| Section labels misspelled / non-canonical | 40 routines | 38 |
| Declarations under the wrong section label | 18 | 6 |
| Locals declared under an argument label | 23 | 17 |
| Trailing whitespace | 625 (310 blank-but-spaces) | 70 |
| Tabs | 5 | 1 |
| Lines > 132 (free) / > 72 code (fixed) | 0 | 0 |
| Lines > 80 (de facto limit) | 167 | 28 |
| CRLF / missing final newline | 0 | 0 |
| Unit statement not in column 2 | 17 | 17 |
| Mixed indent step within a file | — | 55 of 77 free-form files |
| Block end misaligned with its opener | 29 | 13 |
| Body not indented under its opener | 12 | 8 |
| Continuation style mixed within a file | — | 8 |
| Upper-case keywords | 37 | 2 |
| Bare `end` for a procedure | 9 | 2 |
| `enddo` / `endif` / `elseif` | 31 / 25 / 2 | 1 / 3 / 1 |
| `.eq.`-style relational operators | 89 | 41 |
| Declaration without `::` | 14 | 8 |
| `real(8)` / `real*8` / `double precision` | 37 | 3 |
| `character*N` instead of `character(len=N)` | 31 | 9 |
| One-letter variable names (rule 5) | 50 declared (+ all of `exc_off_pbe.f`) | 8 |
| Module naming (`m_` / `_m` / none) | 7 modules, 3 schemes | 4 |
| Commented-out code | ~106 lines | 31 |
| Dead, never-called routines left in place | 2 | 1 |
| Error-message prefix ≠ routine name | 52 | 14 |
| List-directed (`*`) debug output | 4 | 4 |
| Hard-coded unit 6 mixed with `stdout`/`stderr` | — | 3 |
| Sibling (`foo` vs `foo_r`) divergences | 9 items | 14 |

## File naming and contents

Standard: rule 2, one subroutine per file, file named after it; `.f90` for
free form, `.F90` only when the file needs the preprocessor.

- **More than one top-level unit:**
  - `aeo.f90` — `aeo`, `aio`, `aei`, `aii` (lines 23, 36, 49, 62).
  - `exc_libxc.f90:118` — external `derivs` after `exc_libxc`.
  - `exc_off_pbe.f` — `excggc` plus 7 more (`exchpbe:201`, `corpbe:274`,
    `gcor2:417`, `exchpw91:436`, `corlsd:481`, `gcor:516`, `corpw91:536`).
  - `functionals.F90` — three modules: `parser:21`, `xc_messages:41`,
    `functionals_m`.
  - `oncvpsp_r.F90:786,807` — `cmtskp` and `read_error` after the main
    program (see [dead code](#comments-and-dead-code)).
- **File name ≠ unit name:** `exc_libxc_stub.f90:20` (defines `exc_libxc`, as
  intended for a stub), `exc_off_pbe.f:24` (`excggc`), `functionals.F90:21`
  (first unit is `parser`, main module is `functionals_m`).
- **Extension:** `m_libxc_list.F90` uses `.F90` but has no preprocessor
  directives. `exc_off_pbe.f` is the only fixed-form file, which contradicts
  rule 1 ("free-format").
- **Tracked but not compiled:** `fpovlp_r.f90` (only referenced from a
  commented-out call, `oncvpsp_r.F90:522`), `oncpsp_nr.f90`, `oncvpsp_nr.F90`.
  The "Fully-relativistic files" block of `src/CMakeLists.txt` lists
  `fpovlp.f90` and `relatom.f90` a second time (lines 71, 75). `fpovlp.f90` is
  probably meant to be `fpovlp_r.f90`.
- Related, outside this page's scope: `src/CMakeLists.txt` builds
  `oncvpspnr.x` from `oncvpsp.F90` with `RELATIVISTIC=0`. No source references
  `RELATIVISTIC`, and `oncvpsp.F90:114` hard-codes `srel=.true.`, so the "non-
  relativistic" executable looks like a copy of the scalar-relativistic one.
  This belongs in [Known Issues](known_issues.md).

## Header blocks

Standard: `!` + "Copyright (c) 1989-2019 by D. R. Hamann…" + GPLv3 text,
then a description (rule 8), argument comments (rule 9) and grouped
declarations (rule 6).

### License / copyright

- **No copyright line:** `functionals.F90`, `input_text_m.f90`,
  `input_toml_m.f90`. In `run_diag_sr_so_r.f90` the block is truncated: line 1
  is `! University`, so the copyright line is missing.
- **No GPL text:** `input_text_m.f90`, `input_toml_m.f90`, `m_libxc_list.F90`
  (copyright only, line 6).
- **Different license blocks (imported):** `gg1cc.f90`, `gp1cc.f90`,
  `gpp1cc.f90` use the ABINIT robodoc header ("GPL, see ~abinit/COPYING",
  no version). `functionals.F90` has Octopus' GPLv2+ text.
- **Stray first lines:** `run_ghosts.f90:1` (`!:`), `sbf_basis.f90:1` and
  `sbf_basis_con.f90:1` (`!_der`), `run_vkb.f90:1` and `run_vkb_r.f90:1`
  (blank line).
- **Copyright years** (68 files say 1989-2019): `linout.f90:2` 1989-2017
  (sibling `linout_r` says 2019); `ldiracfb.f90:2`, `ldiracfs.f90:2`
  1989-2012; `exc_off_pbe.f:2`, `m_psmlout.f90:2` 1989-2014;
  `sratom.f90:2` reads **`1989-201r`** (typo). The `copyright:` line in
  `mkdocs.yml` says 2020-2025.

### Description comment (rule 8)

20 routines have no description next to the `subroutine`/`function`
statement:

- `aeo.f90:36,49,62` — `aio`, `aei`, `aii` (the file comment covers `aeo`
  only).
- `functionals.F90:28,33,46,54,63,71,76,81` — `real_parser`,
  `logical_parser`, and the six `messages_*` routines.
- `input_toml_m.f90:14,399,784,794,804` — `read_input_toml`,
  `read_input_toml_r`, `toml_table_check_associated`,
  `toml_array_check_associated`, `check_len`.
- `m_psmlout.f90:491,1096,1118,1164` — `do_configuration` (×2),
  `my_add_attribute`, `write_psml_item`.

`vpinteg.f90:24` has its description after `implicit none` and the `dp`
declaration, unlike every other routine.

### Argument comments (rule 9)

233 of 1385 dummy arguments, in 63 of 109 routines across 38 files, have no
`!name  description` line, `!>` doc comment or inline comment. The worst:

- `check_data.f90:19` — 0 of 32 documented ("see annotated input document").
- `m_psmlout.f90:33,515` — `psmlout` 23/46 and `psmlout_r` 27/49
  undocumented; `write_psml_item` (`:1164`) 13/13.
- `exc_libxc.f90:20,118`, `exc_libxc_stub.f90:20` — 0 of 7 each.
- `vpinteg.f90:19` (0/6), `dpnint.f90:19` (1/6),
  `functionals.F90:165` (0/6).
- Single omissions are common in otherwise complete lists:
  `psatom.f90`/`psatom_r.f90` (`fat`, `rhoc`, `rho`, `ierr`),
  `run_config.f90`/`run_config_r.f90` (`facnf`),
  `gnu_script.f90`/`gnu_script_r.f90` (`evkb`), `wellstate.f90` (`nnin`),
  `sbf_basis.f90` (`nconmx`).

**Stale argument comments.** These document a name that is not an argument
(23 lines, 16 files). In most cases the argument was renamed and the comment
was not updated:

- `psatom.f90:27-28`, `psatom_r.f90:27-28` — `!fa` (argument is `fat`),
  `!rpk` (no such argument). `psatom.f90` also documents `!uua`.
- `wellstate.f90:24`, `wellstate_r.f90:25` — `!nn` (argument is `nnin`).
- `vout.f90:29` — `!vc` (argument is `vxc`).
- `vrel.f90:24,26` — `!nn`, `!ierr` (neither is an argument).
- `run_phsft.f90:31,43`, `run_phsft_r.f90:45` — `!ep`, `!irphs`.
- `run_optimize.f90:32,34,40` — `!uu`, `!uu2`, `!np_in`.
- `eresid.f90:31,32,47` — `!dr`, `!dq`, `!eresidevec`.
- `sr_so_r.f90:39`, `run_diag_sr_so_r.f90:35` — `!esol`.
- `relatom.f90:41` (`!srel`), `lschkb.f90:33` (`!zz`),
  `optimize.f90:42` (`!orbasis_it`), `functionals.F90:720` (`m`).

### Section labels (rule 6)

Of 111 procedures (excluding the 2 programs), 44 use exactly the four
canonical labels, 40 use at least one non-canonical label, and 27 have none.

- **None (27):** `aeo.f90` (4), `exc_libxc_stub.f90`, `functionals.F90` (13),
  `input_toml_m.f90:784,794,804`, `m_psmlout.f90` helpers
  (`:491,1096,1118,1128,1147,1164`).
- **Non-canonical spellings** (one example each; 38 files in all):
  - Case: `!Local Variables` (`ldiracfb.f90:52`, `ldiracfs.f90:53`,
    `lschpb.f90:48`, `lschpsbar.f90:55`), `!Input Variables`
    (`lschvkbb.f90:45`, `lschvkbbe.f90:47`), `!Input/Output variables`
    (`run_vkb.f90:53`, `run_vkb_r.f90:53`), `!input variables`
    (`const_basis.f90:41`), `!local variables` (`sbf_rc_der.f90:40`,
    `wf_rc_der.f90:49`).
  - Spacing: `! Input variables` (`check_data.f90:32`, `der2exc.f90:44`,
    `vploc.f90:39`, `input_toml_m.f90:25`).
  - Singular: `!Output variable` (`fpovlp.f90:46`, `vpinteg.f90:33`,
    `renorm_r.f90:39`), `!Local variable` / `!In/Out variable`
    (`input_text_m.f90:468,472`, `oncvpsp_r.F90:790,793`).
  - Typos: `!Local vafriables` (`excwig.f90:44`), `!Input vaiables`
    (`vout.f90:42`), `!OUTPuT` (`qroots.f90:39`).
  - Other schemes: `!INPUT`/`!OUTPUT`/`!IN/OUT` + `!Arguments`
    (`qroots.f90:30`, `pspot.f90:28`, `sbf_basis.f90:28`,
    `sbf_basis_con.f90:29`, `sbf_rc_der.f90:26`, `wf_rc_der.f90:23`); ABINIT
    `!Arguments ----` / `!scalars` (`gg1cc.f90:39`, `gp1cc.f90:44`,
    `gpp1cc.f90:37`); `! local`, `! local vars` (`exc_libxc.f90:35,130`).
  - Three spellings of "printing only": `!Output variables - printing only`
    (12), `!Output -- printing only` (`gnu_script.f90:37`,
    `gnu_script_r.f90:37`), `!Output variables  only printing`
    (`run_config.f90:62`, `run_config_r.f90:62`).
- **Order:** Input > [Input/output] > Output > Local everywhere except
  `renorm_r.f90:34-45` (Input > Output > Input/Output > Local).
  `relatom.f90:62` and `sratom.f90:63` add a `!Local function` block before
  `!Local variables`. `pspot.f90`, `qroots.f90`, `sbf_basis.f90` and
  `sbf_basis_con.f90` put `!INPUT`/`!OUTPUT` labels on the argument-comment
  block and then repeat the canonical labels on the declarations. Prose
  comments that begin with "output" or "input" (`psatom.f90:187` "! output
  potential", `linout.f90:123`) are easily mistaken for section labels.
- **Declarations under the wrong label** (intent ≠ section, 18 lines):
  - `upfout.f90:78-83`, `upfout_r.f90:77-82` — `intent(in)` arrays under
    `!Input/output variables`, after a free-text "additional input" comment.
  - `m_psmlout.f90:72,555` — `intent(inout) :: lmax` under `!Input variables`
    (cause in [Known Issues](known_issues.md), `lmax = lmax`).
  - `optimize.f90:57` — `intent(out) :: eresq` under `!Input variables`.
  - `input_text_m.f90:491`, `oncvpsp_r.F90:813` — `intent(inout) :: inline`
    under `!Input variables`. `input_text_m.f90:470` — `intent(in) :: unit`
    under `!In/Out variable`.
- **Locals under an argument label** (no `intent`, not a dummy argument; 23
  lines):
  - `nlim`: `run_diag.f90:50`, `run_diag_r.f90:48`,
    `run_diag_sr_so_r.f90:50`, `run_plot.f90:55`, `run_plot_r.f90:55`.
  - `srel`: `fphsft_r.f90:47`, `modcore.f90:56`, `modcore2.f90:56`,
    `modcore3.f90:56`, `run_config_r.f90:60`, `run_plot_r.f90:62`.
  - `zz`: `lschkb.f90:43`, `lschps.f90:41` (unused; see Known Issues).
  - `npx`, `lpx`: `run_plot.f90:57`, `run_plot_r.f90:57`.
  - `lwork`: `run_vkb.f90:47`, `run_vkb_r.f90:47`.
  - `orbasis_it`, `eres`: `optimize.f90:55,63`.
  - `vkb`, `evkb`: `run_diag_sr_so_r.f90:54,58`.
  - `xx`: `sbf_rc_der.f90:38`. `nder`: `wf_rc_der.f90:47`.

## Whitespace and layout

### Trailing whitespace, tabs, line endings, length

- **Trailing whitespace:** 625 lines in 70 files. 310 of these are blank lines
  containing only spaces. Worst: `lschfb.f90` (40), `lschpb.f90` (37),
  `ldiracfb.f90` (32), `lschvkbb.f90` (30), `lschpse.f90` (28),
  `exc_libxc.f90` (28).
- **Tabs:** only in `exc_off_pbe.f:225,227,228,245,247`. Lines 245 and 247 are
  code; a leading tab in fixed form is a compiler extension.
- **CRLF / final newline:** none.
- **Length:** no free-form line exceeds 132 characters. `exc_off_pbe.f:2,128,485`
  exceed 72, but they are comments. The Hamann files keep to about 80 columns;
  167 lines in 28 files exceed that, mostly `functionals.F90` (40),
  `input_toml_m.f90` (29), `gnu_script_r.f90` (17), `oncvpsp_r.F90` (10).
  The longest is `functionals.F90:610` (127).

### Indentation

Standard (inferred): unit statement in column 2, top-level statements in
column 2, each block body indented 2 more. Comment lines start in column 1.

- **Unit statement column.** Column 1 in 14 files: the modules and imported
  files, plus `run_phsft.f90`, `run_phsft_r.f90`, `sbf8.f90` and
  `wellstate_r.f90` (sibling `wellstate.f90` uses column 2). Column 3 in
  `run_config.f90` and `run_config_r.f90`; column 4 in `run_optimize.f90`.
- **Step width.** Of 1821 block bodies, 1025 use 2 spaces, 442 use 1 and 319
  use 3. Only 22 of 77 free-form files are uniform. In the legacy files 1-space
  and 2-space bodies alternate within the same routine, for both `do`
  (251×1, 462×2) and `if … then` (191×1, 558×2). For example, at
  `lschfb.f90:77-78` the `if` is in column 2 and its body in column 3, while
  the `do` bodies in the same routine step by 2. Heavily mixed: `modcore3.f90` (33×2, 30×1), `run_vkb_r.f90`
  (28×2, 24×1), `check_data.f90` (46×1, 3×2), `sr_so_r.f90` (29×1, 7×2).
  The module files use 3.
- **Block end not aligned with its opener** (29 lines, 13 files), e.g.
  `eresid.f90:97`, `qroots.f90:97` (closes the line-60 `if` one column off),
  `run_config.f90:119-120` and `run_config_r.f90:119-120`, `psatom_r.f90:194`,
  `wellstate.f90:167`, `wellstate_r.f90:175`, `gnu_script_r.f90:287`,
  `run_plot_r.f90:220`, `m_psmlout.f90:885,901,1045,1079,1084`,
  `functionals.F90:619`, `oncvpsp_r.F90:719,721`. In `oncvpsp.F90:159-198` and
  `oncvpsp_r.F90:160-196` the `case` labels are 3 columns right of
  `select case`.
- **Body not indented under its opener** (step 0): `run_vkb.f90:121,151`,
  `run_vkb_r.f90:153,184`, `run_plot.f90:101`, `run_ghosts.f90:142`,
  `modcore3.f90:151`, `m_psmlout.f90:1003`, `functionals.F90:517`,
  `oncvpsp_r.F90:715,752,763`.
- **Continuation lines.** The legacy convention puts `&` at the end of the line
  and again in column 1 of the next. `functionals.F90`, `input_text_m.f90`,
  `input_toml_m.f90` and `m_libxc_list.F90` use a trailing `&` only. These
  files mix both: `m_psmlout.f90` (21 leading / 47 trailing-only),
  `oncvpsp.F90` (51/27), `oncvpsp_r.F90` (68/24), `ldiracfs.f90` (1/5),
  `run_phsft.f90`, `run_phsft_r.f90` (5/2), `upfout.f90`, `upfout_r.f90`
  (one trailing-only each). In 19 lines the leading `&` is indented rather
  than in column 1.

## Keywords and statement forms

- **Upper-case keywords:** `exc_off_pbe.f`, 34 lines (`SUBROUTINE`,
  `IMPLICIT REAL*8 (A-H,O-Z)`, `CALL`, `RETURN`, `END`, from line 201 on).
  `m_psmlout.f90:31,490,1095` use `CONTAINS`. Upper-case format descriptors
  (`'(A,A)'`) appear in `input_toml_m.f90` (5 lines, e.g. 789).
- **Bare `end` for a procedure** (standard: `end subroutine name`, 110/119):
  `exc_off_pbe.f:270,413,432,477,512,532,629`, `input_toml_m.f90:792,802`.
- **`enddo`** (standard `end do`, 808 of 839): 31 lines, all in
  `m_psmlout.f90` (e.g. 885, 901).
  **`endif`** (standard `end if`, 612 of 637): 25 lines, `m_psmlout.f90` (23,
  e.g. 237, 253), `linout.f90:188`, `exc_off_pbe.f:248`.
  **`elseif`**: `functionals.F90:639,657`.
- **Old relational operators** (`.eq.`, `.ne.`, `.lt.`, `.le.`, `.gt.`,
  `.ge.`): 89 lines in 41 files (91 uses, against 880 symbolic). Most files
  have one or two strays (e.g. `lschfb.f90:180`, `psatom.f90:162`,
  `run_config.f90:184`). Concentrated in `vploc.f90` (14, e.g. 75–125),
  `vout.f90` (12), `exc_off_pbe.f` (5), `m_psmlout.f90` (4), `modcore.f90` (4).
- **Spacing:** `if(` (732) against `if (` (160). The module files use `if (`;
  `oncvpsp.F90`, `oncvpsp_r.F90`, `upfout*.f90`, `vout.f90`, `vploc.f90`,
  `run_vkb_r.f90` and `run_phsft*.f90` mix both. `real(dp),intent(in)`
  without a space appears 6 times, against 715 with one (all in
  `gg1cc`/`gp1cc`/`gpp1cc`).
- **`return` before `end`:** 84 of 119 procedures end with `return`. The
  exceptions are the module/imported routines plus `der2exc.f90:104`. This is
  harmless, but the files disagree.

## Declarations

- **No `::`** (14 lines, 8 files): `dpnint.f90:40`, `gnu_script.f90:40`,
  `qroots.f90:52`, `vkboutwf.f90:54`, `vout.f90:55`,
  `m_psmlout.f90:118,601,1134,1153,1154`, `oncvpsp_r.F90:794`,
  `exc_off_pbe.f:34-36`.
- **Real kind spelling** (standard `real(dp)`, 880 lines): `real(8)` in
  `functionals.F90` (30 lines, e.g. 29, 107, 369–393); `real*8` in
  `exc_off_pbe.f:32-35`, plus implicit `REAL*8 (A-H,O-Z)` in its other 7
  routines; `double precision` in `m_psmlout.f90:1131,1132,1135`. How `dp`
  is defined also varies: a local `integer, parameter :: dp=kind(1.0d0)` in 79
  procedures, against `iso_fortran_env` `real64` in the modules. Real
  literals are written `1.0d0` (2223) rather than `1.0_dp` (38). Which kind
  scheme to standardise on is a `fortran_best_practices.md` question.
- **`character*N`** (31 lines, 9 files) against `character(len=N)` (52):
  `check_data.f90:33-34`, `linout.f90:68,77,87`, `linout_r.f90:73,82,94`,
  `upfout.f90:70,83,93,94`, `upfout_r.f90:69,82,94,95`,
  `m_psmlout.f90:82,91,111,112,118,566,575,594,595,601`,
  `oncvpsp.F90:93-94`, `oncvpsp_r.F90:94-95,794`.
- **Grouped declarations.** The legacy files group several names per line
  (`real(dp) :: de,emax,emin`), and the modules put one name per line with a
  `!>` comment (749 grouped, 874 single lines in all). Neither standard
  addresses this. The two styles are consistent within each file group.

## Naming

- **One-letter variable names** (rule 5): 50 declarations in 8 files, mostly
  imported code: `functionals.F90` (18, e.g. `a, b, c` at 380,
  `n, s, l, t` at 388), `m_psmlout.f90` (19, e.g. `i, j` at 102, `f, r` at
  1131), `input_toml_m.f90` (6, `i`, `l`), `exc_libxc.f90:43,131`
  (`v`, `e`, `i`), `oncvpsp.F90:118`, `oncvpsp_r.F90:119` (`i`),
  `sratom.f90:74` (`u`), `vploc.f90:49` (`x`). `exc_off_pbe.f` uses
  one-letter names throughout.
- **Module naming:** three schemes for 7 modules. Prefix: `m_psmlout`,
  `m_libxc_list`. Suffix: `input_text_m`, `input_toml_m`, `functionals_m`.
  None: `parser`, `xc_messages` (both inside `functionals.F90`).
- **`_r` suffix.** It consistently means "fully relativistic" for external
  routines (`fphsft_r`, `run_vkb_r`, …) and module entry points (`psmlout_r`,
  `read_input_text_r`, `read_input_toml_r`), with these exceptions:
  `m_psmlout.f90:1128` `add_zero_r` uses `_r` to mean "add an r=0 point",
  and the Dirac solvers `ldiracfb`/`ldiracfs` and the SO helper `lschvkbbe`
  are relativistic-only but carry no `_r`. `sr_so_r`, `renorm_r` and
  `run_diag_sr_so_r` have no non-`_r` sibling. That is fine, but it means the
  suffix marks the build target rather than a paired file.

## Comments and dead code

- **Commented-out code:** about 106 lines in 31 files, found by heuristic
  after removing prose false positives. Mostly disabled debug `write`s and
  superseded formulas. Largest:
  - `const_basis.f90:75-85,115` (10) — SVD debug dumps.
  - `sbf_basis.f90:110-111,147-149,177-179` and
    `sbf_basis_con.f90:108-109`.
  - `oncvpsp.F90:115,234-241,596-598` (`!srel=.false.`, alternative
    `amesh` values, `!if(.false.) then`) and
    `oncvpsp_r.F90:116,231-234,249,522`.
  - `lschvkbb.f90:77,116,225-242`, `modcore2.f90:221-229`,
    `modcore3.f90:234,296,328,364,395`, `optimize.f90:105-107,125`,
    `vkboutwf.f90:125-128` (four alternative `if` tests),
    `wellstate_r.f90:156-161`, `check_data.f90:120-124`,
    `psatom.f90:116-118`, `exchdl.f90:57-59`, `exc_libxc.f90:53-54,67-70`,
    `functionals.F90:635,653`.
- **Dead routines kept in the source:** `oncvpsp_r.F90:786` `cmtskp` and
  `:807` `read_error` are never called (input moved to `input_text_m.f90`,
  which has its own private copies at 465/486). See also `fpovlp_r.f90`
  [above](#file-naming-and-contents), and `sbf_basis.f90` (not called; see
  Known Issues).
- **Stale imported boilerplate:** `gg1cc.f90:22-23`, `gp1cc.f90:29-30`,
  `gpp1cc.f90:23-24` list `PARENTS psp1cc`, an ABINIT routine that does not
  exist here. `gg1cc.f90:32-33` says "This section has been created
  automatically by the script Abilint (TD). Do not modify the following lines
  by hand."
- **Open `FIXME`s:** `functionals.F90:217,234,238,252,260` (upstream
  Octopus notes).
- **Error-message prefix ≠ routine name** (the convention is
  `'routine: ERROR …'`; 52 lines, 14 files):
  - `check_data.f90` prints `test_data:` on all 33 messages (e.g. 54, 60, 65).
  - `_r` siblings print the non-`_r` name: `run_vkb_r.f90:177,231,272`
    (`run_vkb:`), `upfout_r.f90:351` (`upfout:`), `oncvpsp_r.F90:202,463`
    (`oncvpsp:`), `run_diag_sr_so_r.f90:170` (`run_diag_r:`).
  - Others: `run_config.f90:186`, `run_config_r.f90:189,204` (`runconfig:`),
    `sratom.f90:127` (`sratom123:`), `modcore.f90:91` (`rhomod:`),
    `sbf_basis.f90:120` (`sbfbasis:`), `sbf_basis_con.f90:104`
    (`sbf_basis:`), `vkboutwf.f90:104` (`vkbout:`),
    `exc_libxc_stub.f90:34` (`exc_libxc_stub:`),
    `functionals.F90:396,640,646`. Line 646 also misspells both words:
    `'functrionals: ERROOR'`.
- **Comment spacing:** `! text` (2097) and `!text` (1705) are both common. The
  legacy files use `!name  description` for argument lists and `! text` for
  prose, so this is not counted as a violation.

## Output and format statements

- **Format style is consistent:** there are no labeled `FORMAT` statements.
  All 951 formatted writes use an inline format string.
- **List-directed debug output** (prints that look left over from
  debugging): `run_plot.f90:198` (`write(6,*) 'lmax',lmax`),
  `oncvpsp.F90:648` and `oncvpsp_r.F90:772` (`print *, 'calling psmlout'`),
  `exc_libxc.f90:77` (`write (*,*)`).
- **Output unit.** 831 writes hard-code unit 6, while the newer code uses the
  named `stdout`/`stderr` from `iso_fortran_env`. `oncvpsp.F90` (40×`6`,
  9×named) and `oncvpsp_r.F90` mix both. `input_text_m.f90:495-496` writes
  to 6 although it imports `stdout`.

## Sibling consistency (`foo` vs `foo_r`)

The 13 compiled pairs (`fphsft`, `gnu_script`, `linout`, `psatom`,
`run_config`, `run_diag`, `run_phsft`, `run_plot`, `run_vkb`, `upfout`,
`wellstate`, `oncvpsp`, and in-module `psmlout`, `read_input_text`,
`read_input_toml`) mostly agree on header, section labels and layout. The
divergences:

- `linout.f90:2` copyright 1989-2017, `linout_r.f90:2` 1989-2019.
- `wellstate.f90` unit in column 2, `wellstate_r.f90` in column 1.
- `oncvpsp_r.F90` carries the dead `cmtskp`/`read_error`; `oncvpsp.F90`
  does not.
- `_r` error messages that name the non-`_r` routine (listed
  [above](#comments-and-dead-code)).
- `run_config_r.f90:60` and `run_plot_r.f90:62` declare a local `srel` under
  `!Input variables`; the non-`_r` versions do not have it.
- `run_vkb.f90` (28×2-space, 16×1-space steps) and `run_vkb_r.f90` (28×2,
  24×1) are indented differently in the parts the two share.
- `psatom.f90` documents `!uua`, which `psatom_r.f90` does not (neither has
  such an argument).
- `fpovlp_r.f90` exists but is not compiled. CMake lists `fpovlp.f90` in its
  place.
- Siblings also diverge in behaviour (`run_vkb_r` dead branch, `psatom`
  `it==400`, `run_diag_sr_so_r` error handling); see
  [Known Issues](known_issues.md).

## Fixing

Changing anything in `write` statements changes program output, and the test
suite compares `.out` files against `tests/ref/` with `fldiff`. Any fix that
touches message text needs the reference outputs regenerated or checked.

**Safe to apply mechanically** (no change to semantics or output):

- Strip trailing whitespace (`sed -i 's/[ \t]*$//' src/*.f90 src/*.F90`).
- Replace `enddo`→`end do`, `endif`→`end if`, `elseif`→`else if`, and
  `CONTAINS`→`contains`, in code lines only.
- Replace `.eq.`/`.ne.`/`.lt.`/`.le.`/`.gt.`/`.ge.` with symbolic operators
  in code (not strings). fprettify's `--enable-replacements --c-relations`
  does this alone.
- `character*N`→`character(len=N)`; add `::` to the 14 declarations without
  it.
- Normalise section-label spelling and case to the four canonical labels
  (comments only, so `sed` is safe), and fix the header typos (`vafriables`,
  `vaiables`, `OUTPuT`, `1989-201r`).
- Remove the stray first lines (`!:`, `!_der`, blank) and restore
  `run_diag_sr_so_r.f90`'s truncated copyright line.
- Name the bare `end`s in `input_toml_m.f90:792,802`.

**Needs manual work or a decision:**

- **Indentation.** fprettify cannot reproduce the column-2 base with a
  2-space step. Running it would reformat every legacy file into a different
  style. First decide whether the standard is "legacy" or "module" style; then
  fix the 29 misaligned ends and 12 unindented bodies by hand, or reformat
  everything in one dedicated commit.
- Line length over 80, if 80 is adopted as the limit.
- Section accuracy: moving the 23 misfiled locals and the 18 mismatched
  intents. This needs checking against the argument list.
- Argument documentation (233 missing) and the 23 stale argument comments.
  These need knowledge of each argument.
- Commented-out code and dead routines: each block needs a keep/delete
  decision. Some (`amesh` alternatives in `oncvpsp*.F90`) document real
  options.
- Error-message prefixes and list-directed debug prints, because they change
  output and reference files.
- One-letter names, `real(8)`/`double precision`, kind definitions and
  module naming, which interact with `fortran_best_practices.md`.
- Imported code (`exc_off_pbe.f`, `functionals.F90`,
  `gg1cc`/`gp1cc`/`gpp1cc`, `m_libxc_list.F90`): decide whether to exempt it
  (easiest to keep in sync with upstream) or convert it. Converting
  `exc_off_pbe.f` to free form includes its two tab-indented code lines.
- Splitting multi-unit files (`aeo.f90`, `exc_libxc.f90`) and renaming
  `functionals.F90`/`exc_off_pbe.f` requires `CMakeLists.txt` and `Makefile`
  updates.
- `src/CMakeLists.txt` duplicate entries (`fpovlp.f90`, `relatom.f90`) and
  the orphaned `fpovlp_r.f90`, `oncpsp_nr.f90` and `oncvpsp_nr.F90`.
