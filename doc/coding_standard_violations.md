# Coding standard violations

Departures from [Coding Standards](coding_standards.md) in the Fortran in
`src/` (the vendored `src/xmlf90-wxml/` is excluded). The audit of
October 2026 and the two passes that followed it (fortitude, then commits
`f5a2ab7` to `f7f76ca`) fixed every violation they found; the details are
in the commit messages. fortitude and `scripts/reindent.py` now check new
code (see [Linting](coding_standards.md#linting)). Suspected bugs are in
[Known Issues](known_issues.md).

## Open

None.

## Kept on purpose

- The imported files (`exc_off_pbe.f`, `functionals.F90`,
  `gg1cc`/`gp1cc`/`gpp1cc` and `m_libxc_list.F90`), which keep their
  upstream form so that they can be kept in sync with their sources.
- Copyright and license text: the years in `linout.f90`, `ldiracf*.f90`
  and `exc_off_pbe.f`, and the missing headers of `input_text_m.f90` and
  `input_toml_m.f90`, whose authors would have to add them.
- The text of messages from 4.0.1, including the `test_data:` prefix in
  `check_data`, the `_r` routines' non-`_r` prefixes, and the
  list-directed `lmax` and `calling psmlout` prints, which the reference
  outputs contain (see [Error messages](coding_standards.md#error-messages)).
- `derivs` in `exc_libxc.f90`, an external routine that `functionals.F90`
  calls too, so it can't become internal, and a file of its own would need
  a change to every user's `make.inc` (`OBJS_LIBXC`).
- Derived types from the libraries and imported files: toml-f's in
  `input_toml_m`, and xmlf90's `xmlf_t` and `m_libxc_list`'s `libxc_t`
  in `psmlout_m` (see rule 3 under
  [Conventions](coding_standards.md#conventions)).
- Unused locals (e.g. `maxset` in `lschpsbar`, `cfloc` in `check_data`),
  which fortitude's default rules don't flag.
