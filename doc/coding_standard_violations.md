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

- The imported files (`exc_off_pbe_m.f`, `gg1cc`/`gp1cc`/`gpp1cc` and
  `m_libxc_list.F90`), which keep their
  upstream form so that they can be kept in sync with their sources.
- Copyright and license text: the years in `linout_m.f90`, `ldiracf*_m.f90`
  and `exc_off_pbe_m.f`, and the missing headers of `input_text_m.f90` and
  `input_toml_m.f90`, whose authors would have to add them.
- The text of messages from 4.0.1, including the `test_data:` prefix in
  `check_data`, the prefixes of the fully-relativistic messages (the
  names of the former `_r` routines, or of the sr ones), and the
  list-directed `lmax` and `calling psmlout` prints, which the reference
  outputs contain (see [Error messages](coding_standards.md#error-messages)).
- Unused locals (e.g. `maxset` in `lschpsbar`, `cfloc` in `check_data`),
  which fortitude's default rules don't flag.
