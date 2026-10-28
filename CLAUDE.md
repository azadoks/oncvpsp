# oncvpsp (smooooth)

## Docs that move with the code

A change is done when the documents describing it are updated in the same
commit. Update every one the change touches:

- **`doc/changes_from_4.0.1.md`**: the running changelog against 4.0.1.
  Any change that moves an output, or changes what an input means, adds
  or edits an entry under Methods or Numerics (the page's preamble says
  which), giving the 4.0.1 behaviour and the new one. Its commit hash
  goes in with the next commit that touches the page.
- **Input documentation**, for a new or changed input value:
  `doc/32_Ge_annotated.dat`, the TOML section of `doc/users_guide.md`,
  `scripts/input_text_to_toml.py`, and the four input echoes (see "Input
  echo" in `doc/coding_standards.md`).
- **`--help` text** of every executable whose options change, and the
  README where it lists them.
