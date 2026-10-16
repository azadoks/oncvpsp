#!/bin/bash
# Checks that the text input can be read from a pipe (as from a Python
# subprocess), which the program cannot backspace, as well as from a file:
# both must give the same output, apart from lines naming a date.
#
# usage: pipe_test.sh EXE INPUT.dat WORKDIR
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); work=$3
rm -rf "$work"
mkdir -p "$work"
cd "$work"
"$exe" < "$dat" > file.out
cat "$dat" | "$exe" > pipe.out
diff <(sed '/date/d' file.out) <(sed '/date/d' pipe.out)
grep -q '^# atsym.*psfile' pipe.out
echo "$(basename "$exe") gives the same output for $(basename "$dat") from a pipe as from a file"
