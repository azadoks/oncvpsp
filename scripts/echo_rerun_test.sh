#!/bin/bash
# Checks the input echoes for one input: with psfile set to FORMAT (psp8, upf
# or psml), the input echoed in the output, and the input in the
# pseudopotential file, run again, must each give the same output
# (compare_output.py --identical: every number to 1e-10, the text exactly).
#
# usage: echo_rerun_test.sh EXE INPUT.dat FORMAT WORKDIR [SED]
#   SED  a sed script applied to the input first (e.g. to set an optional
#        value)
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
exe=$1; dat=$2; format=$3; work=$4; edit=${5:-}
rm -rf "$work"
mkdir -p "$work"
cd "$work"

# psfile is the sixth field of the first data line
sed -e "$edit" "$dat" |
    awk -v f="$format" '!/^#/ && !done {$6=f; done=1} {print}' > in.dat
"$exe" < in.dat > in.out

# rerun NAME: run NAME.dat, which must hold an input, and compare with in.out
rerun() {
    if ! grep -q '^# ATOM AND REFERENCE CONFIGURATION$' "$1.dat"; then
        echo "echo_rerun_test.sh: no input in $1.dat" >&2; exit 1
    fi
    "$exe" < "$1.dat" > "$1.out"
    same in.out "$1.out"
}

awk '/^# ATOM AND REFERENCE CONFIGURATION$/ {on=1}
     /^Reference configufation results/ {exit}
     on' in.out > output_echo.dat
rerun output_echo

case $format in
    psp8)
        awk '/^<\/INPUT>$/ {on=0} on; /^<INPUT>$/ {on=1}' in.out > file_echo.dat ;;
    upf)
        awk '/<\/PP_INPUTFILE>/ {on=0} on; /<PP_INPUTFILE>/ {on=1}' in.out > file_echo.dat ;;
    psml)
        awk '/^Begin PSML$/ {f=1; next} f && /^END_PSP$/ {exit} f' in.out |
            awk '/<input-file/ {on=1; sub(/.*<input-file[^>]*>/, "")}
                 /<\/input-file>/ {sub(/<\/input-file>.*/, ""); print; on=0; next}
                 on' > file_echo.dat ;;
    *)
        echo "echo_rerun_test.sh: unknown format $format" >&2; exit 1 ;;
esac
rerun file_echo
echo "the input echoes in the output and the $format file of $(basename "$dat") give the same output"
