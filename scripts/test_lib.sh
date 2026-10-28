# shellcheck shell=bash
# Sourced by the test scripts (*_test.sh, run_test.sh), first thing:
#   source "$(dirname "$0")/test_lib.sh"
# It sets the strict mode, scripts (this directory), python (the Python of
# the comparisons, $ONCVPSP_PYTHON or python3) and the helpers below.
set -euo pipefail
# shellcheck disable=SC2034 # used by the scripts that source this one
scripts=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
# shellcheck disable=SC2034
python=${ONCVPSP_PYTHON:-python3}

# fail MESSAGE: "<script>: ERROR MESSAGE" on stderr, and exit 1
fail() {
    echo "$(basename "$0"): ERROR $*" >&2
    exit 1
}

# same A.out B.out [OPTION...]: B the same output as A by compare_output.py
# --identical (every number to 1e-10, the text exactly) with the OPTIONs, or
# the report (left in B.diff) and exit 1
same() {
    if ! "$python" "$scripts/compare_output.py" --identical "${@:3}" "$1" "$2" \
            > "${2%.out}.diff"; then
        cat "${2%.out}.diff" >&2
        exit 1
    fi
}
