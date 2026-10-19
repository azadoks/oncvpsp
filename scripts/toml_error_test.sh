#!/bin/bash
# Checks how the TOML reader handles one change to a valid input: it must
# either stop with a non-zero status and the expected message on stdout,
# or, for a change that does not alter the input data, give the output of
# the unchanged input.
#
# usage: toml_error_test.sh EXE DATA_DIR CASE WORKDIR
# DATA_DIR is tests/data; the inputs are DATA_DIR/toml/14_Si.toml and, for
# icmod 6, DATA_DIR/toml/14_Si_bca.toml
set -euo pipefail
exe=$1; data=$2; case=$3; work=$4
mkdir -p "$work"
cd "$work"

# set_key SECTION KEY LINE < IN > OUT: replace the first KEY = ... line of
# [SECTION] (or of the first [[SECTION]]) with LINE, or delete it if LINE
# is empty
set_key() {
    awk -v sec="$1" -v key="$2" -v new="$3" '
        /^\[/ {cur=$0; gsub(/[][ ]/, "", cur)}
        !done && cur == sec && $1 == key {done=1; if (new != "") print new; next}
        {print}
        END {if (!done) {print "set_key: no " key " in [" sec "]" > "/dev/stderr"; exit 1}}'
}

# drop_table SECTION < IN > OUT: delete [SECTION], or every [[SECTION]]
drop_table() {
    awk -v sec="$1" '
        /^\[/ {cur=$0; gsub(/[][ ]/, "", cur)}
        cur != sec {print}'
}

base=$data/toml/14_Si.toml
expect=same
case $case in
    parse_error)
        { cat "$base"; echo 'z = = 14'; } > in.toml
        expect='cannot parse the input' ;;
    missing_table)
        drop_table linear_mesh < "$base" > in.toml
        expect='missing table [linear_mesh]' ;;
    not_a_table)
        { echo 'oncvpsp = 3'; drop_table oncvpsp < "$base"; } > in.toml
        expect='[oncvpsp] must be a table' ;;
    missing_key)
        set_key oncvpsp z '' < "$base" > in.toml
        expect='missing key z in [oncvpsp]' ;;
    missing_string)
        set_key pp_output psfile '' < "$base" > in.toml
        expect='missing key psfile in [pp_output]' ;;
    integer_type)
        set_key oncvpsp nc 'nc = 3.0' < "$base" > in.toml
        expect='nc in [oncvpsp] must be an integer' ;;
    number_type)
        set_key oncvpsp z 'z = "Si"' < "$base" > in.toml
        expect='z in [oncvpsp] must be a number' ;;
    string_type)
        set_key pp_output psfile 'psfile = 8' < "$base" > in.toml
        expect='psfile in [pp_output] must be a string' ;;
    atsym_length)
        set_key oncvpsp atsym 'atsym = "Sil"' < "$base" > in.toml
        expect='atsym in [oncvpsp] must have 1 or 2 characters' ;;
    state_counts)
        set_key oncvpsp nc 'nc = 29' < "$base" > in.toml
        expect='nc and nv in [oncvpsp] must be >= 0 with nc + nv <= 30' ;;
    lmax_range)
        set_key pseudopotentials lmax 'lmax = 6' < "$base" > in.toml
        expect='lmax in [pseudopotentials] must be in 0..5' ;;
    array_length)
        set_key pseudopotentials rc 'rc = [ 2.7, 2.7,]' < "$base" > in.toml
        expect='rc in [pseudopotentials] has 2 instead of 3 elements' ;;
    array_type)
        set_key reference_configuration f 'f = [ 2.0, "x", 6.0, 2.0, 2.0,]' < "$base" > in.toml
        expect='f in [reference_configuration] must hold numbers' ;;
    not_an_array)
        set_key reference_configuration n 'n = 1' < "$base" > in.toml
        expect='n in [reference_configuration] must be an array' ;;
    l_order)
        set_key vkb_projectors l 'l = [ 0, 2, 1,]' < "$base" > in.toml
        expect='l in [vkb_projectors] must be 0, 1, ..., lmax in order' ;;
    icmod3_rcfact)
        set_key model_core_charge icmod 'icmod = 3' < "$base" > in.toml
        expect='missing key rcfact in [model_core_charge]' ;;
    test_config_length)
        set_key test_configurations l 'l = [ 0,]' < "$base" > in.toml
        expect='l in [test_configurations] has 1 instead of 2 elements' ;;
    test_config_count)
        { cat "$base"; printf '\n[[test_configurations]]\nn = [ 3, 3,]\nl = [ 0, 1,]\nf = [ 2.0, 2.0,]\n'; \
          printf '\n[[test_configurations]]\nn = [ 3, 3,]\nl = [ 0, 1,]\nf = [ 2.0, 2.0,]\n'; } > in.toml
        expect='at most 4 [[test_configurations]] are allowed' ;;
    test_config_type)
        { echo 'test_configurations = 3'; drop_table test_configurations < "$base"; } > in.toml
        expect='test_configurations must be an array of tables' ;;
    integer_for_real)
        # an integer where a real is read is the same number
        set_key oncvpsp z 'z = 14' < "$base" > in.toml ;;
    extra_keys)
        # keys and tables the reader does not know are ignored
        { echo 'comment = "extra"'; set_key oncvpsp z 'z = 14.0\nextra = 1' < "$base"; \
          printf '\n[extra_table]\nx = 1\n'; } > in.toml ;;
    icmod6_default_rcfact)
        # rcfact defaults to 1e-5 for icmod 6
        set_key model_core_charge rcfact 'rcfact = 1.0e-5' \
            < "$data/toml/14_Si_bca.toml" > base.toml
        base=base.toml
        set_key model_core_charge rcfact '' < "$base" > in.toml ;;
    icmod6_default_slack)
        # slack defaults to 1 for icmod 6
        base=$data/toml/14_Si_bca.toml
        set_key model_core_charge rcfact 'rcfact = 1.0e-4\nslack = 1.0' \
            < "$base" > in.toml ;;
    *)
        echo "toml_error_test.sh: unknown case $case" >&2
        exit 2 ;;
esac

if [ "$expect" = same ]; then
    "$exe" -t "$base" > base.out
    "$exe" -t in.toml > in.out
    cmp base.out in.out
    echo "$case: same output as $(basename "$base")"
else
    status=0
    "$exe" -t in.toml > in.out 2> in.err || status=$?
    if [ $status -eq 0 ]; then
        echo "$case: the program did not stop" >&2
        exit 1
    fi
    if ! grep -qF "TOML input ERROR: $expect" in.out; then
        echo "$case: expected 'TOML input ERROR: $expect' in the output, got:" >&2
        tail -5 in.out >&2
        cat in.err >&2
        exit 1
    fi
    echo "$case: stopped with status $status and '$expect'"
fi
