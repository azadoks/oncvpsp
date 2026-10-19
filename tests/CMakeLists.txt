enable_testing()

set(test_data "${CMAKE_SOURCE_DIR}/tests/data")
set(test_refs "${CMAKE_SOURCE_DIR}/tests/refs")

# Run TARGET with ARGS (a shell fragment, e.g. "< input.dat") and compare its
# output with REF from the "ATOM AND REFERENCE CONFIGURATION" line on. The
# output goes to tests/data/OUT.out and the comparison to tests/data/OUT.diff.
# Each test runs in its own directory, so that files written to the working
# directory (ONCVPSP.psml, _tmp_input) don't clash under ctest -j. CHECK is
# an optional shell fragment run there after a successful comparison.
function(add_ref_test name target label args out ref)
    set(check "${ARGN}")
    set(test_out "${test_data}/${out}.out")
    set(diff_out "${test_data}/${out}.diff")
    set(workdir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${workdir}")
    add_test(
        NAME ${name}
        COMMAND bash -c "
            TEMP1=${name}.tmp1
            TEMP2=${name}.tmp2
            $<TARGET_FILE:${target}> ${args} > ${test_out}
            if [ \$? -ne 0 ]; then
                echo '${target} execution failed'
                exit 1
            fi
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${ref} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            [ \$success -eq 0 ] || exit \$success
            ${check}
        "
        WORKING_DIRECTORY "${workdir}"
    )
    set_property(TEST ${name} PROPERTY LABELS "${label}")
endfunction()

# Run TARGET with ARGS and require exit status 1 and MESSAGE in the output
# (see "Error messages" in doc/coding_standards.md). SETUP is an optional
# shell fragment run first, e.g. to write a modified input file; text input
# must be a file, not a pipe, because the reader backspaces.
function(add_error_test name target args message)
    set(setup "${ARGN}")
    set(workdir "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    file(MAKE_DIRECTORY "${workdir}")
    add_test(
        NAME ${name}
        COMMAND bash -c "
            ${setup}
            out=\$($<TARGET_FILE:${target}> ${args} 2>&1)
            status=\$?
            echo \"\$out\" | tail -20
            [ \$status -eq 1 ] || { echo \"exit status \$status, expected 1\"; exit 1; }
            echo \"\$out\" | grep -qF '${message}'
        "
        WORKING_DIRECTORY "${workdir}"
    )
    set_property(TEST ${name} PROPERTY LABELS "errors")
endfunction()

set(oncvpspsr_tests
    "03_Li" "07_N" "07_N_dl" "08_O" "14_Si"
    "14_Si_st" "14_Si_UPF" "17_Cl" "19_K"
    "19_K_st" "20_Ca" "22_Ti" "29_Cu" "32_Ge"
    "34_Se" "38_Sr_lxc" "40_Zr" "56_Ba" "57_La"
    "60_Nd_GHOST" "73_Ta" "79_Au_lxc" "83_Bi"
    # Added after 4.0.1; their refs come from this code, not from 4.0.1
    "08_O_wigner"   # iexc=1 (excwig)
    "14_Si_lder"    # rxpsh set, s projector below epsh1 (run_phsft, vkbphsft)
    "40_Zr_grid"    # icmod=4 with an explicit search grid (modcore3)
)
foreach(test_name ${oncvpspsr_tests})
    add_ref_test(${test_name} oncvpspsr "sr" "< ${test_data}/${test_name}.dat"
                 ${test_name} "${test_refs}/${test_name}.out")
endforeach(test_name)

# PSML output with a model core charge (m_psmlout, icmod>=1)
add_ref_test(14_Si_psml oncvpspsr "sr" "< ${test_data}/14_Si_psml.dat"
             14_Si_psml "${test_refs}/14_Si_psml.out"
             "grep -q '<pseudocore-charge' ONCVPSP.psml")

if(ONCVPSP_WITH_TOML)
    # TOML input of 14_Si.dat must give the same output
    add_ref_test(14_Si_toml oncvpspsr "sr" "-t ${test_data}/14_Si.toml"
                 14_Si_toml "${test_refs}/14_Si.out")
endif()

set(oncvpspr_tests
    "52_Te" "74_W" "80_Hg"
)
foreach(test_name ${oncvpspr_tests})
    add_ref_test(${test_name} oncvpspr "r" "< ${test_data}/${test_name}.dat"
                 ${test_name}_r "${test_refs}/${test_name}_r.out")
endforeach(test_name)

set(oncvpspnr_tests
    "08_O" "14_Si"
)
foreach(test_name ${oncvpspnr_tests})
    add_ref_test(${test_name}_nr oncvpspnr "nr" "< ${test_data}/${test_name}.dat"
                 ${test_name}_nr "${test_refs}/${test_name}_nr.out")
endforeach(test_name)

# Error handling
foreach(target oncvpspsr oncvpspr oncvpspnr)
    if(target STREQUAL "oncvpspr")
        set(program_name "oncvpsp_r")
    else()
        set(program_name "oncvpsp")
    endif()
    add_error_test(${target}_missing_input ${target}
                   "-i ${CMAKE_CURRENT_BINARY_DIR}/no_such_file.dat"
                   "${program_name}: ERROR cannot open input file")
endforeach(target)

# rxpsh beyond the end of the radial mesh
add_error_test(rxpsh_beyond_mesh oncvpspsr "< input.dat"
    "run_phsft: ERROR rxpsh beyond radial mesh"
    "sed -e 's/^   -2.0  2.0  0.02$/   -2.0  2.0  0.02  100.0/' ${test_data}/14_Si.dat > input.dat")

# icmod=4 search grid with a zero step
add_error_test(icmod4_bad_grid oncvpspsr "< input.dat"
    "check_data: ERROR icmod=4 needs 0<fcfact_min<=fcfact_max and fcfact_step>0"
    "sed -e 's/^    4   0\\.0$/    4   0.0   0.0   2.0   5.0   0.0   1.1   1.6   0.05/' ${test_data}/40_Zr.dat > input.dat")

if(ONCVPSP_WITH_TOML)
    # TOML input without a required key
    add_error_test(toml_missing_key oncvpspsr "-t input.toml"
        "check_stat: ERROR missing or invalid value, key=iexc"
        "grep -v '^iexc' ${test_data}/14_Si.toml > input.toml")
endif()
