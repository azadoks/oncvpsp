enable_testing()

set(oncvpspsr_tests
    "03_Li" "07_N" "07_N_dl" "08_O" "14_Si"
    "14_Si_bc" "14_Si_bca" "14_Si_st" "14_Si_UPF" "17_Cl" "19_K"
    "19_K_st" "20_Ca" "22_Ti" "29_Cu" "32_Ge"
    "34_Se" "38_Sr_lxc" "40_Zr" "56_Ba" "57_La"
    "60_Nd_GHOST" "73_Ta" "79_Au_lxc" "83_Bi"
)
foreach(test_name ${oncvpspsr_tests})
    set(test_dat "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat")
    set(test_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.out")
    set(ref_out "${CMAKE_SOURCE_DIR}/tests/refs/${test_name}.out")
    set(diff_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.diff")
    add_test(
        NAME ${test_name}
        COMMAND bash -c "
            TEMP1=${test_name}.tmp1
            TEMP2=${test_name}.tmp2
            $<TARGET_FILE:oncvpspsr> < ${test_dat} > ${test_out}
            if [ \$? -ne 0 ]; then
                echo 'oncvpspsr execution failed'
                exit 1
            fi
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${ref_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | \
                sed -e '/^Begin PSML/,/^END_PSP/s/^/-/' >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | \
                sed -e '/^Begin PSML/,/^END_PSP/s/^/-/' >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            exit \$success
        "
    )
    set_property(TEST ${test_name} PROPERTY LABEL "sr")
endforeach(test_name)

set(oncvpspr_tests
    "52_Te" "74_W" "80_Hg"
)
foreach(test_name ${oncvpspr_tests})
    set(test_dat "${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat")
    set(test_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}_r.out")
    set(ref_out "${CMAKE_SOURCE_DIR}/tests/refs/${test_name}_r.out")
    set(diff_out "${CMAKE_SOURCE_DIR}/tests/data/${test_name}_r.diff")
    add_test(
        NAME ${test_name}
        COMMAND bash -c "
            TEMP1=${test_name}.tmp1
            TEMP2=${test_name}.tmp2
            $<TARGET_FILE:oncvpspr> < ${test_dat} > ${test_out}
            if [ \$? -ne 0 ]; then
                echo 'oncvpspr execution failed'
                exit 1
            fi
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${ref_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | \
                sed -e '/^Begin PSML/,/^END_PSP/s/^/-/' >$TEMP1
            awk 'BEGIN{out=0}; {if(out == 1) {print}};\
                /ATOM AND REFERENCE CONFIGURATION/{out=1}' ${test_out} | \
                sed -e /pspd/s/^/-/ | sed -e /date/s/^/-/ | \
                sed -e '/^Begin PSML/,/^END_PSP/s/^/-/' >$TEMP2
            perl ${CMAKE_SOURCE_DIR}/scripts/fldiff.pl -easy $TEMP1 $TEMP2 > ${diff_out}
            rm -f $TEMP1 $TEMP2
            grep -E 'Summary  : no significant difference|Summary  : different lines=\\s*[0-9]+ , max abs_diff= -?0\\.000e[\\+\\-]00 \\(.*\\), max rel_diff= -?0\\.000e[\\+\\-]00 \\(.*\\)' ${diff_out}
            success=\$?
            cat ${diff_out}
            exit \$success
        "
    )
    set_property(TEST ${test_name} PROPERTY LABEL "r")
endforeach(test_name)

# All-electron / pseudopotential split: oncvpsp_ae.x writes the
# all-electron atom to a file and oncvpsp_ps.x reads it (see
# scripts/split_test.sh)
set(oncvpsp_split_tests
    "sr:14_Si" "sr:19_K_st" "sr:29_Cu" "sr:57_La" "nr:14_Si" "nr:29_Cu"
    "r:52_Te" "r:80_Hg"
)
foreach(split_test ${oncvpsp_split_tests})
    string(REPLACE ":" ";" split_test_parts ${split_test})
    list(GET split_test_parts 0 rel)
    list(GET split_test_parts 1 test_name)
    add_test(
        NAME split_${rel}_${test_name}
        COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/split_test.sh
            $<TARGET_FILE:oncvpsp${rel}>
            $<TARGET_FILE:oncvpsp${rel}_ae>
            $<TARGET_FILE:oncvpsp${rel}_ps>
            ${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat
            ${CMAKE_CURRENT_BINARY_DIR}/split_${rel}_${test_name}
    )
    set_property(TEST split_${rel}_${test_name} PROPERTY LABEL "split")
endforeach(split_test)

# TOML input: every input in tests/data/toml (the TOML form of its
# tests/data/*.dat) must give the output of the text input, and the reader
# must stop with a message naming the problem for each of the changes in
# scripts/toml_error_test.sh
if (ONCVPSP_WITH_TOML)
    set(oncvpsp_toml_tests ${oncvpspsr_tests})
    list(TRANSFORM oncvpsp_toml_tests PREPEND "sr:")
    list(APPEND oncvpsp_toml_tests "r:52_Te" "r:74_W" "r:80_Hg" "nr:14_Si" "nr:29_Cu")
    foreach(toml_test ${oncvpsp_toml_tests})
        string(REPLACE ":" ";" toml_test_parts ${toml_test})
        list(GET toml_test_parts 0 rel)
        list(GET toml_test_parts 1 test_name)
        add_test(
            NAME toml_${rel}_${test_name}
            COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/toml_test.sh
                $<TARGET_FILE:oncvpsp${rel}>
                ${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat
                ${CMAKE_SOURCE_DIR}/tests/data/toml/${test_name}.toml
                ${CMAKE_CURRENT_BINARY_DIR}/toml_${rel}_${test_name}
        )
        set_property(TEST toml_${rel}_${test_name} PROPERTY LABEL "toml")
    endforeach(toml_test)

    set(oncvpsp_toml_error_tests
        "parse_error" "missing_table" "not_a_table" "missing_key" "missing_string"
        "integer_type" "number_type" "string_type" "atsym_length" "state_counts"
        "lmax_range" "array_length" "array_type" "not_an_array" "l_order"
        "icmod3_rcfact" "test_config_length" "test_config_count" "test_config_type"
        "integer_for_real" "extra_keys" "icmod6_default_rcfact"
    )
    foreach(case ${oncvpsp_toml_error_tests})
        add_test(
            NAME toml_input_${case}
            COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/toml_error_test.sh
                $<TARGET_FILE:oncvpspsr>
                ${CMAKE_SOURCE_DIR}/tests/data
                ${case}
                ${CMAKE_CURRENT_BINARY_DIR}/toml_input_${case}
        )
        set_property(TEST toml_input_${case} PROPERTY LABEL "toml")
    endforeach(case)
endif ()

# PSML output: written into the output like the psp8 and upf files, with
# no files left in the working directory
foreach(psml_test "sr:14_Si" "r:52_Te" "nr:08_O")
    string(REPLACE ":" ";" psml_test_parts ${psml_test})
    list(GET psml_test_parts 0 rel)
    list(GET psml_test_parts 1 test_name)
    add_test(
        NAME psml_${rel}_${test_name}
        COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/psml_test.sh
            $<TARGET_FILE:oncvpsp${rel}>
            ${CMAKE_SOURCE_DIR}/tests/data/${test_name}.dat
            ${CMAKE_CURRENT_BINARY_DIR}/psml_${rel}_${test_name}
    )
    set_property(TEST psml_${rel}_${test_name} PROPERTY LABEL "psml")
endforeach(psml_test)
