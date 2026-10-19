# The tests run the scripts in scripts/, which need bash and (for the comparison with the references)
# perl. Each test works in its own directory under this one, where its outputs are left.
find_program(ONCVPSP_BASH bash)
find_package(Perl)
if (NOT ONCVPSP_BASH OR NOT PERL_FOUND)
    message(WARNING "bash or perl not found; the tests are not available")
    return()
endif ()

# Unit tests of the input readers
add_executable(test_input_text unit/test_input_text.f90)
target_link_libraries(test_input_text PRIVATE oncvpsp::oncvpsp oncvpsp::compile_options)
set_target_properties(test_input_text PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME unit_input_text COMMAND test_input_text)
set_tests_properties(unit_input_text PROPERTIES LABELS unit)

set(scripts "${PROJECT_SOURCE_DIR}/scripts")
set(data "${CMAKE_CURRENT_SOURCE_DIR}/data")

# The inputs with Libxc functionals (iexc < 0), whose tests are disabled without Libxc
set(libxc_inputs "38_Sr_lxc" "60_Nd_GHOST" "79_Au_lxc" "83_Bi")

# oncvpsp_test(NAME LABEL SCRIPT ARGS...): run SCRIPT with ARGS, labelled LABEL
function(oncvpsp_test name label script)
    add_test(NAME ${name} COMMAND ${ONCVPSP_BASH} "${scripts}/${script}" ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS ${label})
    if (NOT ONCVPSP_WITH_Libxc)
        foreach (input IN LISTS libxc_inputs)
            if (name MATCHES "(^|_)${input}$")
                set_tests_properties(${name} PROPERTIES DISABLED TRUE)
            endif ()
        endforeach ()
    endif ()
endfunction()

# Output against the references in tests/refs
set(oncvpspsr_tests
    "03_Li" "07_N" "07_N_dl" "08_O" "14_Si"
    "14_Si_bc" "14_Si_bca" "14_Si_st" "14_Si_UPF" "17_Cl" "19_K"
    "19_K_st" "20_Ca" "22_Ti" "29_Cu" "32_Ge"
    "34_Se" "38_Sr_lxc" "40_Zr" "56_Ba" "57_La"
    "60_Nd_GHOST" "73_Ta" "79_Au_lxc" "83_Bi"
    "08_O_wigner"   # iexc=1 (excwig)
    "14_Si_lder"    # rxpsh set, s projector below epsh1 (run_phsft, vkbphsft)
    "40_Zr_grid"    # icmod=4 with an explicit search grid (modcore3)
)
foreach (test_name IN LISTS oncvpspsr_tests)
    oncvpsp_test(${test_name} sr ref_test.sh
        $<TARGET_FILE:oncvpspsr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
endforeach ()

set(oncvpspr_tests "52_Te" "74_W" "80_Hg")
foreach (test_name IN LISTS oncvpspr_tests)
    oncvpsp_test(${test_name} r ref_test.sh
        $<TARGET_FILE:oncvpspr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}_r.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}_r")
endforeach ()

set(oncvpspnr_tests "08_O" "14_Si")
foreach (test_name IN LISTS oncvpspnr_tests)
    oncvpsp_test(${test_name}_nr nr ref_test.sh
        $<TARGET_FILE:oncvpspnr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}_nr.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}_nr")
endforeach ()

# All-electron / pseudopotential split: oncvpsp_ae.x writes the all-electron atom to a file and
# oncvpsp_ps.x reads it (see scripts/split_test.sh)
foreach (split_test IN ITEMS "sr:14_Si" "sr:19_K_st" "sr:29_Cu" "sr:57_La" "nr:14_Si" "nr:29_Cu"
                             "r:52_Te" "r:80_Hg")
    string(REPLACE ":" ";" split_test "${split_test}")
    list(GET split_test 0 rel)
    list(GET split_test 1 test_name)
    oncvpsp_test(split_${rel}_${test_name} split split_test.sh
        $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp${rel}_ps>
        "${data}/${test_name}.dat" "${CMAKE_CURRENT_BINARY_DIR}/split_${rel}_${test_name}")
endforeach ()

# Damaged or mismatched all-electron data files: oncvpsp_ps.x must stop with a message naming the
# problem, but accept differences within round-off (see scripts/ae_file_error_test.sh)
foreach (case IN ITEMS "byte_order" "zz" "fa" "mesh" "zz_roundoff" "fa_roundoff" "mesh_roundoff")
    oncvpsp_test(ae_file_${case} split ae_file_error_test.sh
        $<TARGET_FILE:oncvpspsr> $<TARGET_FILE:oncvpspsr_ae> $<TARGET_FILE:oncvpspsr_ps>
        "${data}/14_Si.dat" ${case} "${CMAKE_CURRENT_BINARY_DIR}/ae_file_${case}")
endforeach ()

# Log mesh input (the optional last line of the text input): see scripts/log_mesh_test.sh
foreach (mesh_test IN ITEMS "sr:14_Si" "nr:29_Cu" "r:52_Te")
    string(REPLACE ":" ";" mesh_test "${mesh_test}")
    list(GET mesh_test 0 rel)
    list(GET mesh_test 1 test_name)
    set(mesh_cases "defaults" "explicit" "rerun" "upf" "split" "split_differ")
    if (rel STREQUAL "sr")
        list(APPEND mesh_cases "amesh" "rr1" "rmax" "read")
    endif ()
    if (ONCVPSP_WITH_TOML)
        list(APPEND mesh_cases "toml")
    endif ()
    foreach (case IN LISTS mesh_cases)
        oncvpsp_test(log_mesh_${rel}_${test_name}_${case} mesh log_mesh_test.sh
            $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp${rel}_ps>
            "${data}/${test_name}.dat" ${case} "${CMAKE_CURRENT_BINARY_DIR}/log_mesh_${rel}_${test_name}_${case}"
            "${data}/toml/${test_name}.toml")
    endforeach ()
endforeach ()

# TOML input: every input in tests/data/toml (the TOML form of its tests/data/*.dat) must give the
# output of the text input, and the reader must stop with a message naming the problem for each of the
# changes in scripts/toml_error_test.sh
if (ONCVPSP_WITH_TOML)
    set(toml_tests ${oncvpspsr_tests})
    list(TRANSFORM toml_tests PREPEND "sr:")
    list(APPEND toml_tests "r:52_Te" "r:74_W" "r:80_Hg" "nr:14_Si" "nr:29_Cu")
    foreach (toml_test IN LISTS toml_tests)
        string(REPLACE ":" ";" toml_test "${toml_test}")
        list(GET toml_test 0 rel)
        list(GET toml_test 1 test_name)
        oncvpsp_test(toml_${rel}_${test_name} toml toml_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${data}/toml/${test_name}.toml"
            "${CMAKE_CURRENT_BINARY_DIR}/toml_${rel}_${test_name}")
    endforeach ()

    foreach (case IN ITEMS
            "parse_error" "missing_table" "not_a_table" "missing_key" "missing_string"
            "integer_type" "number_type" "string_type" "atsym_length" "state_counts"
            "lmax_range" "array_length" "array_type" "not_an_array" "l_order"
            "icmod3_rcfact" "test_config_length" "test_config_count" "test_config_type"
            "integer_for_real" "extra_keys" "icmod6_default_rcfact")
        oncvpsp_test(toml_input_${case} toml toml_error_test.sh
            $<TARGET_FILE:oncvpspsr> "${data}" ${case} "${CMAKE_CURRENT_BINARY_DIR}/toml_input_${case}")
    endforeach ()
endif ()

# PSML output: written into the output like the psp8 and upf files, with no files left in the working
# directory
foreach (psml_test IN ITEMS "sr:14_Si" "r:52_Te" "nr:08_O")
    string(REPLACE ":" ";" psml_test "${psml_test}")
    list(GET psml_test 0 rel)
    list(GET psml_test 1 test_name)
    oncvpsp_test(psml_${rel}_${test_name} psml psml_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${CMAKE_CURRENT_BINARY_DIR}/psml_${rel}_${test_name}")
endforeach ()

# Input echoes: the input echoed in the output, and in the psp8, UPF or PSML file, run again, give the
# same output (see scripts/echo_rerun_test.sh). The inputs cover rxpsh (14_Si_lder, and 52_Te with
# rxpsh set), icmod=4 with its search grid (40_Zr_grid, 74_W), a Libxc iexc (83_Bi) and rc on the
# grid (03_Li).
set(te_rxpsh "s/^   -3\\.0  3\\.0  0\\.02$/   -3.0  3.0  0.02  3.5/")
foreach (echo_test IN ITEMS
        "sr:14_Si_lder:psp8" "sr:14_Si_lder:upf" "sr:14_Si_lder:psml"
        "sr:40_Zr_grid:psp8" "sr:40_Zr_grid:upf" "sr:40_Zr_grid:psml"
        "sr:83_Bi:psml" "sr:03_Li:psml" "nr:14_Si:psml"
        "r:74_W:psp8" "r:74_W:psml" "r:52_Te:psp8" "r:52_Te:upf" "r:52_Te:psml")
    string(REPLACE ":" ";" echo_test "${echo_test}")
    list(GET echo_test 0 rel)
    list(GET echo_test 1 test_name)
    list(GET echo_test 2 format)
    set(edit "")
    if (test_name STREQUAL "52_Te")
        set(edit "${te_rxpsh}")
    endif ()
    oncvpsp_test(echo_${format}_${rel}_${test_name} echo echo_rerun_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" ${format}
        "${CMAKE_CURRENT_BINARY_DIR}/echo_${format}_${rel}_${test_name}" "${edit}")
endforeach ()

# icmod=6 with fcfact<=0 (the cutoff taken from the pseudo-waves): the psp8 file must still carry the
# model core (see scripts/psp8_core_test.sh)
oncvpsp_test(psp8_core_sr_14_Si psp8 psp8_core_test.sh
    $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat" "s/^   1    0\\.25$/   6    0.0   1.0d-3/"
    "${CMAKE_CURRENT_BINARY_DIR}/psp8_core_sr_14_Si")
oncvpsp_test(psp8_core_r_74_W psp8 psp8_core_test.sh
    $<TARGET_FILE:oncvpspr> "${data}/74_W.dat" "s/^    4  0\\.0$/    6  0.0  1.0d-3/"
    "${CMAKE_CURRENT_BINARY_DIR}/psp8_core_r_74_W")

# Text input from a pipe (as from a Python subprocess), which cannot be backspaced: the same output
# as from a file, and an input ending early reported with the line the reader stopped at
foreach (pipe_test IN ITEMS "sr:14_Si" "r:52_Te" "nr:08_O")
    string(REPLACE ":" ";" pipe_test "${pipe_test}")
    list(GET pipe_test 0 rel)
    list(GET pipe_test 1 test_name)
    oncvpsp_test(pipe_${rel}_${test_name} pipe pipe_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${CMAKE_CURRENT_BINARY_DIR}/pipe_${rel}_${test_name}")
endforeach ()
oncvpsp_test(pipe_truncated pipe pipe_truncated_test.sh $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat")

# Input errors: exit status 1 and the documented message (see scripts/error_test.sh)
foreach (rel IN ITEMS "sr" "nr" "r")
    if (rel STREQUAL "r")
        set(program_name "oncvpsp_r")
    else ()
        set(program_name "oncvpsp")
    endif ()
    oncvpsp_test(error_${rel}_missing_input errors error_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${program_name}: ERROR cannot open input file"
        "${CMAKE_CURRENT_BINARY_DIR}/error_${rel}_missing_input")
endforeach ()
oncvpsp_test(error_rxpsh_beyond_mesh errors error_test.sh
    $<TARGET_FILE:oncvpspsr> "run_phsft: ERROR rxpsh beyond radial mesh"
    "${CMAKE_CURRENT_BINARY_DIR}/error_rxpsh_beyond_mesh"
    "${data}/14_Si.dat" "s/^   -2.0  2.0  0.02$/   -2.0  2.0  0.02  100.0/")
oncvpsp_test(error_icmod4_bad_grid errors error_test.sh
    $<TARGET_FILE:oncvpspsr> "check_data: ERROR icmod=4 needs 0<fcfact_min<=fcfact_max and fcfact_step>0"
    "${CMAKE_CURRENT_BINARY_DIR}/error_icmod4_bad_grid"
    "${data}/40_Zr.dat" "s/^    4   0\\.0$/    4   0.0   0.0   2.0   5.0   0.0   1.1   1.6   0.05/")
