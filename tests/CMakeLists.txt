# The tests run the scripts in scripts/, which need bash and python3, and those that compare outputs
# (scripts/compare_output.py) a Python with pyoncvpsp and upf-tools (cmake/OncvpspTestPython.cmake),
# without which they are disabled. Each test works in its own directory under this one, where its
# outputs are left.
find_program(ONCVPSP_BASH bash)
find_package(Python3 COMPONENTS Interpreter)
if (NOT ONCVPSP_BASH OR NOT Python3_Interpreter_FOUND)
    message(WARNING "bash or python3 not found; the tests are not available")
    return()
endif ()
include(OncvpspTestPython)
set(ONCVPSP_TEST_PYTHON_DESCRIPTION "${ONCVPSP_TEST_PYTHON_DESCRIPTION}" PARENT_SCOPE)
if (ONCVPSP_COMPARE_PYTHON)
    set(test_python "${ONCVPSP_COMPARE_PYTHON}")
else ()
    set(test_python "${Python3_EXECUTABLE}")
endif ()

# Unit tests of the input readers
add_executable(test_input_text unit/test_input_text.f90)
target_link_libraries(test_input_text PRIVATE oncvpsp::oncvpsp oncvpsp::compile_options)
set_target_properties(test_input_text PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME unit_input_text COMMAND test_input_text)
set_tests_properties(unit_input_text PROPERTIES LABELS unit)

# Unit tests of the rxpsh warning in run_phsft (nr and fr), called directly
add_executable(test_run_phsft unit/test_run_phsft.f90)
target_link_libraries(test_run_phsft PRIVATE oncvpsp::oncvpsp oncvpsp::compile_options)
set_target_properties(test_run_phsft PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME unit_run_phsft COMMAND test_run_phsft WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set_tests_properties(unit_run_phsft PROPERTIES LABELS unit)

set(scripts "${PROJECT_SOURCE_DIR}/scripts")
set(data "${CMAKE_CURRENT_SOURCE_DIR}/data")

# The inputs with Libxc functionals (iexc < 0), whose tests are disabled without Libxc
set(libxc_inputs "38_Sr_lxc" "60_Nd_GHOST" "79_Au_lxc" "83_Bi")

# oncvpsp_test(NAME LABEL SCRIPT ARGS...): run SCRIPT with ARGS, labelled LABEL; not added when it
# runs a fully-relativistic executable and those are not built (ONCVPSP_WITH_FR off)
function(oncvpsp_test name label script)
    if (NOT ONCVPSP_WITH_FR AND "${ARGN}" MATCHES "TARGET_FILE:oncvpspr[_>]")
        return()
    endif ()
    add_test(NAME ${name} COMMAND ${ONCVPSP_BASH} "${scripts}/${script}" ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS ${label} ENVIRONMENT "ONCVPSP_PYTHON=${test_python}")
    if (NOT ONCVPSP_WITH_Libxc)
        foreach (input IN LISTS libxc_inputs)
            if (name MATCHES "(^|_)${input}(_r|_nr)?$")
                set_tests_properties(${name} PROPERTIES DISABLED TRUE)
            endif ()
        endforeach ()
    endif ()
endfunction()

# oncvpsp_compares(NAME): the test NAME (if it was added) compares outputs with
# scripts/compare_output.py, or parses them with pyoncvpsp: disabled without its Python
function(oncvpsp_compares name)
    if (TEST ${name} AND NOT ONCVPSP_COMPARE_PYTHON)
        set_tests_properties(${name} PROPERTIES DISABLED TRUE)
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
    "14_Si_hm"      # icmod=7 with slack (modcore_budget, modcore_probes)
)
foreach (test_name IN LISTS oncvpspsr_tests)
    oncvpsp_test(${test_name} sr ref_test.sh
        $<TARGET_FILE:oncvpspsr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
    oncvpsp_compares(${test_name})
endforeach ()

set(oncvpspr_tests "52_Te" "74_W" "80_Hg"
    "03_Li"         # no core states
    "07_N_dl"       # lloc<4
    "08_O_wigner"   # iexc=1
    "14_Si"         # icmod=1
    "14_Si_bc"      # icmod=6
    "14_Si_hm"      # icmod=7 with slack
    "14_Si_lder"    # rxpsh set
    "20_Ca"         # lpopt=1
    "34_Se"         # icmod=3
    "38_Sr_lxc"     # a single projector for one l
    "40_Zr_grid"    # icmod=4 with an explicit search grid
    "60_Nd_GHOST"   # ghosts
    "73_Ta"         # three projectors
    "79_Au_lxc"     # icmod=2
)
foreach (test_name IN LISTS oncvpspr_tests)
    oncvpsp_test(${test_name}_r r ref_test.sh
        $<TARGET_FILE:oncvpspr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}_r.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}_r")
    oncvpsp_compares(${test_name}_r)
endforeach ()

set(oncvpspnr_tests "08_O" "14_Si" "03_Li" "07_N_dl" "14_Si_lder")
foreach (test_name IN LISTS oncvpspnr_tests)
    oncvpsp_test(${test_name}_nr nr ref_test.sh
        $<TARGET_FILE:oncvpspnr> "${data}/${test_name}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${test_name}_nr.out"
        "${CMAKE_CURRENT_BINARY_DIR}/${test_name}_nr")
    oncvpsp_compares(${test_name}_nr)
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
    oncvpsp_compares(split_${rel}_${test_name})
endforeach ()

# Damaged or mismatched all-electron data files: oncvpsp_ps.x must stop with a message naming the
# problem, but accept differences within round-off (see scripts/ae_file_error_test.sh)
foreach (case IN ITEMS "byte_order" "zz" "fa" "mesh" "c_light" "zz_roundoff" "fa_roundoff"
                       "mesh_roundoff")
    oncvpsp_test(ae_file_${case} split ae_file_error_test.sh
        $<TARGET_FILE:oncvpspsr> $<TARGET_FILE:oncvpspsr_ae> $<TARGET_FILE:oncvpspsr_ps>
        "${data}/14_Si.dat" ${case} "${CMAKE_CURRENT_BINARY_DIR}/ae_file_${case}")
    if (case MATCHES "_roundoff$")
        oncvpsp_compares(ae_file_${case})
    endif ()
endforeach ()

# All-electron SCF convergence on a mesh starting at 1e-8 bohr: see scripts/ae_scf_test.sh
oncvpsp_test(ae_scf_nr_01_H mesh ae_scf_test.sh
    $<TARGET_FILE:oncvpspnr_ae> "${data}/01_H_ae.dat" "${CMAKE_CURRENT_BINARY_DIR}/ae_scf_nr_01_H")

# All-electron atom with a 4f shell: see scripts/ae_4f_test.sh
oncvpsp_test(ae_4f_scf_sr_70_Yb mesh ae_4f_test.sh
    $<TARGET_FILE:oncvpspsr_ae> "${data}/70_Yb_ae.dat" scf "${CMAKE_CURRENT_BINARY_DIR}/ae_4f_scf_sr_70_Yb")
oncvpsp_test(ae_4f_scf_r_70_Yb mesh ae_4f_test.sh
    $<TARGET_FILE:oncvpspr_ae> "${data}/70_Yb_ae.dat" scf "${CMAKE_CURRENT_BINARY_DIR}/ae_4f_scf_r_70_Yb")
oncvpsp_test(ae_4f_mesh_sr_69_Tm mesh ae_4f_test.sh
    $<TARGET_FILE:oncvpspsr_ae> "${data}/69_Tm_ae.dat" mesh "${CMAKE_CURRENT_BINARY_DIR}/ae_4f_mesh_sr_69_Tm")

# The all-electron atom against the NIST atomic reference data, to its stated accuracy, on the light
# atoms, a few heavy ones and the atoms of the cases expected to fail (EXPECTED_FAILURES in the
# script, which fails if one of them passes), neutral and first cations (see scripts/nist_test.py and
# tests/data/nist/README.md). The whole table is scripts/nist_test.py run by hand, e.g.
#   scripts/nist_test.py build/bin/oncvpspr_ae.x RLDA tests/data/nist build/nist_RLDA --jobs 8
foreach (nist_test IN ITEMS "LDA:nr" "ScRLDA:sr" "RLDA:r")
    string(REPLACE ":" ";" nist_test "${nist_test}")
    list(GET nist_test 0 set)
    list(GET nist_test 1 rel)
    set(parts "light:1-18" "heavy:29,54,64,79,92")
    if (set STREQUAL "ScRLDA")
        list(APPEND parts "known:62,65,68-71,87")
    endif ()
    foreach (part IN ITEMS ${parts})
        string(REPLACE ":" ";" part "${part}")
        list(GET part 0 part_name)
        list(GET part 1 zs)
        if (set STREQUAL "RLDA" AND NOT ONCVPSP_WITH_FR)
            continue()
        endif ()
        add_test(NAME nist_${set}_${part_name}
            COMMAND ${Python3_EXECUTABLE} "${scripts}/nist_test.py" $<TARGET_FILE:oncvpsp${rel}_ae>
                ${set} "${data}/nist" "${CMAKE_CURRENT_BINARY_DIR}/nist_${set}_${part_name}" --z ${zs})
        set_tests_properties(nist_${set}_${part_name} PROPERTIES LABELS nist)
    endforeach ()
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
        if (case MATCHES "^(defaults|explicit|rerun|upf|split|toml)$")
            oncvpsp_compares(log_mesh_${rel}_${test_name}_${case})
        endif ()
    endforeach ()
endforeach ()

# UPF file on the log mesh (upf_mesh, the optional third value of the output grid line): see
# scripts/upf_mesh_test.sh
foreach (mesh_test IN ITEMS "sr:14_Si_UPF" "r:52_Te")
    string(REPLACE ":" ";" mesh_test "${mesh_test}")
    list(GET mesh_test 0 rel)
    list(GET mesh_test 1 test_name)
    set(mesh_cases "defaults" "log" "stride")
    if (rel STREQUAL "sr")
        list(APPEND mesh_cases "negative" "read")
    endif ()
    if (ONCVPSP_WITH_TOML)
        list(APPEND mesh_cases "toml")
    endif ()
    foreach (case IN LISTS mesh_cases)
        oncvpsp_test(upf_mesh_${rel}_${test_name}_${case} mesh upf_mesh_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" ${case}
            "${CMAKE_CURRENT_BINARY_DIR}/upf_mesh_${rel}_${test_name}_${case}"
            "${data}/toml/${test_name}.toml")
        if (case MATCHES "^(defaults|log|stride|toml)$")
            oncvpsp_compares(upf_mesh_${rel}_${test_name}_${case})
        endif ()
    endforeach ()
endforeach ()

# TOML input: every input in tests/data/toml (the TOML form of its tests/data/*.dat) must give the
# output of the text input, and the reader must stop with a message naming the problem for each of the
# changes in scripts/toml_error_test.sh
if (ONCVPSP_WITH_TOML)
    set(toml_tests ${oncvpspsr_tests})
    list(TRANSFORM toml_tests PREPEND "sr:")
    set(toml_r_tests ${oncvpspr_tests})
    list(TRANSFORM toml_r_tests PREPEND "r:")
    list(APPEND toml_tests ${toml_r_tests} "nr:14_Si" "nr:29_Cu")
    foreach (toml_test IN LISTS toml_tests)
        string(REPLACE ":" ";" toml_test "${toml_test}")
        list(GET toml_test 0 rel)
        list(GET toml_test 1 test_name)
        oncvpsp_test(toml_${rel}_${test_name} toml toml_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${data}/toml/${test_name}.toml"
            "${CMAKE_CURRENT_BINARY_DIR}/toml_${rel}_${test_name}")
        oncvpsp_compares(toml_${rel}_${test_name})
    endforeach ()

    foreach (case IN ITEMS
            "parse_error" "missing_table" "not_a_table" "missing_key" "missing_string"
            "integer_type" "number_type" "string_type" "atsym_length" "state_counts"
            "lmax_range" "array_length" "array_type" "not_an_array" "l_order"
            "icmod3_rcfact" "test_config_length" "test_config_count" "test_config_type"
            "integer_for_real" "extra_keys" "icmod6_default_rcfact"
            "icmod6_default_slack")
        oncvpsp_test(toml_input_${case} toml toml_error_test.sh
            $<TARGET_FILE:oncvpspsr> "${data}" ${case} "${CMAKE_CURRENT_BINARY_DIR}/toml_input_${case}")
        if (case MATCHES "^(integer_for_real|extra_keys|icmod6_default_rcfact|icmod6_default_slack)$")
            oncvpsp_compares(toml_input_${case})
        endif ()
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
# rxpsh set), icmod=4 with its search grid (40_Zr_grid, 74_W), icmod=7 with slack (14_Si_hm), a
# Libxc iexc (83_Bi), rc on the grid (03_Li), and iexc=7 with c_light on the first line (52_Te_c,
# 14_Si_c).
set(te_rxpsh "s/^   -3\\.0  3\\.0  0\\.02$/   -3.0  3.0  0.02  3.5/")
set(te_c "s/^    Te  52\\.0   8   3   3   upf$/    Te  52.0   8   3   7   upf  137.0359895/")
set(si_c "s/^    Si  14     3   2   3   psp8$/    Si  14     3   2   7   psp8  137.0359895/")
foreach (echo_test IN ITEMS
        "sr:14_Si_lder:psp8" "sr:14_Si_lder:upf" "sr:14_Si_lder:psml"
        "sr:40_Zr_grid:psp8" "sr:40_Zr_grid:upf" "sr:40_Zr_grid:psml"
        "sr:83_Bi:psml" "sr:03_Li:psml" "nr:14_Si:psml" "sr:14_Si_hm:psp8"
        "r:74_W:psp8" "r:74_W:psml" "r:52_Te:psp8" "r:52_Te:upf" "r:52_Te:psml"
        "r:52_Te_c:psp8" "r:52_Te_c:upf" "r:52_Te_c:psml" "sr:14_Si_c:psp8" "sr:14_Si_c:psml")
    string(REPLACE ":" ";" echo_test "${echo_test}")
    list(GET echo_test 0 rel)
    list(GET echo_test 1 test_name)
    list(GET echo_test 2 format)
    set(edit "")
    set(input "${test_name}")
    if (test_name STREQUAL "52_Te")
        set(edit "${te_rxpsh}")
    elseif (test_name STREQUAL "52_Te_c")
        set(edit "${te_c}")
        set(input "52_Te")
    elseif (test_name STREQUAL "14_Si_c")
        set(edit "${si_c}")
        set(input "14_Si")
    endif ()
    oncvpsp_test(echo_${format}_${rel}_${test_name} echo echo_rerun_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${input}.dat" ${format}
        "${CMAKE_CURRENT_BINARY_DIR}/echo_${format}_${rel}_${test_name}" "${edit}")
    oncvpsp_compares(echo_${format}_${rel}_${test_name})
endforeach ()

# icmod=6 with fcfact<=0 (the cutoff taken from the pseudo-waves): the psp8 file must still carry the
# model core (see scripts/psp8_core_test.sh)
oncvpsp_test(psp8_core_sr_14_Si psp8 psp8_core_test.sh
    $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat" "s/^   1    0\\.25$/   6    0.0   1.0d-3/"
    "${CMAKE_CURRENT_BINARY_DIR}/psp8_core_sr_14_Si")
oncvpsp_test(psp8_core_r_74_W psp8 psp8_core_test.sh
    $<TARGET_FILE:oncvpspr> "${data}/74_W.dat" "s/^    4  0\\.0$/    6  0.0  1.0d-3/"
    "${CMAKE_CURRENT_BINARY_DIR}/psp8_core_r_74_W")
oncvpsp_compares(psp8_core_sr_14_Si)
oncvpsp_compares(psp8_core_r_74_W)

# Text input from a pipe (as from a Python subprocess), which cannot be backspaced: the same output
# as from a file, and an input ending early reported with the line the reader stopped at
foreach (pipe_test IN ITEMS "sr:14_Si" "r:52_Te" "nr:08_O")
    string(REPLACE ":" ";" pipe_test "${pipe_test}")
    list(GET pipe_test 0 rel)
    list(GET pipe_test 1 test_name)
    oncvpsp_test(pipe_${rel}_${test_name} pipe pipe_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${CMAKE_CURRENT_BINARY_DIR}/pipe_${rel}_${test_name}")
    oncvpsp_compares(pipe_${rel}_${test_name})
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
foreach (case IN ITEMS "inside:2.0" "at_rc:2.70")
    string(REPLACE ":" ";" case "${case}")
    list(GET case 0 where)
    list(GET case 1 radius)
    oncvpsp_test(error_rxpsh_${where} errors error_test.sh
        $<TARGET_FILE:oncvpspsr> "check_data: ERROR rxpsh not beyond every core radius"
        "${CMAKE_CURRENT_BINARY_DIR}/error_rxpsh_${where}"
        "${data}/14_Si_lder.dat" "s/^   -0\\.3  2\\.0  0\\.02  3\\.5$/   -0.3  2.0  0.02  ${radius}/")
endforeach ()
oncvpsp_test(error_icmod7_slack_range errors error_test.sh
    $<TARGET_FILE:oncvpspsr> "test_data: must have 0.0<slack<=1.0 for icmod= 6 or 7"
    "${CMAKE_CURRENT_BINARY_DIR}/error_icmod7_slack_range"
    "${data}/14_Si_hm.dat" "s/^   7   0   1\\.0e-5   0\\.7$/   7   0   1.0e-5   1.5/")
oncvpsp_test(error_iexc5_r2scan errors error_test.sh
    $<TARGET_FILE:oncvpspsr> "test_data: iexc=5 is r2SCAN, METAPSP only"
    "${CMAKE_CURRENT_BINARY_DIR}/error_iexc5_r2scan"
    "${data}/14_Si.dat" "s/^    Si  14     3   2   3   psp8$/    Si  14     3   2   5   psp8/")
oncvpsp_test(error_icmod4_bad_grid errors error_test.sh
    $<TARGET_FILE:oncvpspsr> "check_data: ERROR icmod=4 needs 0<fcfact_min<=fcfact_max and fcfact_step>0"
    "${CMAKE_CURRENT_BINARY_DIR}/error_icmod4_bad_grid"
    "${data}/40_Zr.dat" "s/^    4   0\\.0$/    4   0.0   0.0   2.0   5.0   0.0   1.1   1.6   0.05/")
