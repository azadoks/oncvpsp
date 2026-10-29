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

set(scripts "${PROJECT_SOURCE_DIR}/scripts")
set(data "${CMAKE_CURRENT_SOURCE_DIR}/data")

# The inputs with Libxc functionals (iexc < 0), whose tests are disabled without Libxc
set(libxc_inputs "07_N_opt" "07_N_flat" "38_Sr_lxc" "60_Nd_GHOST" "79_Au_lxc" "83_Bi")

# oncvpsp_add_test(LABEL FAMILY REL CASE COMMAND ARGS...): every test is registered here. Its name
# is FAMILY_TREATMENT_CASE, or FAMILY_CASE for REL "-" (a test of no one treatment), and its
# labels are LABEL, the kind of test, and TREATMENT: nr, sr or fr (REL r, as in the names of the
# executables, oncvpspr*). CASE is the input, then what varies (14_Si_upf). COMMAND is a script of
# scripts/, run by bash (.sh) or the test Python (.py), or a program; @DIR@ in ARGS is the test's
# own directory. Not added for REL r without ONCVPSP_WITH_FR; disabled when CASE is an input of
# libxc_inputs (with or without a suffix) without ONCVPSP_WITH_Libxc. Returns the name in
# added_test ("" when not added).
function(oncvpsp_add_test label family rel case command)
    set(added_test "" PARENT_SCOPE)
    if (rel STREQUAL "-")
        set(name ${family}_${case})
        set(labels ${label})
    else ()
        if (rel STREQUAL "r")
            if (NOT ONCVPSP_WITH_FR)
                return()
            endif ()
            set(treatment fr)
        else ()
            set(treatment ${rel})
        endif ()
        set(name ${family}_${treatment}_${case})
        set(labels ${label} ${treatment})
    endif ()
    set(args "")
    foreach (arg IN LISTS ARGN)
        string(REPLACE "@DIR@" "${CMAKE_CURRENT_BINARY_DIR}/${name}" arg "${arg}")
        list(APPEND args "${arg}")
    endforeach ()
    if (command MATCHES "\\.sh$")
        add_test(NAME ${name} COMMAND ${ONCVPSP_BASH} "${scripts}/${command}" ${args})
    elseif (command MATCHES "\\.py$")
        add_test(NAME ${name} COMMAND ${test_python} "${scripts}/${command}" ${args})
    else ()
        add_test(NAME ${name} COMMAND ${command} ${args})
    endif ()
    set_tests_properties(${name} PROPERTIES LABELS "${labels}" ENVIRONMENT "ONCVPSP_PYTHON=${test_python}")
    if (NOT ONCVPSP_WITH_Libxc)
        foreach (input IN LISTS libxc_inputs)
            if (case MATCHES "^${input}(_|$)")
                set_tests_properties(${name} PROPERTIES DISABLED TRUE)
            endif ()
        endforeach ()
    endif ()
    set(added_test ${name} PARENT_SCOPE)
endfunction()

# The Fortran unit tests: test_NAME from unit/test_NAME.f90 and any more SOURCES, linked to the
# library and to the assertions they share (unit/unit_check_m.f90)
add_library(unit_check STATIC unit/unit_check_m.f90)
target_link_libraries(unit_check PUBLIC oncvpsp::oncvpsp PRIVATE oncvpsp::compile_options)
set_target_properties(unit_check PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/unit_modules")
target_include_directories(unit_check PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/unit_modules")
function(oncvpsp_unit_executable name)
    add_executable(test_${name} unit/test_${name}.f90 ${ARGN})
    target_link_libraries(test_${name} PRIVATE unit_check oncvpsp::oncvpsp oncvpsp::compile_options)
    set_target_properties(test_${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()

# Unit tests of the input readers
oncvpsp_unit_executable(input_text)
oncvpsp_add_test(unit unit - input_text test_input_text)

# Unit tests of the rxpsh warning in run_phsft (nr and fr), called directly
oncvpsp_unit_executable(run_phsft)
oncvpsp_add_test(unit unit - run_phsft test_run_phsft)
set_tests_properties(${added_test} PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Unit tests of oncvpsp_run called as from the library: its test results against the DIAGNOSTICS
# block it writes, and a second run in the same process (test_run.f90); nr, sr, fr, and the weighted
# cutoff of the model core's budget
oncvpsp_unit_executable(run)
foreach (run_test IN ITEMS "0:nr:08_O" "1:sr:14_Si" "2:r:14_Si" "1:sr:14_Si_bcw")
    string(REPLACE ":" ";" run_test "${run_test}")
    list(GET run_test 0 irel)
    list(GET run_test 1 rel)
    list(GET run_test 2 name)
    oncvpsp_add_test(unit unit_run ${rel} ${name} test_run ${irel} "${data}/${name}.dat")
endforeach ()

# Unit test of the stages of oncvpsp_run on their own (test_stages.f90): the all-electron atom once,
# every channel twice from it, bit for bit the same and the same as the whole run's; nr, sr, fr,
# and a local channel (lloc 2, 07_N_dl)
oncvpsp_unit_executable(stages)
foreach (stage_test IN ITEMS "0:nr:08_O" "1:sr:14_Si" "2:r:14_Si" "0:nr:07_N_dl")
    string(REPLACE ":" ";" stage_test "${stage_test}")
    list(GET stage_test 0 irel)
    list(GET stage_test 1 rel)
    list(GET stage_test 2 name)
    oncvpsp_add_test(unit unit_stages ${rel} ${name} test_stages ${irel} "${data}/${name}.dat")
endforeach ()

# Unit test of run_coupled called again on one run's stages (test_coupled.f90): with the input's
# model core, with none, and with the input's again, each the same as a whole run, bit for bit;
# the budget (icmod 6) and the polynomial (icmod 1), sr and fr
oncvpsp_unit_executable(coupled)
foreach (coupled_test IN ITEMS "1:sr:14_Si_bc" "1:sr:14_Si" "2:r:14_Si")
    string(REPLACE ":" ";" coupled_test "${coupled_test}")
    list(GET coupled_test 0 irel)
    list(GET coupled_test 1 rel)
    list(GET coupled_test 2 name)
    oncvpsp_add_test(unit unit_coupled ${rel} ${name} test_coupled ${irel} "${data}/${name}.dat")
endforeach ()

# Unit test of the failure record (test_failure.f90): an ifail split into its kind, message and
# status, and two runs that fail, at the input and in a channel, with the message and the stage
oncvpsp_unit_executable(failure)
oncvpsp_add_test(unit unit_failure nr 08_O test_failure "${data}/08_O.dat")

# Test of the C interface (c/test_c_api.c, liboncvpsp): two handles at once running one TOML input,
# with the same sizes and no failure, and a broken input failing at the input stage
if (ONCVPSP_WITH_C_API AND ONCVPSP_WITH_TOML)
    add_executable(test_c_api c/test_c_api.c)
    target_link_libraries(test_c_api PRIVATE oncvpsp::c)
    set_target_properties(test_c_api PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    foreach (c_test IN ITEMS "0:nr:08_O" "1:sr:14_Si" "2:r:14_Si")
        string(REPLACE ":" ";" c_test "${c_test}")
        list(GET c_test 0 irel)
        list(GET c_test 1 rel)
        list(GET c_test 2 name)
        oncvpsp_add_test(unit c_api ${rel} ${name} test_c_api ${irel} "${data}/toml/${name}.toml"
            "${CMAKE_CURRENT_BINARY_DIR}/c_api_${rel}_${name}.out")
    endforeach ()
endif ()

# Unit test of oncvpsp_run in several threads at once (test_run_threads.f90): each case alone, then
# all of them twice over concurrently, which must give the same results and text bit for bit, so the
# library keeps no state between or across runs; the cases differ in irel, c_light and the
# functional (with Libxc, its GGA in icmod 4 and in icmod 6, whose budget calls LAPACK). A threaded
# BLAS runs on all its threads outside the parallel region and on one inside, whose sums differ in
# their last bits (07_N_opt with MKL), so it is held to one thread. Only the test is compiled with
# OpenMP, and it is left out without it
find_package(OpenMP COMPONENTS Fortran)
if (OpenMP_Fortran_FOUND)
    oncvpsp_unit_executable(run_threads)
    target_link_libraries(test_run_threads PRIVATE OpenMP::OpenMP_Fortran)
    set(thread_cases "0:${data}/08_O.dat:0" "1:${data}/14_Si.dat:0" "1:${data}/14_Si.dat:137.0359991"
                     "2:${data}/14_Si.dat:0" "2:${data}/14_Si.dat:100")
    if (ONCVPSP_WITH_Libxc)
        list(APPEND thread_cases "1:${data}/38_Sr_lxc.dat:0" "1:${data}/07_N_opt.dat:0")
    endif ()
    oncvpsp_add_test(unit unit - run_threads test_run_threads ${thread_cases})
    set_property(TEST unit_run_threads APPEND PROPERTY ENVIRONMENT
                 "OPENBLAS_NUM_THREADS=1" "MKL_NUM_THREADS=1" "BLIS_NUM_THREADS=1")
endif ()

# Unit tests of the operations at any radius on the log mesh (logmesh_m), the closed-form Bessel
# overlaps (sbf_overlap) and the derivatives at rc from the radial equation (ode_rc_der), phase 1 of
# doc/smoothness_plan.md
foreach (unit_test logmesh ode_rc_der solver_path piecewise matching ugrid quadbasis pieces vploc ae_delta outwf outwf_index sbf free_tail qroots sr_so_r phsft adams fxc)
    oncvpsp_unit_executable(${unit_test})
    oncvpsp_add_test(unit unit - ${unit_test} test_${unit_test})
endforeach ()

# Unit tests of the pieces of the model core charges (modcore_m), phase 5 of doc/smoothness_plan.md
# the jump of the basis wave vectors where ulgd passes l/rc (test_qroots.f90):
# known, so the probe passes while it fails; drop WILL_FAIL once fixed
oncvpsp_add_test(unit unit - qroots_lcross test_qroots lcross)
set_tests_properties(${added_test} PROPERTIES WILL_FAIL TRUE)
oncvpsp_unit_executable(modcore unit/gg1cc.f90 unit/gp1cc.f90 unit/gpp1cc.f90)
# imported, so not held to the Debug build's -std=f2008 and -Wconversion
set_source_files_properties(unit/gg1cc.f90 unit/gp1cc.f90 unit/gpp1cc.f90 PROPERTIES
    COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-std=gnu;-Wno-conversion>")
oncvpsp_add_test(unit unit - modcore test_modcore)

# oncvpsp_compares(NAME): the test NAME (if it was added) compares outputs with
# scripts/compare_output.py, or parses them with pyoncvpsp: disabled without its Python
function(oncvpsp_compares name)
    if (TEST "${name}" AND NOT ONCVPSP_COMPARE_PYTHON)
        set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    endif ()
endfunction()

# oncvpsp_uses_run(NAME REL INPUT [AE]): the test NAME (if it was added) reads the output of
# oncvpsp<REL>.x for tests/data/INPUT.dat from the fixture run_TREATMENT_INPUT ($ONCVPSP_RUN), or
# with AE the all-electron file of oncvpsp<REL>_ae.x from run_ae_TREATMENT_INPUT
# ($ONCVPSP_AE_FILE): each input is run once for all the tests that read it (scripts/run_test.sh)
function(oncvpsp_uses_run name rel input)
    if (NOT TEST "${name}")
        return()
    endif ()
    if (rel STREQUAL "r")
        set(treatment fr)
    else ()
        set(treatment ${rel})
    endif ()
    if (ARGN STREQUAL "AE")
        set(family run_ae)
        set(fixture run_ae_${treatment}_${input})
        set(file "${CMAKE_CURRENT_BINARY_DIR}/${fixture}/atom.ae")
        set(variable ONCVPSP_AE_FILE)
        set(exe $<TARGET_FILE:oncvpsp${rel}_ae>)
        set(ae_file "${file}")
        set(out "${CMAKE_CURRENT_BINARY_DIR}/${fixture}/ae.out")
    else ()
        set(family run)
        set(fixture run_${treatment}_${input})
        set(file "${CMAKE_CURRENT_BINARY_DIR}/${fixture}/${input}.out")
        set(variable ONCVPSP_RUN)
        set(exe $<TARGET_FILE:oncvpsp${rel}>)
        set(ae_file "")
        set(out "${file}")
    endif ()
    if (NOT TEST ${fixture})
        oncvpsp_add_test(run ${family} ${rel} ${input} run_test.sh ${exe} "${data}/${input}.dat" "${out}"
            ${ae_file})
        set_tests_properties(${fixture} PROPERTIES FIXTURES_SETUP ${fixture})
    endif ()
    set_property(TEST ${name} APPEND PROPERTY FIXTURES_REQUIRED ${fixture})
    set_property(TEST ${name} APPEND PROPERTY ENVIRONMENT "${variable}=${file}")
endfunction()

# Smooth dependence on the inputs (doc/smoothness_plan.md, T2): smooth_scan.py steps one input of a
# TOML input through 21 to 41 values (by the fit's degree) and fits each feature of the UPF file on the log mesh (nonlocal
# kernels, V_loc, model core) with a polynomial; the residuals must be within the plan's goals and
# no step between neighbouring values a jump. Label scan (ctest -LE scan leaves them out). The
# harness checks itself (unit).
# The coding conventions that fortitude and gfortran do not check (scripts/check_conventions.py,
# doc/coding_standards.md), on all of src/: no finding beyond scripts/conventions_baseline.txt
oncvpsp_add_test(lint lint - conventions check_conventions.py --all)

oncvpsp_add_test(unit unit - smooth_scan smooth_scan.py --self-test)
oncvpsp_compares(smooth_scan_self_test)

# The weights of logmesh_m's stencils (src/logmesh_weights.inc) are those scripts/logmesh_weights.py writes
oncvpsp_add_test(unit unit - logmesh_weights logmesh_weights.py --check
    "${PROJECT_SOURCE_DIR}/src/logmesh_weights.inc")

# kernel_diff.py and compare_output.py give the same kernel differences: 60_Nd_GHOST_r's reference
# output against itself with its psp8 projectors changed (scripts/kernel_agreement_test.py)
oncvpsp_add_test(unit unit - kernel_agreement kernel_agreement_test.py
    "${CMAKE_CURRENT_SOURCE_DIR}/refs/60_Nd_GHOST_r.out" "@DIR@")
oncvpsp_compares("${added_test}")

# The live docs name only files and modules that exist (scripts/docrefs.py)
oncvpsp_add_test(lint lint - docrefs docrefs.py "${PROJECT_SOURCE_DIR}")

# goals G1 (an input that was rounded or searched) and G2 (one that never was: the noise floor)
set(smooth_g1 "kernel=1e-10,vloc=1e-10,core=1e-12,tests=1e-10")
set(smooth_g2 "kernel=1e-12,vloc=1e-12,core=1e-12,tests=1e-10")
# tests: the test configurations' eigenvalues and energy differences from the DIAGNOSTICS block
# (phase 6), at G1 in every row whose input has test configurations (not 07_N_opt): <= 3.9e-12
# measured but in 08_O rc (3.5e-11: PZ81's jump at r_s = 1 crossing nodes near the origin in the
# pseudo 2s of 2s1 2p4, F-1) and 20_Ca rcloc (9.0e-11: a kink where rcloc crosses 0.95 rc(s), two
# breakpoints of outwf_m's outward integration, R-16; doc/smoothness_plan.md). Not
# in the rows that move r(r_s = 1) through the pseudo-atom's density (14_Si and 34_Se fcfact, 34_Se
# rcfact, 52_Te rc[*]): PZ81's jump there (F-1, kept, D-6) steps the energy differences, residuals
# 2.8e-7 to 7.0e-7 (14_Si fcfact with PW92 or PBE: 7e-15, 1.1e-14)
# smooth_scan_test(CASE REL INPUT KEY GROUPS TOL WINDOW...): the test scan_TREATMENT_CASE, the scan of KEY of tests/data/toml/INPUT.toml
# by oncvpsp<REL>.x over WINDOW (smooth_scan.py's --center/--cells/--span), features GROUPS. Each row
# runs ONCVPSP_SMOOTH_JOBS oncvpsp processes at once and reserves as many ctest processors: 4 for
# the whole label (ctest -j N runs N/4 rows), the number of cores to run one row in the dev loop
set(ONCVPSP_SMOOTH_JOBS 4 CACHE STRING "Concurrent oncvpsp runs of each smooth_scan test row")
function(smooth_scan_test case rel input key groups tol)
    if (NOT ONCVPSP_WITH_TOML)
        return()
    endif ()
    oncvpsp_add_test(scan scan ${rel} ${case} smooth_scan.py $<TARGET_FILE:oncvpsp${rel}>
        "${data}/toml/${input}.toml" "${key}" ${ARGN} --groups ${groups} --tol ${tol}
        --max-jump 3 --jobs ${ONCVPSP_SMOOTH_JOBS} --workdir "@DIR@")
    if (added_test)
        set_tests_properties(${added_test} PROPERTIES PROCESSORS ${ONCVPSP_SMOOTH_JOBS})
        oncvpsp_compares("${added_test}")
    endif ()
endfunction()
# core radii, local-potential radius (lpopt 5 and 1), qcut (phase 3). KEY[*] moves every channel's
# value by the same offset (each channel's kernel its own feature), one row for what were one per
# channel; with the psp8 and PSML files (phase 4, O-1: the writers interpolated the functions on
# the log mesh across rc, so the points next to it moved with the nodes past it) and the UPF file
# on the linear mesh, a row of its own (one run writes one UPF file). 14_Si's rc[1] and rcloc bend
# the kernel enough over the window that degree 6 leaves 1e-9 and 1e-5 (degree 10: 3e-11 for rc[1]).
# 52_Te's rc[*] moves r(r_s = 1), where PZ81's V_xc jumps (F-1, kept, D-6), across the nodes the
# linear meshes' interpolants reach (1e-8 in V_loc there): the linear meshes with rc[1] alone
smooth_scan_test(14_Si_rc sr 14_Si "pseudopotentials.rc[*]" kernel,vloc,core,tests ${smooth_g1} --center 0.1 --degree 14
    --output log,psp8,psml)
smooth_scan_test(14_Si_rc_upf sr 14_Si "pseudopotentials.rc[*]" kernel,vloc,tests ${smooth_g1} --center 0.1 --degree 14
    --output upf)
smooth_scan_test(14_Si_rcloc sr 14_Si "local_potential.rcloc" kernel,vloc,core,tests ${smooth_g1} --center 2.6 --degree 14)
smooth_scan_test(20_Ca_rcloc sr 20_Ca "local_potential.rcloc" kernel,vloc,tests ${smooth_g1} --center 1.8)
smooth_scan_test(14_Si_qcut sr 14_Si "pseudopotentials.qcut[*]" kernel,vloc,core,tests ${smooth_g1} --span 0.04)
smooth_scan_test(08_O_rc nr 08_O "pseudopotentials.rc[*]" kernel,vloc,tests ${smooth_g1} --center 0.05 --output log,psp8)
smooth_scan_test(52_Te_rc r 52_Te "pseudopotentials.rc[*]" kernel,vloc ${smooth_g1} --center 0.05)
smooth_scan_test(52_Te_rc1_linear r 52_Te "pseudopotentials.rc[1]" kernel,vloc,tests ${smooth_g1} --center 3.15
    --output psp8,upf)
# inputs never rounded: the noise floor (phase 2), fitted with degree 10 (the kernel bends enough in
# dvloc0 that degree 6 leaves 2e-9 over the window; 10 still sees a jump of 1e-7 at 2e-8)
smooth_scan_test(14_Si_dvloc0 sr 14_Si "local_potential.dvloc0" kernel,vloc,core,tests ${smooth_g2} --span 0.02 --degree 10)
smooth_scan_test(14_Si_debl sr 14_Si "vkb_projectors.debl[*]" kernel,vloc,core,tests ${smooth_g2} --span 0.02 --degree 10)
smooth_scan_test(14_Si_ep2 sr 14_Si "pseudopotentials.ep[2]" kernel,vloc,core,tests ${smooth_g2} --span 0.02 --degree 10)
smooth_scan_test(52_Te_dvloc0 r 52_Te "local_potential.dvloc0" kernel,vloc,tests ${smooth_g2} --span 0.02 --degree 10)
# model core, icmod 1, 2, 3, 4, 6, 7 (phase 5), and the model core of the psp8 file with its four
# derivatives, written from its pieces on each side of the joint (r*, r_s; modcore_m core_pieces).
# V_loc of the psp8 file is left out: the core moves r(r_s = 1), as rc does above (1e-8 at 1 bohr,
# 14_Si; 9e-8, 34_Se)
smooth_scan_test(14_Si_fcfact sr 14_Si "model_core_charge.fcfact" kernel,vloc.log,core ${smooth_g1} --span 0.02
    --output log,psp8)
smooth_scan_test(79_Au_lxc_fcfact sr 79_Au_lxc "model_core_charge.fcfact" kernel,vloc,core,tests ${smooth_g1} --span 0.02)
smooth_scan_test(34_Se_fcfact sr 34_Se "model_core_charge.fcfact" kernel,vloc,core ${smooth_g1} --span 0.02)
smooth_scan_test(34_Se_rcfact sr 34_Se "model_core_charge.rcfact" kernel,vloc.log,core ${smooth_g1} --span 0.02
    --output log,psp8)
# icmod 4: the optimum does not move with the grid, so the scan sees only where Newton stops on
# the noise of der2exc (~1e-17 on f = 1.2e-4) over the soft curvature of the valley (5e-4):
# ~8e-13 in ln(amplitude), ~1e-12 in the core (M-6)
smooth_scan_test(40_Zr_grid_fcfact_min sr 40_Zr_grid "model_core_charge.fcfact_min" kernel,vloc,core,tests
    "kernel=1e-10,vloc=1e-10,core=3e-12,tests=1e-10" --span 0.02)
# icmod 6/7 (M-8..M-12): r_s is a root of the error of the model's density cutoffs, whose
# integral over [0, r_s] is taken on a log mesh ending at r_s (on the fixed mesh its error rippled
# with r_s, ~1e-9 in r_s, the residuals of 14_Si_bc's scans until M-12 was fixed). What is left is
# the fit: degree 18 for the curvature over the wide fcfact span (14 leaves 1e-10), and 14_Si_hm's
# xc round-off (~2e-11 in r_s). The psp8 file's 3rd and 4th columns of the core jump at r_s and are
# steep just inside it: at 0.36 bohr (r_s 0.367) round-off in r_s leaves 1.2e-10 in the 3rd (the
# same at span 0.02, any degree), and the span of 0.2 that moves r_s across mesh nodes needs degree
# 24 (18 leaves 1.8e-9 of truncation; the interpolant across the joint: 3.0e-2, kinked). The 3rd
# column's round-off there is 1.1e-10 to 4.3e-10 with 37 to 45 samples and either q grid of
# run_optimize (dq = 2 pi/(150 rc) or 2 pi/(300 rc)): its tolerance above it, 1e-9
smooth_scan_test(14_Si_bc_ecut sr 14_Si_bc "model_core_charge.fcfact" kernel,vloc,core,tests
    "kernel=1e-10,vloc=1e-11,core=2e-11,core.psp8=1e-9,tests=1e-10" --span 0.2 --degree 18,psp8=24 --output log,psp8)
# the onset of the well (well_onset) where 14_Si's well states keep their n (l=1 takes another
# node between 1.06 and 1.1)
smooth_scan_test(14_Si_well_onset sr 14_Si "pseudopotentials.well_onset" kernel,vloc,core,tests ${smooth_g1} --center 1.02 --span 0.01)
# icmod 6/7: the xc round-off (14_Si_hm), as above; 14_Si_bc at degree 10 (6 leaves 3e-9). 14_Si_hm
# at its budget of before, over every wave ("max_all", 33.9 Ha): over the ground-state waves (11.5 Ha)
# r_s is not the root of the budget there and steps by up to 1e-3 (smoothness_plan, bugs found)
smooth_scan_test(14_Si_hm_slack sr 14_Si_hm "model_core_charge.slack" kernel,vloc,core,tests
    "kernel=1e-10,vloc=1e-10,core=1e-10,tests=1e-10" --span 0.02 --degree 14
    --set "model_core_charge.budget_cutoff=\"max_all\"")
smooth_scan_test(14_Si_bc_slack sr 14_Si_bc "model_core_charge.slack" kernel,vloc,core,tests
    "kernel=1e-10,vloc=3e-12,core=1e-11,tests=1e-10" --center 0.9 --span 0.02 --degree 10)
# the optimizer's setup (07_N_opt: PBE, icmod 6 at fcfact 65, lpopt 5): the inputs it varies and
# the outputs it reads. ld and ecut from the DIAGNOSTICS block (O-6, 17 digits; up to 1.2e-10 and
# 3.7e-10 measured; ld 3.0e-10 in rcloc since the exact Adams-Moulton step, a smooth feature, jump
# ratio 1.03, l=2 at 21 Ha, that falls with the order: 2.7e-10 at 6, 1.9e-10 at 8); the core at the noise of icmod 6's search (~1e-11); the ghost tests' margins
# from the block too (X-3, up to 9.5e-13 measured); ecut has the local potential's cutoffs too
# (vloc_cutoff, up to 2.4e-12, 1.6e-11 in rc1_wcut). At the default rc[0] = rcloc, and
# rc < rcloc is not allowed, so the rc are scanned 0.01 above their values (rc[0] = 1.21) and rcloc
# below
set(smooth_opt "kernel=1e-10,vloc=1e-10,core=5e-11,ld=4e-10,ecut=1e-9,ghost=3e-12")
set(smooth_opt_groups kernel,vloc,core,ld,ecut,ghost)
smooth_scan_test(07_N_opt_rc sr 07_N_opt "pseudopotentials.rc[*]" ${smooth_opt_groups} ${smooth_opt} --center 0.01)
smooth_scan_test(07_N_opt_qcut sr 07_N_opt "pseudopotentials.qcut[*]" ${smooth_opt_groups} ${smooth_opt} --span 0.04)
smooth_scan_test(07_N_opt_debl sr 07_N_opt "vkb_projectors.debl[*]" ${smooth_opt_groups} ${smooth_opt} --span 0.02)
smooth_scan_test(07_N_opt_rcloc sr 07_N_opt "local_potential.rcloc" ${smooth_opt_groups} ${smooth_opt} --center 1.19)
# rxpsh, straddling the node 6e-6 above it (X-2): the phases are taken at rxpsh itself, no longer
# at the first node beyond it (a staircase). At 11-12 Ha the l=2 angle turns by nearly pi within
# ~0.01 bohr (k r ~ 19, so (phi, r phi') traces a thin ellipse), smooth but steep: a narrow window
# and degree 10 (1.8e-11)
smooth_scan_test(07_N_opt_rxpsh sr 07_N_opt "log_derivative_analysis.rxpsh" ld "ld=3e-10" --span 0.004 --degree 10)
# the free tail of the well states (well_tail = "free", phase 8): the optimizer's rc[1], whose well
# state's convergence table bends more (its 1e-5 crossing near 75 Ha, normalized inside rc), so
# degree 14 (5e-12; degree 6 leaves 5e-7 of truncation); 14_Si's rc[1] (3.8e-12); and ep[2] of
# 14_Si's l = 2, whose first projector is a well state, with ep_as_given, from 0 up: the free tail
# is for ep >= 0 only (D-7, 2026-10-23), below 0 an input error (this row crossed 0 before, with no
# branch there, D-7b). The ecut group is left out of the ep[2] row: its 1e-5 crossing slid 0.7 Ha
# between -0.014 and -0.012 Ha, a shoulder of the curve passing the level (G5)
set(smooth_free "pseudopotentials.well_tail=\"free\"")
smooth_scan_test(07_N_opt_rc1_free sr 07_N_opt "pseudopotentials.rc[1]" ${smooth_opt_groups} ${smooth_opt}
    --degree 14 --set ${smooth_free})
# the weighted cutoff of the pseudo-waves for the model core's budget (budget_cutoff = "weighted",
# fcfact 0): the occupation-weighted mean of the convergence curves, smooth where two channels'
# crossings pass each other (the largest is not); ecut:weighted to 1.2e-11 measured; the core
# to 2e-10 (3.0e-11 at -O2, 1.06e-10 under -ffast-math since 43f2c6d, with no jump)
set(smooth_wcut "model_core_charge.budget_cutoff=\"weighted\"")
smooth_scan_test(07_N_opt_rc1_wcut sr 07_N_opt "pseudopotentials.rc[1]" kernel,vloc,core,ecut
    "kernel=1e-10,vloc=1e-10,core=2e-10,ecut=1e-9"
    --set model_core_charge.fcfact=0.0 --set ${smooth_wcut})
smooth_scan_test(14_Si_rc1_free sr 14_Si "pseudopotentials.rc[1]" kernel,vloc,core,tests ${smooth_g1} --center 2.8
    --degree 14 --set ${smooth_free})
smooth_scan_test(14_Si_ep2_free sr 14_Si "pseudopotentials.ep[2]" kernel,vloc,core,tests ${smooth_g1} --center 0.02
    --span 0.04 --degree 10 --set ${smooth_free} --set pseudopotentials.ep_as_given=true)
# the confined well (well_tail = "confined"): 14_Si's l = 0 well state across debl[0] = 0.8425, where
# Hamann's well moves it from n = 4 to n = 5 (its kernel residual 6e-3, a step 825 times the median);
# the confined well keeps n = 4 (1.2e-13)
set(smooth_confined "pseudopotentials.well_tail=\"confined\"")
smooth_scan_test(14_Si_debl0_confined sr 14_Si "vkb_projectors.debl[0]" kernel,vloc,core,tests ${smooth_g1}
    --center 0.8425 --span 0.01 --set ${smooth_confined})

# The projectors of the psp8 and UPF files on the linear mesh are 0 beyond rc, and the fully-
# relativistic psp8 file has every term (phase 4, O-1, O-2; see scripts/projector_output_test.py)
foreach (projector_test IN ITEMS "sr:14_Si" "r:14_Si" "nr:08_O" "r:52_Te")
    string(REPLACE ":" ";" projector_test "${projector_test}")
    list(GET projector_test 0 rel)
    list(GET projector_test 1 test_name)
    if (NOT ONCVPSP_WITH_TOML)
        continue()
    endif ()
    oncvpsp_add_test(output projector_output ${rel} ${test_name} projector_output_test.py
        $<TARGET_FILE:oncvpsp${rel}> "${data}/toml/${test_name}.toml" "@DIR@")
    oncvpsp_compares("${added_test}")
endforeach ()

# Output against the references in tests/refs. The lists here and below are the fewest inputs, and
# on them the cheapest tests, that keep the coverage of each kind of test (the lines and branches of
# src/ that the reference tests run, those the TOML tests run, ...) as it was with every input of
# tests/data (2026-10-22, gcov; doc/smoothness_plan.md, "Optimize the test suite")
set(oncvpspsr_tests
    "03_Li"         # no core states
    "14_Si"         # icmod=1
    "14_Si_bc"      # icmod=6
    "14_Si_bca"     # icmod=6 with fcfact=0 (the cutoff from the pseudo-waves)
    "20_Ca"         # lpopt=1
    "32_Ge"         # psfile both
    "34_Se"         # icmod=3
    "40_Zr"         # icmod=4
    "60_Nd_GHOST"   # ghosts
    "79_Au_lxc"     # icmod=2
    "73_Ta"         # the example of doc/users_guide.md (replot.sh, doc/73_Ta_plots.pdf)
    "40_Zr_grid"    # icmod=4 with an explicit search grid (modcore3)
    "14_Si_hm"      # icmod=7 with slack (modcore_budget, modcore_probes)
    "14_Si_bcw"     # icmod=6 with the weighted cutoff of the waves (modcore_budget_table)
    "07_N_opt"      # the optimizer's setup (oncvpsp-opt): PBE (Libxc), icmod 6, rxpsh, its own
                    # log mesh, where the AE SCF stops on its residual's floor (ae_converged)
    "07_N_flat"     # the same at qcut[0] = 22, where E_resid(qcut) is flat to round-off and the
                    # wave had a node: wave_regularization 1e-8 (optimize)
)
foreach (input IN LISTS oncvpspsr_tests)
    oncvpsp_add_test(ref ref sr ${input} ref_test.sh
        $<TARGET_FILE:oncvpspsr> "${data}/${input}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${input}.out" "@DIR@")
    oncvpsp_compares("${added_test}")
    oncvpsp_uses_run("${added_test}" sr ${input})
endforeach ()

set(oncvpspr_tests
    "52_Te" "80_Hg" # the fr inputs the docs measure; the fr mode tests run cheaper ones
    "07_N_dl"       # lloc<4
    "08_O_wigner"   # iexc=1
    "14_Si_bc"      # icmod=6
    "14_Si_hm"      # icmod=7 with slack
    "20_Ca"         # lpopt=1
    "38_Sr_lxc"     # a single projector for one l
    "60_Nd_GHOST"   # ghosts
    "79_Au_lxc"     # icmod=2
)
foreach (input IN LISTS oncvpspr_tests)
    oncvpsp_add_test(ref ref r ${input} ref_test.sh
        $<TARGET_FILE:oncvpspr> "${data}/${input}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${input}_r.out" "@DIR@")
    oncvpsp_compares("${added_test}")
    oncvpsp_uses_run("${added_test}" r ${input})
endforeach ()

set(oncvpspnr_tests "07_N_dl")
foreach (input IN LISTS oncvpspnr_tests)
    oncvpsp_add_test(ref ref nr ${input} ref_test.sh
        $<TARGET_FILE:oncvpspnr> "${data}/${input}.dat" "${CMAKE_CURRENT_SOURCE_DIR}/refs/${input}_nr.out" "@DIR@")
    oncvpsp_compares("${added_test}")
    oncvpsp_uses_run("${added_test}" nr ${input})
endforeach ()

# All-electron / pseudopotential split: oncvpsp_ae.x writes the all-electron atom to a file and
# oncvpsp_ps.x reads it (see scripts/split_test.sh)
foreach (split_test IN ITEMS "sr:14_Si" "sr:57_La" "nr:14_Si" "nr:29_Cu" "r:03_Li")
    string(REPLACE ":" ";" split_test "${split_test}")
    list(GET split_test 0 rel)
    list(GET split_test 1 test_name)
    oncvpsp_add_test(split split ${rel} ${test_name} split_test.sh
        $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp${rel}_ps>
        "${data}/${test_name}.dat" "@DIR@")
    oncvpsp_compares("${added_test}")
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name})
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name} AE)
endforeach ()

# Log derivatives from the all-electron and pseudopotential files (oncvpsp_ld.x): the run's with the
# defaults, each radius of a list as alone, and --no-log-derivatives (see scripts/ld_test.sh)
foreach (ld_test IN ITEMS "sr:14_Si" "nr:14_Si" "r:03_Li")
    string(REPLACE ":" ";" ld_test "${ld_test}")
    list(GET ld_test 0 rel)
    list(GET ld_test 1 test_name)
    oncvpsp_add_test(split ld ${rel} ${test_name} ld_test.sh
        $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp_ld>
        "${data}/${test_name}.dat" "@DIR@")
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name})
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name} AE)
endforeach ()

# The confined atom from the all-electron and pseudopotential files (oncvpsp_confine.x): the band edges
# of a Wigner-Seitz sphere against the all-electron atom's, each radius as alone, independent of the scan
# step, and the cutoffs in the sphere (see scripts/confine_test.sh)
foreach (confine_test IN ITEMS "sr:14_Si" "nr:14_Si" "r:03_Li")
    string(REPLACE ":" ";" confine_test "${confine_test}")
    list(GET confine_test 0 rel)
    list(GET confine_test 1 test_name)
    oncvpsp_add_test(split confine ${rel} ${test_name} confine_test.sh
        $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp_confine>
        "${data}/${test_name}.dat" "@DIR@")
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name} AE)
endforeach ()
# a core state's band (4f) in the window, an f wave function that gains a node inside the projectors,
# and a bound f ghost at -50 Ha
oncvpsp_add_test(split confine sr 60_Nd_GHOST confine_test.sh
    $<TARGET_FILE:oncvpspsr> $<TARGET_FILE:oncvpspsr_ae> $<TARGET_FILE:oncvpsp_confine>
    "${data}/60_Nd_GHOST.dat" "@DIR@" ghost)
oncvpsp_uses_run("${added_test}" sr 60_Nd_GHOST AE)

# Damaged or mismatched all-electron data files: oncvpsp_ps.x must stop with a message naming the
# problem, but accept differences within round-off (see scripts/ae_file_error_test.sh)
foreach (case IN ITEMS "byte_order" "zz" "fa" "mesh" "c_light" "zz_roundoff" "fa_roundoff"
                       "mesh_roundoff")
    oncvpsp_add_test(split ae_file sr 14_Si_${case} ae_file_error_test.sh
        $<TARGET_FILE:oncvpspsr> $<TARGET_FILE:oncvpspsr_ae> $<TARGET_FILE:oncvpspsr_ps>
        "${data}/14_Si.dat" ${case} "@DIR@")
    oncvpsp_uses_run("${added_test}" sr 14_Si AE)
    if (case MATCHES "_roundoff$")
        oncvpsp_compares("${added_test}")
        oncvpsp_uses_run("${added_test}" sr 14_Si)
    endif ()
endforeach ()

# All-electron SCF convergence on a mesh starting at 1e-8 bohr: see scripts/ae_scf_test.sh
oncvpsp_add_test(scf ae_scf nr 01_H ae_scf_test.sh
    $<TARGET_FILE:oncvpspnr_ae> "${data}/01_H_ae.dat" "@DIR@")

# A test configuration whose occupation does not round back from the ramp:
# see scripts/occupation_ramp_test.sh
oncvpsp_add_test(scf occupation_ramp sr 14_Si occupation_ramp_test.sh
    $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat" "@DIR@")

# All-electron atom with a 4f shell: see scripts/ae_4f_test.sh
oncvpsp_add_test(scf ae_4f_scf sr 70_Yb ae_4f_test.sh
    $<TARGET_FILE:oncvpspsr_ae> "${data}/70_Yb_ae.dat" scf "@DIR@")
oncvpsp_add_test(scf ae_4f_scf r 70_Yb ae_4f_test.sh
    $<TARGET_FILE:oncvpspr_ae> "${data}/70_Yb_ae.dat" scf "@DIR@")
oncvpsp_add_test(scf ae_4f_mesh sr 69_Tm ae_4f_test.sh
    $<TARGET_FILE:oncvpspsr_ae> "${data}/69_Tm_ae.dat" mesh "@DIR@")

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
        oncvpsp_add_test(nist nist ${rel} ${set}_${part_name} nist_test.py $<TARGET_FILE:oncvpsp${rel}_ae>
            ${set} "${data}/nist" "@DIR@" --z ${zs})
    endforeach ()
endforeach ()

# Log mesh input (the optional last line of the text input): see scripts/log_mesh_test.sh
foreach (mesh_test IN ITEMS "sr:14_Si:upf,split,order,amesh,rr1,rmax,read,adams_order"
                            "nr:29_Cu:split,split_differ"
                            "r:14_Si:split_differ")
    string(REPLACE ":" ";" mesh_test "${mesh_test}")
    list(GET mesh_test 0 rel)
    list(GET mesh_test 1 test_name)
    list(GET mesh_test 2 mesh_cases)
    string(REPLACE "," ";" mesh_cases "${mesh_cases}")
    if (ONCVPSP_WITH_TOML)
        list(APPEND mesh_cases "toml")
    endif ()
    foreach (case IN LISTS mesh_cases)
        oncvpsp_add_test(mesh log_mesh ${rel} ${test_name}_${case} log_mesh_test.sh
            $<TARGET_FILE:oncvpsp${rel}> $<TARGET_FILE:oncvpsp${rel}_ae> $<TARGET_FILE:oncvpsp${rel}_ps>
            "${data}/${test_name}.dat" ${case} "@DIR@"
            "${data}/toml/${test_name}.toml")
        if (case MATCHES "^(defaults|explicit|rerun|upf|split|toml|order)$")
            oncvpsp_compares("${added_test}")
        endif ()
    endforeach ()
endforeach ()

# UPF file on the log mesh (upf_mesh, the optional third value of the output grid line): see
# scripts/upf_mesh_test.sh
foreach (mesh_test IN ITEMS "sr:14_Si_UPF:negative,read" "r:14_Si_UPF:stride")
    string(REPLACE ":" ";" mesh_test "${mesh_test}")
    list(GET mesh_test 0 rel)
    list(GET mesh_test 1 test_name)
    list(GET mesh_test 2 mesh_cases)
    string(REPLACE "," ";" mesh_cases "${mesh_cases}")
    if (ONCVPSP_WITH_TOML)
        list(APPEND mesh_cases "toml")
    endif ()
    foreach (case IN LISTS mesh_cases)
        oncvpsp_add_test(mesh upf_mesh ${rel} ${test_name}_${case} upf_mesh_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" ${case}
            "@DIR@"
            "${data}/toml/${test_name}.toml")
        if (case MATCHES "^(defaults|log|stride|toml)$")
            oncvpsp_compares("${added_test}")
        endif ()
    endforeach ()
endforeach ()

# TOML input: every input in tests/data/toml (the TOML form of its tests/data/*.dat) must give the
# output of the text input, and the reader must stop with a message naming the problem for each of the
# changes in scripts/toml_error_test.sh
if (ONCVPSP_WITH_TOML)
    set(toml_tests
        "sr:03_Li" "sr:07_N_dl" "sr:07_N_opt" "sr:08_O" "sr:08_O_wigner" "sr:14_Si_bc" "sr:14_Si_bca"
        "sr:14_Si_bcw" "sr:20_Ca" "sr:32_Ge" "sr:34_Se" "sr:38_Sr_lxc" "sr:40_Zr" "sr:40_Zr_grid"
        "sr:60_Nd_GHOST" "sr:79_Au_lxc" "sr:07_N_flat"
        "r:07_N_dl" "r:14_Si_bc" "r:14_Si_hm" "r:20_Ca" "r:40_Zr_grid" "r:60_Nd_GHOST" "r:79_Au_lxc"
        "nr:14_Si")
    foreach (toml_test IN LISTS toml_tests)
        string(REPLACE ":" ";" toml_test "${toml_test}")
        list(GET toml_test 0 rel)
        list(GET toml_test 1 test_name)
        oncvpsp_add_test(toml toml ${rel} ${test_name} toml_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "${data}/toml/${test_name}.toml"
            "@DIR@")
        oncvpsp_compares("${added_test}")
        oncvpsp_uses_run("${added_test}" ${rel} ${test_name})
    endforeach ()

    foreach (case IN ITEMS
            "parse_error" "missing_table" "not_a_table" "missing_key" "missing_string"
            "integer_type" "number_type" "string_type" "atsym_length" "state_counts"
            "lmax_range" "array_length" "array_type" "not_an_array" "l_order"
            "icmod3_rcfact" "test_config_length" "test_config_count" "test_config_type"
            "integer_for_real" "extra_keys" "icmod6_default_rcfact"
            "icmod6_default_slack" "well_onset_range" "well_onset_default"
            "well_tail_value" "well_tail_default" "well_n_negative" "well_norm_value"
            "budget_cutoff_value"
            "budget_cutoff_default" "wave_regularization_value"
            "wave_regularization_default" "nan_real" "inf_array")
        oncvpsp_add_test(toml toml_input sr ${case} toml_error_test.sh
            $<TARGET_FILE:oncvpspsr> "${data}" ${case} "@DIR@")
        if (case MATCHES "^(integer_for_real|extra_keys|icmod6_default_rcfact|icmod6_default_slack|well_onset_default|well_tail_default|budget_cutoff_default|wave_regularization_default)$")
            oncvpsp_compares("${added_test}")
        endif ()
    endforeach ()
endif ()

# invariance_test(NAME REL INPUT ARGS...): scripts/invariance_test.py on oncvpsp<REL>.x and
# tests/data/toml/INPUT.toml with ARGS, labelled toml and REL; not added without TOML or, for REL r,
# without the fully-relativistic executables, and disabled without pyoncvpsp
function(invariance_test family rel case input)
    if (NOT ONCVPSP_WITH_TOML)
        return()
    endif ()
    oncvpsp_add_test(well ${family} ${rel} ${case} invariance_test.py $<TARGET_FILE:oncvpsp${rel}>
        "${data}/toml/${input}.toml" ${ARGN} --workdir "@DIR@")
    oncvpsp_compares("${added_test}")
endfunction()

# The onset of the well of the well states (well_onset in [pseudopotentials]): the same well states on
# two log meshes, and the input sets it (see scripts/well_onset_test.sh)
if (ONCVPSP_WITH_TOML)
    foreach (well_test IN ITEMS "sr:14_Si" "nr:08_O" "r:14_Si")
        string(REPLACE ":" ";" well_test "${well_test}")
        list(GET well_test 0 rel)
        list(GET well_test 1 test_name)
        oncvpsp_add_test(well well_onset ${rel} ${test_name} well_onset_test.sh
            $<TARGET_FILE:oncvpsp${rel}> "${data}/toml/${test_name}.toml"
            "@DIR@")
    endforeach ()
endif ()

# The free tail of the well states (well_tail = "free", doc/smoothness_plan.md, phase 8): nothing
# depends on the well, so the outputs are the same for any well_onset (bit for bit, measured; with
# the well's tail the kernels move by 5 % in 14_Si and 0.5 % in 08_O). Not the fr code, whose
# j = l-1/2 energy is the eigenvalue of the j = l+1/2 well (D-7c)
if (ONCVPSP_WITH_TOML)
    foreach (tail_test IN ITEMS "sr:14_Si" "nr:08_O")
        string(REPLACE ":" ";" tail_test "${tail_test}")
        list(GET tail_test 0 rel)
        list(GET tail_test 1 test_name)
        invariance_test(well_tail ${rel} ${test_name} ${test_name}
            pseudopotentials.well_onset 1.0306 1.0 1.09 2.5
            --set "pseudopotentials.well_tail=\"free\"" --groups kernel,vloc,ld,ecut,ghost --tol 1e-12)
    endforeach ()
endif ()

# The end of the log mesh (rmax) is not a box for the well states (well_mesh_m): 14_Si's d well state
# at 0.02 Ha, half-point radius 75 bohr, is solved on the mesh continued past rmax as far as it needs,
# so nothing moves with rmax but the SCF's round-off floor, ~1e-9 (bounded by the mesh's end at 45
# bohr, its kernel moved by 2e-2 and its 3e-5 Ha cutoff by 95 %)
if (ONCVPSP_WITH_TOML)
    foreach (rmax_test IN ITEMS "sr:14_Si" "r:14_Si" "nr:14_Si")
        string(REPLACE ":" ";" rmax_test "${rmax_test}")
        list(GET rmax_test 0 rel)
        list(GET rmax_test 1 test_name)
        invariance_test(well_rmax ${rel} ${test_name} ${test_name}
            log_mesh.rmax 45.0 60.0 80.0
            --set "pseudopotentials.ep=[0.0,0.0,0.02]" --groups kernel,vloc,tests,ld,ecut,ghost
            --tol 1e-7)
    endforeach ()
    # the confined well (well_tail = "confined"), whose d state at 0.02 Ha is held by a wall 1 Ha
    # high at 39 bohr
    invariance_test(well_rmax sr 14_Si_confined 14_Si log_mesh.rmax 45.0 60.0 80.0
        --set "pseudopotentials.ep=[0.0,0.0,0.02]" --set "pseudopotentials.well_tail=\"confined\""
        --groups kernel,vloc,tests,ld,ecut,ghost --tol 1e-7)
endif ()

# The normalization of the well states (well_norm) scales each state, which the kernels do not see:
# all space or inside rc give the same pseudopotential to round-off (6e-13 measured), for the well's
# tail and the confined well; only the convergence tables move (49 % at 1e-5 in 14_Si's d)
if (ONCVPSP_WITH_TOML)
    foreach (norm_test IN ITEMS "sr:well" "r:confined")
        string(REPLACE ":" ";" norm_test "${norm_test}")
        list(GET norm_test 0 rel)
        list(GET norm_test 1 tail)
        invariance_test(well_norm ${rel} 14_Si_${tail} 14_Si
            pseudopotentials.well_norm "\"space\"" "\"rc\""
            --set "pseudopotentials.well_tail=\"${tail}\"" --groups kernel,vloc,core,tests,ld,ghost
            --tol 1e-11)
    endforeach ()
endif ()

# PSML output: written into the output like the psp8 and upf files, with no files left in the working
# directory
foreach (psml_test IN ITEMS "sr:14_Si" "r:14_Si" "nr:08_O")
    string(REPLACE ":" ";" psml_test "${psml_test}")
    list(GET psml_test 0 rel)
    list(GET psml_test 1 test_name)
    oncvpsp_add_test(psml psml ${rel} ${test_name} psml_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "@DIR@")
endforeach ()

# Input echoes: the input echoed in the output, and in the psp8, UPF or PSML file, run again, give the
# same output (see scripts/echo_rerun_test.sh). The inputs cover rxpsh (14_Si_lder), icmod=4 with its
# search grid (40_Zr_grid), icmod=7 with slack (14_Si_hm), a Libxc iexc (83_Bi), rc on the grid
# (03_Li), iexc=7 with c_light on the first line (14_Si_c), ep_as_given, well_onset and well_tail =
# "free" on the lmax line (14_Si_wt), well_tail = "confined" with well_n on a projector line
# (14_Si_cw), well_norm on the lmax line (14_Si_wn), budget_cutoff on the icmod line (14_Si_bcw;
# "max_all" in 14_Si_hma), ecut_weight on a projector line (14_Si_bcwe) and
# wave_regularization on the lmax line (07_N_flat); each in sr or fr, the fr ones on inputs that run in seconds.
set(si_c "s/^    Si  14     3   2   3   psp8$/    Si  14     3   2   7   psp8  137.0359895/")
set(si_wt "14s/^    2$/    2    1    1.2    1/")
set(si_wn "14s/^    2$/    2    0    0.0    0    0.0    2/")
set(si_hma "s/^   7   0   1\\.0e-5   0\\.7$/   7   0   1.0e-5   0.7   2/")
set(si_ew "s/^    2    2    0\\.75$/    2    2    0.75    0    0.3/")
set(si_cw "14s/^    2$/    2    1    1.2    2/\ns/^    1    2    0\\.75$/    1    2    0.75    5/")
foreach (echo_test IN ITEMS
        "sr:14_Si_lder:psp8" "sr:40_Zr_grid:psp8" "sr:40_Zr_grid:upf"
        "sr:83_Bi:psml" "sr:03_Li:psml" "nr:14_Si:psml" "sr:14_Si_hm:psp8"
        "r:40_Zr_grid:psml" "r:14_Si_lder:psp8" "r:14_Si_c:upf" "r:14_Si_c:psml" "sr:14_Si_c:psp8"
        "sr:14_Si_wt:upf" "sr:14_Si_wn:psp8" "sr:14_Si_cw:psp8" "r:14_Si_cw:psml" "sr:14_Si_bcw:psml" "r:14_Si_bcw:psp8"
        "sr:14_Si_hma:upf" "r:14_Si_bcwe:psml" "sr:07_N_flat:upf")
    string(REPLACE ":" ";" echo_test "${echo_test}")
    list(GET echo_test 0 rel)
    list(GET echo_test 1 test_name)
    list(GET echo_test 2 format)
    set(edit "")
    set(input "${test_name}")
    if (test_name STREQUAL "14_Si_c")
        set(edit "${si_c}")
        set(input "14_Si")
    elseif (test_name STREQUAL "14_Si_wt")
        set(edit "${si_wt}")
        set(input "14_Si")
    elseif (test_name STREQUAL "14_Si_bcwe")
        set(edit "${si_ew}")
        set(input "14_Si_bcw")
    elseif (test_name STREQUAL "14_Si_hma")
        set(edit "${si_hma}")
        set(input "14_Si_hm")
    elseif (test_name STREQUAL "14_Si_wn")
        set(edit "${si_wn}")
        set(input "14_Si")
    elseif (test_name STREQUAL "14_Si_cw")
        set(edit "${si_cw}")
        set(input "14_Si")
    endif ()
    oncvpsp_add_test(echo echo_${format} ${rel} ${test_name} echo_rerun_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${input}.dat" ${format}
        "@DIR@" "${edit}")
    oncvpsp_compares("${added_test}")
endforeach ()

# icmod=6 with fcfact<=0 (the cutoff taken from the pseudo-waves): the psp8 file must still carry the
# model core (see scripts/psp8_core_test.sh)
oncvpsp_add_test(psp8 psp8_core sr 14_Si psp8_core_test.sh
    $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat" "s/^   1    0\\.25$/   6    0.0   1.0d-3/"
    "@DIR@")
oncvpsp_compares("${added_test}")
oncvpsp_add_test(psp8 psp8_core r 60_Nd_GHOST psp8_core_test.sh
    $<TARGET_FILE:oncvpspr> "${data}/60_Nd_GHOST.dat"
    "s/^    1   0\\.25000   0\\.00000$/    6    0.0   1.0d-3/" "@DIR@")
oncvpsp_compares("${added_test}")

# Text input from a pipe (as from a Python subprocess), which cannot be backspaced: the same output
# as from a file, and an input ending early reported with the line the reader stopped at
foreach (pipe_test IN ITEMS "sr:14_Si" "r:14_Si" "nr:08_O")
    string(REPLACE ":" ";" pipe_test "${pipe_test}")
    list(GET pipe_test 0 rel)
    list(GET pipe_test 1 test_name)
    oncvpsp_add_test(pipe pipe ${rel} ${test_name} pipe_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${data}/${test_name}.dat" "@DIR@")
    oncvpsp_compares("${added_test}")
    oncvpsp_uses_run("${added_test}" ${rel} ${test_name})
endforeach ()
oncvpsp_add_test(pipe pipe_truncated sr 14_Si pipe_truncated_test.sh $<TARGET_FILE:oncvpspsr> "${data}/14_Si.dat")

# Input errors: exit status 1 and the documented message (see scripts/error_test.sh)
foreach (rel IN ITEMS "sr" "nr" "r")
    if (rel STREQUAL "sr")
        set(program_name "oncvpsp")
    else ()
        set(program_name "oncvpsp${rel}")
    endif ()
    oncvpsp_add_test(errors error ${rel} missing_input error_test.sh
        $<TARGET_FILE:oncvpsp${rel}> "${program_name}: ERROR cannot open input file"
        "@DIR@")
endforeach ()
oncvpsp_add_test(errors error sr rxpsh_beyond_mesh error_test.sh
    $<TARGET_FILE:oncvpspsr> "run_phsft: ERROR rxpsh beyond radial mesh"
    "@DIR@"
    "${data}/14_Si.dat" "s/^   -2.0  2.0  0.02$/   -2.0  2.0  0.02  100.0/")
foreach (case IN ITEMS "inside:2.0" "at_rc:2.70")
    string(REPLACE ":" ";" case "${case}")
    list(GET case 0 where)
    list(GET case 1 radius)
    oncvpsp_add_test(errors error sr rxpsh_${where} error_test.sh
        $<TARGET_FILE:oncvpspsr> "check_data: ERROR rxpsh not beyond every core radius"
        "@DIR@"
        "${data}/14_Si_lder.dat" "s/^   -0\\.3  2\\.0  0\\.02  3\\.5$/   -0.3  2.0  0.02  ${radius}/")
endforeach ()
oncvpsp_add_test(errors error sr icmod7_slack_range error_test.sh
    $<TARGET_FILE:oncvpspsr> "test_data: must have 0.0<slack<=1.0 for icmod= 6 or 7"
    "@DIR@"
    "${data}/14_Si_hm.dat" "s/^   7   0   1\\.0e-5   0\\.7$/   7   0   1.0e-5   1.5/")
oncvpsp_add_test(errors error sr budget_cutoff_value error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  34"
    "@DIR@"
    "${data}/14_Si_bcw.dat" "s/^   6   0   1\\.0e-4   1\\.0   1$/   6   0   1.0e-4   1.0   3/")
oncvpsp_add_test(errors error sr wave_regularization_value error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  12"
    "@DIR@"
    "${data}/07_N_flat.dat" "s/      1e-08$/      -1e-08/")
oncvpsp_add_test(errors error sr iexc5_r2scan error_test.sh
    $<TARGET_FILE:oncvpspsr> "test_data: iexc=5 is r2SCAN, METAPSP only"
    "@DIR@"
    "${data}/14_Si.dat" "s/^    Si  14     3   2   3   psp8$/    Si  14     3   2   5   psp8/")
# well_n (the n of the first well state of an l) is for well_tail = "confined" only, and not below
# that state's n; well_tail is 0, 1 or 2
oncvpsp_add_test(errors error sr well_n_not_confined error_test.sh
    $<TARGET_FILE:oncvpspsr> "oncvpsp_run: ERROR well_n is for well_tail = \"confined\" only"
    "@DIR@"
    "${data}/14_Si.dat" "s/^    1    2    0\\.75$/    1    2    0.75    5/")
oncvpsp_add_test(errors error sr well_n_below error_test.sh
    $<TARGET_FILE:oncvpspsr> "below the first well state's n"
    "@DIR@"
    "${data}/14_Si.dat" "14s/^    2$/    2    0    0.0    2/\ns/^    1    2    0\\.75$/    1    2    0.75    3/")
# well_norm (the normalization of the well states) is 0, 1 or 2, and "space" (1) is not defined for
# the free tail
oncvpsp_add_test(errors error sr well_norm_free_space error_test.sh
    $<TARGET_FILE:oncvpspsr> "oncvpsp_run: ERROR well_norm = \"space\" (1) is not defined for well_tail = \"free\""
    "@DIR@"
    "${data}/14_Si.dat" "14s/^    2$/    2    0    0.0    1    0.0    1/")
oncvpsp_add_test(errors error sr well_norm_value error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  14"
    "@DIR@"
    "${data}/14_Si.dat" "14s/^    2$/    2    0    0.0    0    0.0    3/")
# ecut_weight is for a channel with no occupied state
oncvpsp_add_test(errors error sr ecut_weight_occupied error_test.sh
    $<TARGET_FILE:oncvpspsr> "ecut_weights: ERROR ecut_weight for l=1 is for a channel with no occupied state"
    "@DIR@"
    "${data}/14_Si.dat" "s/^    1    2    0\\.75$/    1    2    0.75    0    1.0/")
oncvpsp_add_test(errors error sr well_tail_value error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  14"
    "@DIR@"
    "${data}/14_Si.dat" "14s/^    2$/    2    0    0.0    3/")
# a real that is not finite (NaN or Inf, which list-directed input reads)
oncvpsp_add_test(errors error sr nan_real error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  33"
    "@DIR@"
    "${data}/14_Si.dat" "s/^   1    0\\.25$/   1    nan/")
oncvpsp_add_test(errors error sr inf_real error_test.sh
    $<TARGET_FILE:oncvpspsr> "Read ERROR, input data file line  18"
    "@DIR@"
    "${data}/14_Si.dat" "s/^    1    2\\.70 /    1    inf /")
oncvpsp_add_test(errors error sr icmod4_bad_grid error_test.sh
    $<TARGET_FILE:oncvpspsr> "check_data: ERROR icmod=4 needs 0<fcfact_min<=fcfact_max and fcfact_step>0"
    "@DIR@"
    "${data}/40_Zr.dat" "s/^    4   0\\.0$/    4   0.0   0.0   2.0   5.0   0.0   1.1   1.6   0.05/")
