###########################################################
# ONCVPSP library
###########################################################
set(src_oncvpsp
    # Common files
    ae_scf.f90
    aeatom.f90
    boundstate.f90
    build_info_m.f90
    check_data.f90
    config_ae_m.f90
    const_basis.f90
    constants_m.f90
    der2exc.f90
    derivs.f90
    dpnint.f90
    eresid.f90
    exc_off_pbe.f
    exc_relx.f90
    excvwn.f90
    exchdl.f90
    excpzca.f90
    excwig.f90
    fphsft.f90
    gg1cc.f90
    gp1cc.f90
    gpp1cc.f90
    inward_start.f90
    ionic_tail.f90
    ldiracfb.f90
    ldiracfs.f90
    lschfb.f90
    lschfs.f90
    lschkb.f90
    lschpb.f90
    lschps.f90
    lschpsbar.f90
    lschvkbb.f90
    lschvkbbe.f90
    lschvkbs.f90
    logmesh_m.f90
    match_m.f90
    m_libxc_list.F90
    psmlout_m.f90
    modcore.f90
    modcore2.f90
    modcore3.f90
    modcore_budget.f90
    ode_rc_der.f90
    ae_delta.f90
    oncvpsp_main.F90
    optimize.f90
    psatom.f90
    pspot.f90
    qroots.f90
    relatom.f90
    relativity_m.f90
    renorm_r.f90
    run_config.f90
    run_diag.f90
    run_ghosts.f90
    run_optimize.f90
    run_phsft.f90
    run_plot.f90
    run_vkb.f90
    sbf_basis_con.f90
    sbf_rc_der.f90
    sbf8.f90
    sbf_overlap.f90
    scatterstate.f90
    sr_flux_m.f90
    sbf_q_m.f90
    pieces_m.f90
    vploc_m.f90
    qlinalg_m.f90
    ugrid_m.f90
    sratom.f90
    tfapot.f90
    vkboutwf.f90
    wf_rc_der.f90
    vkbphsft.f90
    vout.f90
    vpinteg.f90
    wellstate.f90
    input_echo_m.f90
    input_text_m.f90
    # Non-/scalar-relativistic files
    ae_file_m.f90
    fpovlp.f90
    gnu_script.f90
    linout.f90
    upfout.f90
    # Fully-relativistic files
    run_diag_sr_so_r.f90
    sr_so_r.f90
)
if (ONCVPSP_WITH_Libxc)
    list(APPEND src_oncvpsp functionals.F90 exc_libxc.f90)
else ()
    list(APPEND src_oncvpsp exc_libxc_stub.f90)
endif ()
if (ONCVPSP_WITH_TOML)
    list(APPEND src_oncvpsp input_toml_m.f90)
endif ()

# Git commit, source state and build time, regenerated at every build (see cmake/BuildInfo.cmake)
add_custom_target(oncvpsp_build_info
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${PROJECT_SOURCE_DIR} -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/build_info.inc
            -P ${PROJECT_SOURCE_DIR}/cmake/BuildInfo.cmake
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/build_info.inc
    COMMENT "Checking git commit and source state for build_info.inc"
    VERBATIM
)

add_library(oncvpsp STATIC ${src_oncvpsp})
add_library(oncvpsp::oncvpsp ALIAS oncvpsp)
add_dependencies(oncvpsp oncvpsp_build_info)
set_target_properties(oncvpsp PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(oncvpsp
    PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/modules"
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
)
target_compile_definitions(oncvpsp PUBLIC $<$<BOOL:${ONCVPSP_WITH_TOML}>:WITH_TOML>)
target_link_libraries(oncvpsp
    PUBLIC
        oncvpsp::lapack
        xmlf90-lib
        $<$<BOOL:${ONCVPSP_WITH_Libxc}>:oncvpsp::libxc>
        $<$<BOOL:${ONCVPSP_WITH_TOML}>:toml-f::toml-f>
    PRIVATE
        oncvpsp::compile_options
)

###########################################################
# Executables
###########################################################
set(ONCVPSP_EXECUTABLES "")
function(oncvpsp_executable target source output_name)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE oncvpsp::oncvpsp oncvpsp::compile_options)
    set_target_properties(${target} PROPERTIES OUTPUT_NAME ${output_name})
    list(APPEND ONCVPSP_EXECUTABLES ${target})
    set(ONCVPSP_EXECUTABLES ${ONCVPSP_EXECUTABLES} PARENT_SCOPE)
endfunction()

# Scalar-relativistic, non-relativistic and fully-relativistic drivers
oncvpsp_executable(oncvpspsr oncvpsp.F90 oncvpsp.x)
oncvpsp_executable(oncvpspnr oncvpsp_nr.F90 oncvpspnr.x)
if (ONCVPSP_WITH_FR)
    oncvpsp_executable(oncvpspr oncvpsp_r.F90 oncvpspr.x)
endif ()

# All-electron and pseudopotential stages of the three drivers, which pass the all-electron atom through a
# file (see oncvpsp_main.F90 and ae_file_m.f90)
oncvpsp_executable(oncvpspsr_ae oncvpsp_ae.F90 oncvpsp_ae.x)
oncvpsp_executable(oncvpspsr_ps oncvpsp_ps.F90 oncvpsp_ps.x)
oncvpsp_executable(oncvpspnr_ae oncvpspnr_ae.F90 oncvpspnr_ae.x)
oncvpsp_executable(oncvpspnr_ps oncvpspnr_ps.F90 oncvpspnr_ps.x)
if (ONCVPSP_WITH_FR)
    oncvpsp_executable(oncvpspr_ae oncvpspr_ae.F90 oncvpspr_ae.x)
    oncvpsp_executable(oncvpspr_ps oncvpspr_ps.F90 oncvpspr_ps.x)
endif ()

set(ONCVPSP_EXECUTABLES ${ONCVPSP_EXECUTABLES} PARENT_SCOPE)
