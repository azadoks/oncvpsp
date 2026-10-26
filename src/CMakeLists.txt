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
    const_basis_m.f90
    constants_m.f90
    der2exc_m.f90
    derivs_m.f90
    diagout.f90
    vloc_cutoff.f90
    dpnint_m.f90
    eresid.f90
    free_tail_m.f90
    exc_off_pbe_m.f
    exc_relx_m.f90
    excvwn_m.f90
    exchdl_m.f90
    excpzca_m.f90
    excwig_m.f90
    inward_start.f90
    inward2_m.f90
    cinf_taper_m.f90
    ionic_tail.f90
    ldiracfb.f90
    ldiracfs.f90
    lschfb.f90
    lschfs.f90
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
    modcore_m.f90
    ode_rc_der_m.f90
    ae_delta.f90
    oncvpsp_main.F90
    optimize.f90
    psatom.f90
    pspot.f90
    qroots_m.f90
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
    sbf_basis_con_m.f90
    sbf_rc_der_m.f90
    sbf8_m.f90
    sbf_m.f90
    sbf_overlap_m.f90
    scatterstate.f90
    sr_flux_m.f90
    dirac_m.f90
    sbf_q_m.f90
    phsft_m.f90
    ws_basis_m.f90
    pieces_m.f90
    outmesh_m.f90
    vploc_m.f90
    qlinalg_m.f90
    ugrid_m.f90
    quad_m.f90
    xc_dispatch_m.f90
    boundsolve_m.f90
    kappa_m.f90
    projbasis_m.f90
    main_setup_m.F90
    ae_projectors_m.f90
    ps_valence_m.f90
    cli_m.f90
    ae_ps_pair_m.f90
    confine_summary_m.f90
    confine_ws_m.f90
    confine_wall_m.f90
    outwf_m.F90
    sratom.f90
    tfapot_m.f90
    vkboutwf.f90
    wf_rc_der_m.f90
    vout_m.f90
    wellstate.f90
    well_mesh_m.f90
    czt_m.f90
    input_echo_m.f90
    input_text_m.f90
    # Non-/scalar-relativistic files
    ae_file_m.f90
    ps_file_m.f90
    fpovlp_m.f90
    gnu_script.f90
    linout.f90
    upfout.f90
    # Fully-relativistic files
    sr_so_r.f90
)
if (ONCVPSP_WITH_Libxc)
    list(APPEND src_oncvpsp functionals.F90 exc_libxc_m.f90)
else ()
    list(APPEND src_oncvpsp exc_libxc_stub_m.f90)
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
target_compile_definitions(oncvpsp PUBLIC $<$<BOOL:${ONCVPSP_WITH_TOML}>:WITH_TOML>
    $<$<BOOL:${ONCVPSP_BLAS_DOT}>:WITH_BLAS_DOT>)
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

# Log derivatives at any radii and energies from the all-electron file and the pseudopotential file of
# --ps-file (see oncvpsp_ld.f90 and ps_file_m.f90), for every relativistic treatment
oncvpsp_executable(oncvpsp_ld oncvpsp_ld.f90 oncvpsp_ld.x)

# The confined atom: the band edges of a Wigner-Seitz sphere and the cutoffs of the pseudopotential's states
# in it, from the same two files (see oncvpsp_confine.f90 and ws_basis_m.f90)
oncvpsp_executable(oncvpsp_confine oncvpsp_confine.f90 oncvpsp_confine.x)

set(ONCVPSP_EXECUTABLES ${ONCVPSP_EXECUTABLES} PARENT_SCOPE)
