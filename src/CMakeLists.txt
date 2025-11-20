###########################################################
# ONCVPSP library
###########################################################
set(src_oncvpsp
    input_text_m.f90
    check_data_m.f90
    adams_m.f90
    dpnint_m.f90
    der2exc_m.f90
    eresid_m.f90
    exc_libxc_list_m.F90
    exc_internal_m.f90
    exc_off_pbe_m.f90
    teter_m.f90
    dirac_m.f90
    schroedinger_m.f90
    vout_m.f90
    wellstate_m.f90
    model_core_m.f90
    dft_m.f90
    sbf_m.f90
    vsl_m.f90
    vploc_m.f90
    vkb_m.f90
    thomas_fermi_m.f90
    vpinteg_m.f90
    wf_rc_der_m.f90
    vrel_m.f90
    psp8_m.f90
    upf_m.f90
    psml_m.f90
    test_ghosts_m.f90
    test_log_der_m.f90
    test_pseudization_m.f90
    test_config_m.f90
    gnuplot_m.f90
    plot_m.f90
    renormalize_m.f90
    overlap_m.f90
)
# Append Libxc interface or stub
if(Libxc_FOUND AND ONCVPSP_WITH_Libxc)
    list(APPEND src_oncvpsp exc_libxc_m.F90 exc_libxc_parser_m.f90 exc_libxc_xc_messages_m.f90)
else()
    list(APPEND src_oncvpsp exc_libxc_m_stub.f90)
endif()
# Append TOML input parser if enabled
if(ONCVPSP_WITH_TOML)
    list(APPEND src_oncvpsp input_toml_m.f90)
endif()
# Create ONCVPSP library
add_library(oncvpsp STATIC ${src_oncvpsp})
target_link_libraries(oncvpsp wxml ${ONCVPSP_DEPENDENCY_LIBRARIES})
if(ONCVPSP_WITH_TOML)
    target_compile_definitions(oncvpsp PRIVATE WITH_TOML)
endif()

###########################################################
# oncvpsp.x
###########################################################
add_executable(oncvpspsr oncvpsp.F90)
target_link_libraries(oncvpspsr PRIVATE oncvpsp)
set_target_properties(oncvpspsr
    PROPERTIES
    OUTPUT_NAME oncvpsp.x
)

##########################################################
# oncvpspnr.x
##########################################################
add_executable(oncvpspnr oncvpsp_nr.F90)
target_link_libraries(oncvpspnr PRIVATE oncvpsp)
set_target_properties(oncvpspnr
    PROPERTIES
    OUTPUT_NAME oncvpspnr.x
)

##########################################################
# oncvpspr.x
##########################################################
add_executable(oncvpspfr oncvpsp_r.F90)
target_link_libraries(oncvpspfr PRIVATE oncvpsp)
set_target_properties(oncvpspfr
    PROPERTIES
    OUTPUT_NAME oncvpspr.x
)

###########################################################
# Subdirectories
###########################################################
add_subdirectory(xmlf90-wxml)
