include_directories(SYSTEM)

###########################################################
# oncvpsp.x
###########################################################
add_executable(oncvpsp_exe oncvpsp.f90)
add_dependencies(oncvpsp_exe oncvpsp wxml)
target_link_libraries(oncvpsp_exe oncvpsp wxml Libxc::xcf03 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(oncvpsp_exe
    PROPERTIES
        OUTPUT_NAME oncvpsp.x
)

##########################################################
# oncvpspnr.x
##########################################################
add_executable(oncvpsp_nr_exe oncvpsp_nr.f90)
add_dependencies(oncvpsp_nr_exe oncvpsp)
target_link_libraries(oncvpsp_nr_exe oncvpsp wxml Libxc::xcf03 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(oncvpsp_nr_exe
    PROPERTIES
        OUTPUT_NAME oncvpspnr.x
)

##########################################################
# oncvpspr.x
##########################################################
add_executable(oncvpsp_r_exe oncvpsp_r.f90)
add_dependencies(oncvpsp_r_exe oncvpsp)
target_link_libraries(oncvpsp_r_exe oncvpsp wxml Libxc::xcf03 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(oncvpsp_r_exe
    PROPERTIES
        OUTPUT_NAME oncvpspr.x
)
