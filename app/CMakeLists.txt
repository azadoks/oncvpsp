include_directories(SYSTEM)

###########################################################
# oncvpsp.x
###########################################################
add_executable(oncvpsp_exe oncvpsp.f90)
add_dependencies(oncvpsp_exe oncvpsp wxml)
target_link_libraries(oncvpsp_exe oncvpsp wxml Libxc::xcf03 ${LAPACKBLAS_LIBRARIES})
set_target_properties(oncvpsp_exe
    PROPERTIES
        OUTPUT_NAME oncvpsp.x
)

###########################################################
# oncvpsp_nr.x
###########################################################
# add_executable(oncvpsp_nr_exe oncvpsp_nr.f90)
# add_dependencies(oncvpsp_nr_exe oncvpsp)
# target_link_libraries(oncvpsp_nr_exe
#     oncvpsp
# )
# set_target_properties(oncvpsp_nr_exe
#     PROPERTIES
#         OUTPUT_NAME oncvpsp_nr.x
# )

###########################################################
# oncvpsp_r.x
###########################################################
# add_executable(oncvpsp_r_exe oncvpsp_r.f90)
# add_dependencies(oncvpsp_r_exe oncvpsp)
# target_link_libraries(oncvpsp_r_exe
#     oncvpsp
# )
# set_target_properties(oncvpsp_r_exe
#     PROPERTIES
#         OUTPUT_NAME oncvpsp_r.x
# )
