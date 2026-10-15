# Apply a patch to the current directory unless it is already applied, so that the patch step of a
# dependency (which runs again whenever its download step does) can be repeated, also after the patch
# has changed.
#
# Run in script mode, as the PATCH_COMMAND of CPMAddPackage:
#   cmake -DPATCH=<file.patch> -P ApplyPatch.cmake

if (NOT PATCH)
    message(FATAL_ERROR "ApplyPatch.cmake: PATCH must be defined")
endif ()
find_package(Git REQUIRED QUIET)

execute_process(
    COMMAND "${GIT_EXECUTABLE}" apply --reverse --check "${PATCH}"
    RESULT_VARIABLE applied
    OUTPUT_QUIET ERROR_QUIET
)
if (applied EQUAL 0)
    message(STATUS "${PATCH} is already applied")
    return()
endif ()
# Another version of the patch may have been applied: a git checkout (of the dependency's pinned
# version) is restored first
execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --show-prefix
    RESULT_VARIABLE not_git
    OUTPUT_VARIABLE prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if (not_git EQUAL 0 AND prefix STREQUAL "")
    execute_process(COMMAND "${GIT_EXECUTABLE}" checkout -- .)
endif ()
execute_process(
    COMMAND "${GIT_EXECUTABLE}" apply "${PATCH}"
    RESULT_VARIABLE result
)
if (NOT result EQUAL 0)
    message(FATAL_ERROR "Could not apply ${PATCH}")
endif ()
