# The Python that compares outputs in the tests (scripts/compare_output.py), which needs pyoncvpsp
# (https://github.com/azadoks/pyoncvpsp), to parse the outputs and their psp8 and PSML files, and
# upf-tools, to parse their UPF files.
#
# ONCVPSP_TEST_PYTHON, if set, is that Python. Otherwise the configuration makes a virtual
# environment in the build tree (test-venv/) with ONCVPSP_PYONCVPSP and ONCVPSP_UPF_TOOLS, with uv if
# it is on PATH, else with the venv module and pip of python3, once: it is made again only when they
# change. A pyoncvpsp source directory is installed editable, so that its changes need no
# reconfiguration; one from a repository is not updated until test-venv/ is removed.
#
# Sets ONCVPSP_COMPARE_PYTHON to the Python, or to "" (with a warning) if there is none, and
# ONCVPSP_TEST_PYTHON_DESCRIPTION for the summary.

set(ONCVPSP_TEST_PYTHON "" CACHE FILEPATH
    "Python with pyoncvpsp and upf-tools for the output comparisons of the tests (empty: make a venv in the build tree)")
set(ONCVPSP_PYONCVPSP "git+https://github.com/azadoks/pyoncvpsp@oncvpsp-tests" CACHE STRING
    "pyoncvpsp for that venv: a source directory (installed editable) or a pip requirement")
set(ONCVPSP_UPF_TOOLS "upf-tools==0.2.0" CACHE STRING "upf-tools for that venv: a pip requirement")

# what compare_output.py imports
set(_oncvpsp_python_check
    "import upf_tools; from pyoncvpsp.io import OncvpspTextParser, parse_psml, parse_psp8")

# _oncvpsp_python_works(PYTHON RESULT): whether PYTHON imports what compare_output.py needs
function(_oncvpsp_python_works python result)
    set(${result} FALSE PARENT_SCOPE)
    if (EXISTS "${python}")
        execute_process(COMMAND "${python}" -c "${_oncvpsp_python_check}"
            RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
        if (status EQUAL 0)
            set(${result} TRUE PARENT_SCOPE)
        endif ()
    endif ()
endfunction()

set(ONCVPSP_COMPARE_PYTHON "")
if (ONCVPSP_TEST_PYTHON)
    _oncvpsp_python_works("${ONCVPSP_TEST_PYTHON}" works)
    if (works)
        set(ONCVPSP_COMPARE_PYTHON "${ONCVPSP_TEST_PYTHON}")
        set(ONCVPSP_TEST_PYTHON_DESCRIPTION "${ONCVPSP_TEST_PYTHON}")
    else ()
        set(_oncvpsp_python_problem
            "ONCVPSP_TEST_PYTHON=${ONCVPSP_TEST_PYTHON} cannot import pyoncvpsp (with parse_psp8 and parse_psml) and upf_tools")
    endif ()
else ()
    set(venv "${CMAKE_BINARY_DIR}/test-venv")
    set(venv_python "${venv}/bin/python")
    set(stamp "${venv}/oncvpsp-requirements.txt")
    set(pyoncvpsp "${ONCVPSP_PYONCVPSP}")
    if (IS_DIRECTORY "${pyoncvpsp}")
        get_filename_component(pyoncvpsp "${pyoncvpsp}" ABSOLUTE)
        set(pyoncvpsp_args -e "${pyoncvpsp}")
    else ()
        set(pyoncvpsp_args "${pyoncvpsp}")
    endif ()
    set(requirements "${pyoncvpsp_args};${ONCVPSP_UPF_TOOLS}")
    set(installed "")
    if (EXISTS "${stamp}")
        file(READ "${stamp}" installed)
    endif ()
    if (NOT EXISTS "${venv_python}" OR NOT installed STREQUAL requirements)
        find_program(ONCVPSP_UV uv)
        message(STATUS "Making the Python of the output comparisons in ${venv}: ${ONCVPSP_PYONCVPSP}, ${ONCVPSP_UPF_TOOLS}")
        file(REMOVE_RECURSE "${venv}")
        if (ONCVPSP_UV)
            execute_process(COMMAND "${ONCVPSP_UV}" venv --quiet "${venv}"
                RESULT_VARIABLE status ERROR_VARIABLE error)
            if (status EQUAL 0)
                execute_process(
                    COMMAND "${ONCVPSP_UV}" pip install --quiet --python "${venv_python}"
                            ${pyoncvpsp_args} "${ONCVPSP_UPF_TOOLS}"
                    RESULT_VARIABLE status ERROR_VARIABLE error)
            endif ()
        else ()
            execute_process(COMMAND "${Python3_EXECUTABLE}" -m venv "${venv}"
                RESULT_VARIABLE status ERROR_VARIABLE error)
            if (status EQUAL 0)
                execute_process(
                    COMMAND "${venv_python}" -m pip install --quiet ${pyoncvpsp_args} "${ONCVPSP_UPF_TOOLS}"
                    RESULT_VARIABLE status ERROR_VARIABLE error)
            endif ()
        endif ()
        if (status EQUAL 0)
            file(WRITE "${stamp}" "${requirements}")
        else ()
            string(STRIP "${error}" error)
            set(_oncvpsp_python_problem "making ${venv} failed: ${error}")
        endif ()
    endif ()
    _oncvpsp_python_works("${venv_python}" works)
    if (works)
        set(ONCVPSP_COMPARE_PYTHON "${venv_python}")
        set(ONCVPSP_TEST_PYTHON_DESCRIPTION "${venv} (${ONCVPSP_PYONCVPSP}, ${ONCVPSP_UPF_TOOLS})")
    elseif (NOT DEFINED _oncvpsp_python_problem)
        set(_oncvpsp_python_problem
            "${venv} cannot import pyoncvpsp (with parse_psp8 and parse_psml) and upf_tools: is ONCVPSP_PYONCVPSP=${ONCVPSP_PYONCVPSP} recent enough?")
    endif ()
endif ()

if (NOT ONCVPSP_COMPARE_PYTHON)
    message(WARNING
        "No Python for the output comparisons of the tests: ${_oncvpsp_python_problem}. "
        "The tests that compare outputs (scripts/compare_output.py) are disabled. Set "
        "ONCVPSP_TEST_PYTHON to a Python with pyoncvpsp and upf-tools, or ONCVPSP_PYONCVPSP to a "
        "pyoncvpsp source directory or pip requirement (see doc/developer_setup.md).")
    set(ONCVPSP_TEST_PYTHON_DESCRIPTION "none: the output comparisons are disabled")
endif ()
