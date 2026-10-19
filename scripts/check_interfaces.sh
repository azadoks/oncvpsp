#!/bin/bash
#
# Check every call in ONCVPSP against the definition of the routine it calls.
#
# Most routines are external procedures with no explicit interface, so the
# compiler can't check their calls. This script concatenates the drivers and
# the external procedures into one file per executable, so that gfortran sees
# both sides of every call and reports argument mismatches as errors.
#
# Usage: scripts/check_interfaces.sh [module-dir]
#
# module-dir holds the .mod files of the ONCVPSP modules (input_text_m,
# m_psmlout, ...) from a CMake build; the default is build/modules. Build
# with CMake first.
#
# Not checked: the modules themselves (their procedures already have
# explicit interfaces), exc_off_pbe.f (fixed form), and the libxc interface
# (exc_libxc_stub.f90 is used in its place, so calls into it are checked
# against the stub's signature).

set -euo pipefail

FC=${FC:-gfortran}
top=$(cd "$(dirname "$0")/.." && pwd)
moddir=${1:-$top/build/modules}
moddir=$(cd "$moddir" && pwd)

if [ ! -f "$moddir/m_psmlout.mod" ]; then
    echo "check_interfaces: ERROR no ONCVPSP modules in $moddir (build with CMake first)"
    exit 1
fi

# Files in the src_oncvpsp list of src/CMakeLists.txt, minus the modules and
# the fixed-form file
sources=$(awk '/^set\(src_oncvpsp/{on=1; next} on && /^\)/{exit} on' \
              "$top/src/CMakeLists.txt" |
          sed -e 's/#.*//' -e 's/[[:space:]]//g' |
          grep -E '\.[fF]90$' |
          grep -vE '^(input_text_m\.f90|m_psmlout\.f90|m_libxc_list\.F90)$' | sort -u)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

status=0
for driver in oncvpsp.F90:RELATIVISTIC=1 oncvpsp.F90:RELATIVISTIC=0 oncvpsp_r.F90:; do
    file=${driver%%:*}
    define=${driver#*:}
    all="$work/all.F90"
    : > "$all"
    for src in $sources exc_libxc_stub.f90; do
        cat "$top/src/$src" >> "$all"
        echo >> "$all"
    done
    cat "$top/src/$file" >> "$all"

    label="$file${define:+ ($define)}"
    if $FC -fsyntax-only -cpp ${define:+-D$define} -ffree-line-length-none \
           -J "$work" -I "$moddir" "$all" > "$work/log" 2>&1; then
        echo "check_interfaces: $label: OK"
    else
        echo "check_interfaces: $label: FAILED"
        # Point the messages back at the original files
        awk -v all="$all" -v top="$top/src/" '
            FNR==NR { start[++n]=$1; name[n]=$2; next }
            {
                if (match($0, /all\.F90:[0-9]+/)) {
                    line=substr($0, RSTART+8, RLENGTH-8)+0
                    for (i=n; i>=1; i--) if (line>=start[i]) break
                    sub(/[^ ]*all\.F90:[0-9]+/, name[i] ":" line-start[i]+1)
                }
                print
            }' <(line=1
                 for src in $sources exc_libxc_stub.f90 $file; do
                     echo "$line $src"
                     line=$((line + $(wc -l < "$top/src/$src") + 1))
                 done) "$work/log"
        status=1
    fi
done
exit $status
