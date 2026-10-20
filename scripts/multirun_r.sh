#!/bin/bash
#runs ONCVPSP for each command-line argument <prefix>.dat
#launches up to NPROC jobs at a time
#uses the fully-relativistic all-electron atom calculation
#does not display plots
#

NPROC=4

# the executables: $ONCVPSP_BIN_DIR, or build/bin of the CMake build of
# this tree
BIN_DIR=${ONCVPSP_BIN_DIR:-$(cd "$(dirname "$0")/.." && pwd)/build/bin}

nn=0


for ii do

	nn=$(( $nn + 1))

	OUTFILE=`echo $ii | sed -e s/.dat// `_r.out

	$BIN_DIR/oncvpspr.x <$ii >$OUTFILE &


       if [ $(($nn%$NPROC)) == 0 ]
                then
                wait
        fi

done
wait

echo '*** ALL DONE ***'
