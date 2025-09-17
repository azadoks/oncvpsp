#!/bin/bash
#runs ONCVPSP with the command-line argument <prefix> and the graphics
#which review the results
#_r is appended to the prefix of the output file
#uses the fully-relativistic all-electron atom calculation

PREFIX=/path/to/repository

INFILE=$1.dat

OUTFILE=$1_r.out

GNUFILE=$$.scr

PLOTFILE=$1_r.plot

TEMP=$$.tmp

$PREFIX/build/bin/src/oncvpspr.x <$INFILE >$OUTFILE  #Edit if your executable is
                                            #in another directory

awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}};\
	/DATA FOR PLOTTING/{out=1}' $OUTFILE >$PLOTFILE

awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}};\
	/GNUSCRIPT/{out=1}' $OUTFILE >$TEMP

sed -e s/t1/$PLOTFILE/ $TEMP | sed -e s/t2/$1_r/ >$GNUFILE

if [ "$2" != "-np" ]
then
	gnuplot $GNUFILE
fi

rm  $GNUFILE $TEMP $PLOTFILE
