# Makefile for ONCVPSP
#
# Copyright (c) 1989-2016 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

ifeq ($(wildcard ../make.inc),)
$(error ../make.inc not found. Copy ../make.inc.example to ../make.inc and fill it in)
endif
include ../make.inc

ifeq ($(strip $(F90)),)
$(error F90 is not set. Fill in ../make.inc (see ../make.inc.example))
endif
ifeq ($(strip $(OBJS_LIBXC)),)
$(error OBJS_LIBXC is not set. Fill in ../make.inc (see ../make.inc.example))
endif

OBJS_ALL =	aeatom.o boundstate.o build_info_m.o check_data.o const_basis.o constants_m.o dpnint.o eresid.o exchdl.o \
		exc_off_pbe.o excpzca.o excwig.o input_text_m.o ldiracfb.o ldiracfs.o \
		lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
		lschvkbbe.o lschvkbb.o lschvkbs.o modcore.o optimize.o \
		pspot.o qroots.o run_optimize.o sbf8.o \
		sbf_rc_der.o tfapot.o vkboutwf.o vkbphsft.o vout.o vpinteg.o \
		vploc.o wf_rc_der.o modcore2.o gg1cc.o gp1cc.o gpp1cc.o \
		der2exc.o modcore3.o modcore_budget.o sbf_basis_con.o run_ghosts.o \
		lschpsbar.o vrel.o ae_file_m.o config_ae_m.o sratom.o input_echo_m.o \
		relatom.o fphsft.o run_phsft.o scatstate.o wellstate.o \
		psatom.o run_config.o run_diag.o renorm_r.o run_plot.o \
		run_vkb.o fpovlp.o gnu_script.o linout.o upfout.o \
		run_diag_sr_so_r.o sr_so_r.o \
		psmlout_m.o oncvpsp_main.o

OBJS_SR =	oncvpsp.o

OBJS_NR =	oncvpsp_nr.o

OBJS_SR_AE =	oncvpsp_ae.o

OBJS_SR_PS =	oncvpsp_ps.o

OBJS_NR_AE =	oncvpspnr_ae.o

OBJS_NR_PS =	oncvpspnr_ps.o

OBJS_R =	oncvpsp_r.o

OBJS_R_AE =	oncvpspr_ae.o

OBJS_R_PS =	oncvpspr_ps.o

OBJS_OPSP =     $(OBJS_SR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPNR = 	$(OBJS_NR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPR = 	$(OBJS_R) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_SPLIT =	$(OBJS_ALL) $(OBJS_LIBXC)

all:	oncvpsp oncvpspnr oncvpspr oncvpsp_ae oncvpsp_ps oncvpspnr_ae oncvpspnr_ps \
	oncvpspr_ae oncvpspr_ps

# WXML_LIB and INCLUDE_WXML, for xmlf90, are set in make.inc
psmlout_m.o: m_libxc_list.o
oncvpsp_main.o: psmlout_m.o input_text_m.o
oncvpsp_main.o run_config.o: config_ae_m.o
oncvpsp_main.o: ae_file_m.o
oncvpsp_main.o: build_info_m.o
oncvpsp_main.o linout.o upfout.o psmlout_m.o: input_echo_m.o
build_info_m.o: build_info.inc

# Git commit, source state and build time, regenerated at every build (see ../cmake/BuildInfo.cmake)
build_info.inc: FORCE
	cmake -DSOURCE_DIR=.. -DOUTPUT=build_info.inc -P ../cmake/BuildInfo.cmake
FORCE:
.PHONY: FORCE

ldiracfb.o ldiracfs.o lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
lschpsbar.o lschvkbb.o vkboutwf.o vout.o: adams.inc

exc_libxc.o: functionals.o

# Every routine uses constants_m
$(filter-out constants_m.o,$(sort $(OBJS_ALL) \
	$(OBJS_LIBXC))): constants_m.o



oncvpsp:	$(OBJS_OPSP)
		$(FLINKER) -o oncvpsp.x $(OBJS_OPSP) $(LIBS) $(WXML_LIB)

oncvpspnr:	$(OBJS_OPSPNR)
		$(FLINKER) -o oncvpspnr.x $(OBJS_OPSPNR) $(LIBS) $(WXML_LIB)

oncvpspr:	$(OBJS_OPSPR)
		$(FLINKER) -o oncvpspr.x $(OBJS_OPSPR) $(LIBS) $(WXML_LIB)

oncvpsp_ae:	$(OBJS_SR_AE) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpsp_ae.x $(OBJS_SR_AE) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

oncvpsp_ps:	$(OBJS_SR_PS) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpsp_ps.x $(OBJS_SR_PS) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

oncvpspnr_ae:	$(OBJS_NR_AE) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpspnr_ae.x $(OBJS_NR_AE) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

oncvpspnr_ps:	$(OBJS_NR_PS) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpspnr_ps.x $(OBJS_NR_PS) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

oncvpspr_ae:	$(OBJS_R_AE) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpspr_ae.x $(OBJS_R_AE) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

oncvpspr_ps:	$(OBJS_R_PS) $(OBJS_SPLIT)
		$(FLINKER) -o oncvpspr_ps.x $(OBJS_R_PS) $(OBJS_SPLIT) $(LIBS) $(WXML_LIB)

.SUFFIXES:
.SUFFIXES: .c .o .f .f90 .F90

#	$(FCCPP) $(FCCPPFLAGS) $*.F90  > $*_cpp.f90
.F90.o:
	$(F90) $(FCCPPFLAGS) -cpp $(FFLAGS) $(INC) $(INCLUDE_WXML) -o $*.o -c $*.F90

.f90.o:
	$(F90) $(FFLAGS) $(INC) $(INCLUDE_WXML) -c $*.f90

.f.o:
	$(F90) $(FFLAGS) $(INC) $(INCLUDE_WXML) -c $*.f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c

clean:
	 /bin/rm -f *.o  *.d  *.mod *.log *_cpp.f90 build_info.inc
