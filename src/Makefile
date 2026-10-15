# Makefile for ONCVPSP
#
# Copyright (c) 1989-2016 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
#

include ../make.inc

OBJS_ALL =	build_info_m.o check_data.o const_basis.o dpnint.o eresid.o exchdl.o \
		exc_off_pbe.o excpzca.o excwig.o input_text_m.o ldiracfb.o ldiracfs.o \
		lschfb.o lschfs.o lschkb.o lschpb.o lschpse.o lschps.o \
		lschvkbbe.o lschvkbb.o lschvkbs.o modcore.o optimize.o \
		pspot.o qroots.o run_optimize.o sbf8.o sbf_basis.o \
		sbf_rc_der.o tfapot.o vkboutwf.o vkbphsft.o vout.o vpinteg.o \
		vploc.o wf_rc_der.o modcore2.o gg1cc.o gp1cc.o gpp1cc.o \
		der2exc.o modcore3.o modcore_budget.o sbf_basis_con.o run_ghosts.o \
		lschpsbar.o vrel.o ae_file_m.o config_ae_m.o sratom.o

OBJS_NRSR =	fphsft.o fpovlp.o gnu_script.o linout.o psatom.o \
		relatom.o run_config.o run_diag.o run_phsft.o run_plot.o \
		run_vkb.o upfout.o wellstate.o

OBJS_SR =	m_psmlout.o oncvpsp.o

OBJS_NR =	m_psmlout.o oncvpsp_nr.o

OBJS_SR_AE =	m_psmlout.o oncvpsp_ae.o

OBJS_SR_PS =	m_psmlout.o oncvpsp_ps.o

OBJS_NR_AE =	m_psmlout.o oncvpspnr_ae.o

OBJS_NR_PS =	m_psmlout.o oncvpspnr_ps.o

OBJS_RLIB =	fphsft_r.o fpovlp.o gnu_script_r.o linout_r.o \
		psatom_r.o relatom.o renorm_r.o run_config_r.o run_diag_r.o \
		run_diag_sr_so_r.o run_phsft_r.o run_plot_r.o run_vkb_r.o \
		sr_so_r.o upfout_r.o wellstate_r.o

OBJS_R =	m_psmlout.o oncvpsp_r.o $(OBJS_RLIB)

OBJS_R_AE =	m_psmlout.o oncvpspr_ae.o $(OBJS_RLIB)

OBJS_R_PS =	m_psmlout.o oncvpspr_ps.o $(OBJS_RLIB)

OBJS_OPSP =     $(OBJS_SR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPNR = 	$(OBJS_NR) $(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_OPSPR = 	$(OBJS_R) $(OBJS_ALL) $(OBJS_LIBXC)

OBJS_SPLIT =	$(OBJS_NRSR) $(OBJS_ALL) $(OBJS_LIBXC)

all:	oncvpsp oncvpspnr oncvpspr oncvpsp_ae oncvpsp_ps oncvpspnr_ae oncvpspnr_ps \
	oncvpspr_ae oncvpspr_ps

WXML_LIB=xmlf90-wxml/libwxml.a
wxml:
	(cd xmlf90-wxml; $(MAKE) module)
#
m_psmlout.o: wxml m_libxc_list.o
DRIVERS_NRSR =	oncvpsp.o oncvpsp_nr.o oncvpsp_ae.o oncvpsp_ps.o \
		oncvpspnr_ae.o oncvpspnr_ps.o
DRIVERS_R =	oncvpsp_r.o oncvpspr_ae.o oncvpspr_ps.o
$(DRIVERS_NRSR) $(DRIVERS_R): m_psmlout.o input_text_m.o
$(DRIVERS_NRSR) $(DRIVERS_R) run_config.o run_config_r.o ae_file_m.o: config_ae_m.o
$(DRIVERS_NRSR) $(DRIVERS_R): ae_file_m.o
oncvpsp_nr.o oncvpsp_ae.o oncvpsp_ps.o oncvpspnr_ae.o oncvpspnr_ps.o: oncvpsp.F90
oncvpspr_ae.o oncvpspr_ps.o: oncvpsp_r.F90
$(DRIVERS_NRSR) $(DRIVERS_R): build_info_m.o
build_info_m.o: build_info.inc

# Git commit, source state and build time, regenerated at every build (see ../cmake/BuildInfo.cmake)
build_info.inc: FORCE
	cmake -DSOURCE_DIR=.. -DOUTPUT=build_info.inc -P ../cmake/BuildInfo.cmake
FORCE:
.PHONY: FORCE

ldiracfb.o ldiracfs.o lschfb.o lschfs.o lschkb.o lschpb.o lschps.o \
lschpsbar.o lschpse.o lschvkbb.o vout.o: adams.inc
INCLUDE_WXML=-I./xmlf90-wxml



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

oncvpspr_ae:	$(OBJS_R_AE) $(OBJS_ALL) $(OBJS_LIBXC)
		$(FLINKER) -o oncvpspr_ae.x $(OBJS_R_AE) $(OBJS_ALL) $(OBJS_LIBXC) $(LIBS) $(WXML_LIB)

oncvpspr_ps:	$(OBJS_R_PS) $(OBJS_ALL) $(OBJS_LIBXC)
		$(FLINKER) -o oncvpspr_ps.x $(OBJS_R_PS) $(OBJS_ALL) $(OBJS_LIBXC) $(LIBS) $(WXML_LIB)

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

