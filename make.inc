# System-dependent makefile options for ONCVPSP
# This must be carefully edited before executing "make" in src
#
# Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
 
##### Edit the following lines to correspond to your compilers ####

F77        = /opt/homebrew/bin/gfortran	
F90        = /opt/homebrew/bin/gfortran	
CC         = /opt/homebrew/bin/gcc	
FCCPP      = /opt/homebrew/bin/cpp-13

FLINKER     = $(F90)

FCCPPFLAGS = -ansi 

##### Edit the following optimization flags for your system ####


FFLAGS     = -O3 -fallow-argument-mismatch -ffast-math -funroll-loops
CFLAGS     = -O3

##### Edit the following LAPACK and BLAS library paths for your system ####

LIBS = -L/Users/mverstra/CODES/LAPACK/ -llapack -lblas

##### Edit the following for to use libxc if available #####

OBJS_LIBXC =	exc_libxc_stub.o

# To build oncvpsp with libxc, uncomment the following 3 lines and edit
# the paths to point to your libxc library and include directories
# make clean in src before rebuilding after changing this

#LIBS += -L/Users/mverstra/CODES/LIBXC/5.1.7/lib -lxcf90 -lxc
#FFLAGS += -I/Users/mverstra/CODES/LIBXC/5.1.7/include
#FCCPPFLAGS += -I/Users/mverstra/CODES/LIBXC/5.1.7/include
#OBJS_LIBXC =    functionals.o exc_libxc.o

LIBS += -L/Users/mverstra/CODES/LIBXC/6.0.0/lib -lxcf03 -lxc
FFLAGS += -I/Users/mverstra/CODES/LIBXC/6.0.0/include
FCCPPFLAGS += -I/Users/mverstra/CODES/LIBXC/6.0.0/include
OBJS_LIBXC =	functionals.o exc_libxc.o 

# xmlf90 1.6.3 (https://gitlab.com/siesta-project/libraries/xmlf90), whose
# wxml writes the PSML output, with cmake/patches/xmlf90-1.6.3.patch
# applied (the CMake build fetches and patches it itself):
#   git clone --branch 1.6.3 https://gitlab.com/siesta-project/libraries/xmlf90.git
#   cd xmlf90 && git apply <oncvpsp>/cmake/patches/xmlf90-1.6.3.patch
#   cmake -B build -DCMAKE_INSTALL_PREFIX=$(XMLF90_ROOT) && cmake --build build --target install
XMLF90_ROOT = /usr/local/xmlf90
# (its modules are installed in a directory named after the compiler, which
# its pkg-config file gives)
XMLF90_PC = PKG_CONFIG_PATH=$(XMLF90_ROOT)/lib/pkgconfig pkg-config xmlf90
WXML_LIB = $(shell $(XMLF90_PC) --libs)
INCLUDE_WXML = $(shell $(XMLF90_PC) --cflags)
