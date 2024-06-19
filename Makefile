########################################################################
#                                                                      #
#     FINFLO® Makefile                                                 #
#                                                                      #
#     Elomatic Consulting & Engineering                                #
#                                                                      #
#     Vaisalantie 2             Tel.  +358-40-7751010                  #
#     FIN-02130 Espoo           E-mail esa.salminen@elomatic.com       #
#     Finland                   URL    http://www.finflo.fi            #
#                                                                      #
########################################################################

SHELL = /bin/sh

FC    = mpif90
CC    = mpicc

MpCCI = NO

OBJ   = base3c.o multi.o turb3.o rotdia.o reynol.o fscal.o messp.o \
        tgas.o ugas.o outp3.o state.o bound.o chim.o sec.o imps4.o \
        scimp.o earsm.o charac.o surpro.o fresur.o sixdof.o transm.o \
	fsitools.o

ARCH  = $(shell getconf LONG_BIT)

ifeq "$(MpCCI)" "YES"
        OBJ       += mpccitools.o
        COBJ       = mpcciadapter.o
        MPCCI_HOME = $(shell mpcci home)
        MPCCI_ARCH = $(shell mpcci arch)
        HEADERS    = -I. -I$(MPCCI_HOME)/include 
        LIBRARIES  = $(MPCCI_HOME)/lib/$(MPCCI_ARCH)/libmpcci-64.a
else
        OBJ       += emptympccitools.o
        COBJ       =
        HEADERS    = 
        LIBRARIES  =
endif


GRID  = griddeform.o gridmodules.o

FIN   = main.o ns3c.o res3c.o imps3.o gluebm.o precor.o ns3fs.o \
        $(OBJ) $(COBJ) $(GRID)


ifeq "$(OPT)" ""
  OPT := 3
endif

OPTMZ := $(shell if [ $(OPT) -le 0 ] ; then echo 0 ; \
               elif [ $(OPT) -eq 1 ] ; then echo 1 ; \
               elif [ $(OPT) -eq 2 ] ; then echo 2 ; \
               elif [ $(OPT) -ge 3 ] ; then echo 3 ; fi)


ifeq "$(FORTRAN)" ""
  FORTRAN := ifort
endif

ifeq "$(FORTRAN)" "ifort"
        IFCFLAGS = -r8 -O$(OPTMZ) -xHost -w # -check bounds -traceback
#        IFCFLAGS = -r8 -O$(OPTMZ) -xHost -assume buffered_io -w  # -check bounds -traceback
        FFLAGS   = $(IFCFLAGS)
        CFLAGS   = -Wall 
        GFLAGS   = -O3 -w
        LFLAGS   = 
endif

ifeq "$(FORTRAN)" "gfortran"
#        IFCFLAGS = -O$(OPTMZ) -fdefault-real-8 -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal
        IFCFLAGS = -O$(OPTMZ) -fdefault-real-8 -w -fallow-argument-mismatch -ffpe-summary=none
#        IFCFLAGS = -O$(OPTMZ) -fdefault-real-8 -w  -fbounds-check -fbacktrace
#        IFCFLAGS = -O$(OPTMZ) -fdefault-real-8 -fdefault-double-8 -w  # -fbounds-check -fbacktrace
        FFLAGS   = $(IFCFLAGS)
        LFLAGS   = 
        GFLAGS   = -O0 -fcray-pointer -w
endif


install: all

all: finflo

clean:
	rm -f *.o *.mod finflo

finflo: $(FIN)
ifeq "$(MpCCI)" "YES"
	@rm -f emptympccitools.o
else
	@rm -f mpccitools.o
endif
	$(FC) $(LFLAGS) $(FIN) $(LIBRARIES) -o finflo

main.o ns3c.o ns3fs.o: %.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

res3c.o imps3.o gluebm.o precor.o: %.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJ): %.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(COBJ): %.o: %.c
	$(CC) -c $(HEADERS) $(CFLAGS) $< -o $@

$(GRID): %.o: %.f
	$(FC) -c $(GFLAGS) $< -o $@

# Module dependencies

main.o ns3c.o reynol.o outp3.o earsm.o turb3.o state.o surpro.o: charac.o
res3c.o imps3.o imps4.o gluebm.o precor.o rotdia.o base3c.o ns3fs.o: charac.o
messp.o multi.o sixdof.o transm.o chim.o mpccitools.o fscal.o: charac.o
fresur.o bound.o fsitools.o: charac.o
griddeform.o fresur.o ns3fs.o main.o ns3c.o: gridmodules.o
