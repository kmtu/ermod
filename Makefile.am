ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod slvfe

dist_pkglibexec_SCRIPTS = vmdplugins/libexec/*.so

ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 mpiproc.F90 engproc.F90 recpcal.F90 insertion.F90 setconf.F90 trajectory.F90 utility.F90 hash_real.c

EXTRA_DIST = parameters_er.example 
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/gromacs/*
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
EXTRA_DIST += realcal.F90.php
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = -lfftw3
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_SOURCES = $(ersources) vmdfio_if.F90 vmdfio.c
ermod_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_FCFLAGS = -DGROMACS -DVMDPLUGINS $(FFTFLAGS)
ermod_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

extra_flags = 

# module dependencies
ermod-engproc.$(OBJEXT): engproc.F90 ermod-insertion.$(OBJEXT) ermod-recpcal.$(OBJEXT) ermod-recpcal.$(OBJEXT) ermod-realcal.$(OBJEXT) ermod-utility.$(OBJEXT)
ermod-enganal.$(OBJEXT): enganal.F90 ermod-engproc.$(OBJEXT) ermod-setconf.$(OBJEXT) ermod-vmdfio_if.$(OBJEXT)
ermod-recpcal.$(OBJEXT): recpcal.F90 ermod-spline.$(OBJEXT) ermod-fft_iface.$(OBJEXT)
ermod-insertion.$(OBJEXT): insertion.F90 ermod-setconf.$(OBJEXT)
ermod-setconf.$(OBJEXT): setconf.F90 ermod-engmain.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-trajectory.$(OBJEXT)

slvfe_SOURCES = slvfe.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)

.phony: generate
generate: 
	php5 realcal.F90.php > realcal.F90
