ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod_gromacs ermod_namd slvfe

dist_pkglibexec_SCRIPTS = vmdplugins/libexec/*.so

ersources = engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 mpiproc.f90
EXTRA_DIST = enganal.f insertion.f setconf.f stranal.f param_eng parameters_er.example gromacsfio.c
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/gromacs/*
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = -lfftw3
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_gromacs_SOURCES = $(ersources) vmdfio.c
ermod_gromacs_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_gromacs_FCFLAGS = $(FFTFLAGS)
ermod_gromacs_LDADD = ermod_gromacs-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm

extra_flags = 

# extra dependencies added
ermod_gromacs-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DGROMACS -DVMDPLUGINS $(extra_flags) -c -o $@ $<

ermod_namd_SOURCES = $(ersources)
ermod_namd_FCFLAGS = $(FFTFLAGS)
ermod_namd_LDADD = ermod_namd-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS)

ermod_namd-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DNAMD -c -o $@ $<


slvfe_SOURCES = slvfe.f
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
