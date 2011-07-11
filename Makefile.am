ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod_gromacs ermod_namd slvfe

dist_pkglibexec_SCRIPTS = vmdplugins/libexec/*.so

ersources = enganal.f90 engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 mpiproc.f90 engproc.f90 recpcal.f90
EXTRA_DIST = insertion.f setconf.f param_eng parameters_er.example gromacsfio.c
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
ermod_gromacs_FCFLAGS = -DGROMACS -DVMDPLUGINS $(FFTFLAGS)
ermod_gromacs_LDADD = ermod_gromacs-insertion.$(OBJEXT) ermod_gromacs-setconf.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm

extra_flags = 

# extra dependencies added
ermod_gromacs-engproc.$(OBJEXT): engproc.f90 param_eng ermod_gromacs-insertion.$(OBJEXT) ermod_gromacs-recpcal.$(OBJEXT)

ermod_gromacs-insertion.$(OBJEXT): insertion.f param_eng ermod_gromacs-setconf.$(OBJEXT)
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DGROMACS -DVMDPLUGINS $(extra_flags) -c -o $@ $<

ermod_gromacs-setconf.$(OBJEXT): setconf.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DGROMACS -DVMDPLUGINS $(extra_flags) -c -o $@ $<

ermod_namd_SOURCES = $(ersources)
ermod_namd_FCFLAGS = -DNAMD $(FFTFLAGS)
ermod_namd_LDADD = ermod_namd-insertion.$(OBJEXT) ermod_namd-setconf.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS)

ermod_namd-engproc.$(OBJEXT): engproc.f90 param_eng ermod_namd-insertion.$(OBJEXT) ermod_namd-recpcal.$(OBJEXT)

ermod_namd-insertion.$(OBJEXT): insertion.f param_eng ermod_namd-setconf.$(OBJEXT)
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DNAMD -c -o $@ $<

ermod_namd-setconf.$(OBJEXT): setconf.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) $(FFTFLAGS) -DNAMD -c -o $@ $<



slvfe_SOURCES = slvfe.f
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
