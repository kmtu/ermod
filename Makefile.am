ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod_gromacs ermod_namd slvfe

dist_pkglibexec_SCRIPTS = vmdplugins/libexec/*.so

ersources = enganal.f90 engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 mpiproc.f90 engproc.f90 recpcal.f90 insertion.f90 setconf.f90 utility.f90 hash_real.c
EXTRA_DIST = parameters_er.example 
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/gromacs/*
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = -lfftw3
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_gromacs_SOURCES = $(ersources) vmdfio_if.f90 vmdfio.c
ermod_gromacs_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_gromacs_FCFLAGS = -DGROMACS -DVMDPLUGINS $(FFTFLAGS)
ermod_gromacs_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

extra_flags = 

# extra dependencies added
ermod_gromacs-engproc.$(OBJEXT): engproc.f90 ermod_gromacs-insertion.$(OBJEXT) ermod_gromacs-recpcal.$(OBJEXT) ermod_gromacs-recpcal.$(OBJEXT) ermod_gromacs-realcal_blk.$(OBJEXT)
ermod_gromacs-enganal.$(OBJEXT): enganal.f90 ermod_gromacs-engproc.$(OBJEXT) ermod_gromacs-setconf.$(OBJEXT) ermod_gromacs-vmdfio_if.$(OBJEXT)
ermod_gromacs-recpcal.$(OBJEXT): recpcal.f90 ermod_gromacs-spline.$(OBJEXT) ermod_gromacs-fft_iface.$(OBJEXT)
ermod_gromacs-insertion.$(OBJEXT): insertion.f90 ermod_gromacs-setconf.$(OBJEXT)
ermod_gromacs-setconf.$(OBJEXT): setconf.f90 ermod_gromacs-engmain.$(OBJEXT) ermod_gromacs-mpiproc.$(OBJEXT)

ermod_namd_SOURCES = $(ersources)
ermod_namd_FCFLAGS = -DNAMD $(FFTFLAGS)
ermod_namd_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS)

ermod_namd-engproc.$(OBJEXT): engproc.f90 param_eng ermod_namd-insertion.$(OBJEXT) ermod_namd-recpcal.$(OBJEXT) ermod_namd-recpcal.$(OBJEXT) ermod_namd-realcal_blk.$(OBJEXT)
ermod_namd-enganal.$(OBJEXT): enganal.f90 ermod_namd-engproc.$(OBJEXT) ermod_namd-setconf.$(OBJEXT)
ermod_namd-recpcal.$(OBJEXT): recpcal.f90 ermod_namd-spline.$(OBJEXT) ermod_namd-fft_iface.$(OBJEXT)
ermod_namd-insertion.$(OBJEXT): insertion.f90 ermod_namd-setconf.$(OBJEXT)
ermod_namd-setconf.$(OBJEXT): setconf.f90 ermod_namd-engmain.$(OBJEXT) ermod_namd-mpiproc.$(OBJEXT)




slvfe_SOURCES = slvfe.f
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
