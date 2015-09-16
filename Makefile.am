ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod slvfe
nobase_pkgdata_SCRIPTS = \
	tools/NAMD/gen_structure tools/NAMD/gen_input \
	tools/gromacs/gen_structure tools/gromacs/gen_input \
	tools/AMBER/gen_structure tools/AMBER/gen_input \
	tools/MARBLE/gen_structure tools/MARBLE/gen_input \
	tools/LAMMPS/gen_structure tools/LAMMPS/gen_input \
	tools/CHARMM/gen_structure tools/CHARMM/gen_input

dist_pkglibexec_SCRIPTS = vmdplugins/libexec/*.so

ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 mpiproc.F90 engproc.F90 recpcal.F90 insertion.F90 setconf.F90 trajectory.F90 utility.F90 bestfit.F90

EXTRA_DIST = parameters_er.example parameters_fe.example
EXTRA_DIST += extensions/*
EXTRA_DIST += trajectory-formats/*
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/gromacs/*
EXTRA_DIST += tools/AMBER/*
EXTRA_DIST += tools/MARBLE/*
EXTRA_DIST += tools/CHARMM/*
EXTRA_DIST += tools/LAMMPS/*
EXTRA_DIST += tools/mol_dissection
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = $(FFTWLIBS)
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_SOURCES = $(ersources)
ermod_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_FCFLAGS = $(VMDPLUGINFLAGS) $(FFTFLAGS)
ermod_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

if USE_VMDPLUGIN
ermod_SOURCES += vmdfio.c
endif

if USE_TRANSFER
ermod_FCFLAGS += -DHAVE_TRANSFER
else
ermod_SOURCES += hash_real.c
endif

extra_flags = 

# module dependencies
ermod-enganal.$(OBJEXT): enganal.F90 ermod-engproc.$(OBJEXT) ermod-setconf.$(OBJEXT) ermod-insertion.$(OBJEXT) ermod-trajectory.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-engmain.$(OBJEXT)
ermod-engproc.$(OBJEXT): engproc.F90 ermod-insertion.$(OBJEXT) ermod-realcal.$(OBJEXT) ermod-recpcal.$(OBJEXT) ermod-bestfit.$(OBJEXT) ermod-utility.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-engmain.$(OBJEXT)
ermod-realcal.$(OBJEXT): realcal.F90 ermod-engmain.$(OBJEXT)
ermod-recpcal.$(OBJEXT): recpcal.F90 ermod-spline.$(OBJEXT) ermod-fft_iface.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-engmain.$(OBJEXT)
ermod-insertion.$(OBJEXT): insertion.F90 ermod-setconf.$(OBJEXT) ermod-bestfit.$(OBJEXT) ermod-trajectory.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-engmain.$(OBJEXT)
ermod-setconf.$(OBJEXT): setconf.F90 ermod-utility.$(OBJEXT) ermod-trajectory.$(OBJEXT) ermod-mpiproc.$(OBJEXT) ermod-engmain.$(OBJEXT)
ermod-bestfit.$(OBJEXT): bestfit.F90 ermod-mpiproc.$(OBJEXT)
ermod-trajectory.$(OBJEXT): trajectory.F90 ermod-utility.$(OBJEXT)
ermod-utility.$(OBJEXT): utility.F90 ermod-engmain.$(OBJEXT)
ermod-mpiproc.$(OBJEXT): mpiproc.F90 ermod-engmain.$(OBJEXT)

slvfe_SOURCES = slvfe.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC_CANONICAL) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
