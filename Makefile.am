ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod_gromacs ermod_namd ermod_charmm slvfe

ersources = engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90
EXTRA_DIST = enganal.f insertion.f setconf.f stranal.f param_eng parameters_er.example
EXTRA_DIST += tools/NAMD/*

ermod_gromacs_SOURCES = $(ersources)
ermod_gromacs_LDADD = ermod_gromacs-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS)

# extra dependencies added
ermod_gromacs-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DGROMACS -c -o $@ $<

ermod_namd_SOURCES = $(ersources)
ermod_namd_LDADD = ermod_namd-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS)

ermod_namd-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DNAMD -c -o $@ $<

#
ermod_charmm_SOURCES = $(ersources)
ermod_charmm_LDADD = ermod_charmm-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) 

ermod_charmm-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DCHARMM -c -o $@ $<

slvfe_SOURCES = slvfe.f
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
