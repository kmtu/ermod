ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = er_gromacs er_namd er_charmm

ersources = engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90
extra_DIST = enganal.f insertion.f setconf.f stranal.f param_eng
CLEANFILES = er_gromacs-enganal.$(OBJEXT) er_namd-enganal.$(OBJEXT) er_charmm-enganal.$(OBJEXT)

er_gromacs_SOURCES = $(ersources)
er_gromacs_LDADD = er_gromacs-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS)

# extra dependencies added
er_gromacs-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DGROMACS -c -o $@ $<

er_namd_SOURCES = $(ersources)
er_namd_LDADD = er_namd-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS)

er_namd-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DNAMD -c -o $@ $<

#
er_charmm_SOURCES = $(ersources)
er_charmm_LDADD = er_charmm-enganal.$(OBJEXT) $(LAPACK_LIBS) $(BLAS_LIBS) 

er_charmm-enganal.$(OBJEXT): enganal.f insertion.f setconf.f stranal.f param_eng
	$(F77) $(AM_FFLAGS) $(FFLAGS) -DCHARMM -c -o $@ $<
