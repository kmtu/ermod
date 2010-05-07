ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = er_gromacs er_namd er_charmm

ersources = engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 enganal.f

er_gromacs_SOURCES = $(ersources)
er_gromacs_FFLAGS = -DGROMACS
er_gromacs_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)

er_namd_SOURCES = $(ersources)
er_namd_FFLAGS = -DNAMD
er_namd_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)

er_charmm_SOURCES = $(ersources)
er_charmm_FFLAGS = -DCHARMM
er_charmm_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) 




