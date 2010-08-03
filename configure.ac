#                                               -*- Autoconf -*-

AC_PREREQ([2.59])
#AC_INIT(enerepr, 0.0.%rev%, Matubayasi Laboratory)
AC_INIT(enerepr, 0.0.91, Matubayasi Laboratory)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_MACRO_DIR(m4)

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC
AC_PROG_F77([ifort gfortran xlf f90 g77 xlf90])
AC_PROG_FC([ifort gfortran xlf f77 g95 xlf77])

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran 77)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

AC_ARG_WITH([mkl],
  AS_HELP_STRING([--with-mkl[[=(MKL version)]]],
    [Automatically set environmental variables for MKL; for a specific version, use this option like --with-mkl=10.2.3]))

AS_IF([test "x$with_mkl" != "x"],
  [case x`uname -m` in # (
  xx86_64) :
  mkl_plat=em64t
  ;; # (
  xi386) :
  mkl_plat=32
  ;; # (
  *) :
  ;;
  esac
  AS_IF([test "x$with_mkl" = "xyes"],
    AS_IF([test "x$MKLROOT" = "x"],
      [AC_MSG_FAILURE([could not get MKL version; set MKLROOT or use --with-mkl=(version number)])],
      mklroot="$MKLROOT"),
    mklroot="/opt/intel/mkl/$with_mkl")
  with_blas="-L$mklroot/lib/$mkl_plat -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread"
  with_lapack="-L$mklroot/lib/$mkl_plat -lmkl_lapack"])

# due to bug in ACX_BLAS/LAPACK
AC_LANG_PUSH(C)
ACX_BLAS
ACX_LAPACK([lapack_exist=yes],
[lapack_exist=no])
AC_LANG_POP(C)

# Fortran/FORTRAN compatibility issues
AC_DEFUN([MY_F77_EXTEND_SOURCE],[
AC_CACHE_CHECK([how to extend column width],
my_f77_extend_source_cv_flag,
[AC_LANG_PUSH(Fortran 77)
save_FFLAGS=$FFLAGS
for flag in "" "-ffixed-line-length-132" "-extend-source"; do
FFLAGS="$save_FFLAGS $flag"
  AC_COMPILE_IFELSE([      program main
                                                                                test
      end program], [my_check_f77_extend=no], [my_check_f77_extend=yes])
  if test $my_check_f77_extend = yes; then
    break
  fi
done
if test $my_check_f77_extend = no; then
  AC_MSG_ERROR([unable to find compiler flag for extended source (longer column length)])
fi
FFLAGS="$save_FFLAGS $flag"
echo -n $flag
AC_LANG_POP(Fortran 77)
])])

MY_F77_EXTEND_SOURCE

AC_DEFUN([MY_FC_REAL_DOUBLE],[
AC_CACHE_CHECK([how to read real as real(8) as default],
my_fc_real_double_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
for flag in "" "-fdefault-real-8" "-autodouble"; do
FCFLAGS="$save_FCFLAGS $flag"
  AC_RUN_IFELSE(
    [AC_LANG_PROGRAM([], [dnl
      real :: a
      real(8) :: b
      a = 1.000000001d0
      b = 1.000000001d0
      if(dble(a) /= b) then
        stop 1
      end if
])], [my_check_realasdouble=yes], [my_check_realasdouble=no])
  if test $my_check_realasdouble = yes; then
    break
  fi
done
if test $my_check_realasdouble = no; then
  AC_MSG_ERROR([unable to find compiler flag to read real as real(8)])
fi
FCFLAGS="$save_FCFLAGS $flag"
FFLAGS="$FFLAGS $flag"
if test x$flag = x; then
  echo -n "(none)"
else
  echo -n $flag
fi
AC_LANG_POP(Fortran)
])])

MY_FC_REAL_DOUBLE

# set to f90 as default
AC_FC_SRCEXT(f90)

# select FFT library
AC_ARG_WITH([fft],
	AS_HELP_STRING([--with-fft=(fftw|mkl)],
		[specify which FFT library to be used]),
	[],
	[with_fft=no]
)

# Search for fftw3 - libfftw3.a
AC_CHECK_LIB([fftw3], [dfftw_execute], [have_fftw3=yes], [have_fftw3=no])

AS_IF([test "x$with_fft" == xfftw], 
            [AS_IF([test "x$have_fftw3" == xno],
                         [AC_MSG_FAILURE(
                                 [--with-fft=fftw is specified, but could not find fftw3])])])
AS_IF([test "x$with_fft" == xno],
	    [AS_IF([test "x$have_fftw3" == xyes],
	    		 [with_fft=fftw])
             AS_IF([test "x$with_mkl" != "x"],
                         [with_fft=mkl])])

AS_IF([test "x$with_fft" != xfftw && test "x$with_fft" != xmkl],
	    [AC_MSG_FAILURE(
		[invalid --with-fft=<fft-type> specification, please check with --help])])


# check whether mkl_dfti.f90 is accesible
# FIXME: This may fail to detect with gfortran 4.4.0-4.4.3
AS_IF([test "x$with_fft" = xmkl], [dnl
  AC_FC_SRCEXT(F90)
  AC_MSG_CHECKING([mkl_dfti.f90 is accessible])
  AC_COMPILE_IFELSE([#include "mkl_dfti.f90"
program main 
end program], [echo yes], [echo no
  AC_MSG_FAILURE([could not find mkl_dfti.f90 in include path])])
  AC_FC_SRCEXT(f90)])

AC_ARG_ENABLE([mpi],
	AS_HELP_STRING([--disable-mpi],
		[disable MPI even if exists]),
	AS_IF([test x$enable_mpi = xno],
		     [use_mpi=no]),
	[])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# FIXME: check cpp flags
FFLAGS="$FFLAGS -g -cpp -Dtrjctry"
FCFLAGS="$FCFLAGS -g -cpp -Dtrjctry"

AC_ARG_ENABLE([perf],
	AS_HELP_STRING([--enable-perf],
		[(internal use only)]),
	AS_IF([test x$enable_perf = xyes],
		     [FFLAGS="$FFLAGS -DPERF"
	              FCFLAGS="$FCFLAGS -DPERF"
                      LDFLAGS="$LDFLAGS -lprofiler"]))

AS_IF([test "x$use_mpi" = "xyes"],
	    [mpi_result="enabled"
	     F77="$MPIF77"
	     FC="$MPIFC"],
	    [mpi_result="disabled"
	     FFLAGS="$FFLAGS -DnoMPI"
	     FCFLAGS="$FCFLAGS -DnoMPI"])

# optimization flags
AC_DEFUN([MY_FC_OPT_ACCEPT_NATIVE],[
AC_CACHE_CHECK([fortran compiler accepts -march=native],
my_fc_opt_accept_native_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
flag="-march=native"
FCFLAGS="$save_FCFLAGS $flag"
AC_COMPILE_IFELSE([      program main
      end program], [my_fc_opt_accept_native=yes], [my_fc_opt_accept_native=no])
echo -n $my_fc_opt_accept_native
FCFLAGS=$save_FCFLAGS
AC_LANG_POP(Fortran)
])])

AC_ARG_ENABLE([opt],
	AS_HELP_STRING([--disable-opt],
		[disable automatic edit to optimization flags]))

AS_IF([test x$enable_opt = xyes || test x$enable_opt = x],
  [FFLAGS="$FFLAGS"
   FCFLAGS="$FCFLAGS"
   AS_IF([test x$ac_compiler_gnu = xyes],
        [optflag="-O3 -ffast-math"
         MY_FC_OPT_ACCEPT_NATIVE
         AS_IF([test x$my_fc_opt_accept_native = xyes],
                optflag="$optflag -march=native")],
        [case x$F77 in #(
           xmpiifort|xifort)
             optflag="-O3 -no-prec-div -xHost"
             ;; #(
           *)
             optflag="-O"
             ;;
         esac])
   FFLAGS="$FFLAGS $optflag"
   FCFLAGS="$FCFLAGS $optflag"])

# AM_CONDITIONAL([MPI], [test x$use_mpi = xyes])

case "x$with_fft" in #(
 xfftw) :
 FFLAGS="$FFLAGS -DFFTW"
 FCFLAGS="$FCFLAGS -DFFTW"
 LDFLAGS="$LDFLAGS -lfftw3"
 ;; #(
 xmkl) :
 FFLAGS="$FFLAGS -DMKL"
 FCFLAGS="$FCFLAGS -DMKL"
 ;; #(
 *)
 ;;
esac

AS_IF([test "x$lapack_exist" != xyes],
	    [AC_MSG_ERROR(
		[dnl
LAPACK library was not found.
Try specifying --with-lapack=<lib> options or set FCFLAGS/FFLAGS/FLIBS environment variables to fix the problem.])])

lapack_result="$LAPACK_LIBS $BLAS_LIBS"


AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_MSG_RESULT([dnl
LAPACK: $lapack_result
FFT: $with_fft
MPI: $mpi_result
FFLAGS=$FFLAGS
FCFLAGS=$FCFLAGS

Configuration complete; type \"make\" to start compiling])

