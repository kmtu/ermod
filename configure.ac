#                                               -*- Autoconf -*-

AC_PREREQ([2.59])
#AC_INIT(enerepr, 0.1.%rev%, Matubayasi Laboratory)
AC_INIT(enerepr, 0.1.221, Matubayasi Laboratory)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_MACRO_DIR(m4)

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC([icc gcc cc fcc])
AC_PROG_CC_C_O
AC_PROG_F77([ifort gfortran xlf f77 frt g77 xlf77])
AC_PROG_FC([ifort gfortran xlf f95 f90 g95 xlf90])

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran 77)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

AS_IF([test "x$F77" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 77). Set environment variable F77 to work around")]
      [])
AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around")]
      [])

AC_OPENMP

AC_DEFUN([MY_CHECK_BLAS], [
  AC_LANG_PUSH(Fortran 77)
  AC_MSG_CHECKING([whether $1 is available])
  save_LIBS=$LIBS
  LIBS="$LIBS $2"
  AC_LINK_IFELSE([AC_LANG_CALL([], [sgemm])],
                 [blas_ok=yes], [blas_ok=no])
  LIBS=$save_LIBS
  AC_MSG_RESULT([$blas_ok])
  AC_LANG_POP
])

AC_ARG_WITH([mkl],
  AS_HELP_STRING([--with-mkl[[=(MKL version)]]],
    [Automatically set environmental variables for MKL; for a specific version, use this option like --with-mkl=10.2.3]))

AS_IF([test "x$with_mkl" != "x"],
  [case x`uname -m` in # (
  xx86_64) :
  mkl_plat=em64t
  mkl_platnew=intel64
  ;; # (
  xi386) :
  mkl_plat=32
  mkl_platnew=ia32
  ;; # (
  *) :
  ;;
  esac
  AS_IF([test "x$with_mkl" = "xyes"],
    AS_IF([test "x$MKLROOT" = "x"],
      [AC_MSG_FAILURE([could not get MKL version; set MKLROOT or use --with-mkl=(version number)])],
      mklroot="$MKLROOT"),
    mklroot="/opt/intel/mkl/$with_mkl")
  # Here comes super-duper crazily dirty MKL versioning checker
  # check 10.0 - 10.2
  with_blas="-L$mklroot/lib/$mkl_plat -lmkl_solver_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -lguide -lpthread"
  MY_CHECK_BLAS([MKL 10.0-10.2], 
                [$with_blas])
  AS_IF([test "x$blas_ok" = "xyes"],
        [],
	# check for 10.3
	[with_blas="-L$mklroot/lib/$mkl_platnew -lmkl_rt $OPENMP_FCFLAGS -lpthread"
         MY_CHECK_BLAS([MKL 10.3+],
                       [$with_blas])
         AS_IF([test "x$blas_ok" = "xyes"],
               [],
               [AC_MSG_WARN([Failed to find appropriate MKL. If you are using old mkl, please update])])])
  with_lapack="-L$mklroot/lib/$mkl_plat -lmkl_lapack"
  FCFLAGS="-I$mklroot/include $FCFLAGS"
])

AC_ARG_WITH([fftw],
  AS_HELP_STRING([--with-fftw[[=(fftw root, default to /usr)]]],
    [Set fftw root directory]))

case x$with_fftw in # (
  xyes) :
  with_fftw=/usr
  ;; # (
  xno) :
  with_fftw=
  ;; # (
  *) :
  # do nothing
  ;;
esac

# FIXME: AC_SUBST should be used
AS_IF([test "x$with_fftw" != "x"],
  [FFTWROOT=$with_fftw
   CFLAGS="$CFLAGS -I$FFTWROOT/include"
   FCFLAGS="$FCFLAGS -I$FFTWROOT/include"
   LDFLAGS="$LDFLAGS -L$FFTWROOT/lib"])

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
# following switches for: default, gfortran, ifort, xlf, Fujitsu f90
for flag in "" "-ffixed-line-length-132" "-extend-source" "-qfixed=132" "-w"; do
FFLAGS="$save_FFLAGS $flag"
  # some compiler cuts line in 72 column, others 80 column
  AC_COMPILE_IFELSE([      program main
      implicit none
      integer :: foobarbazbaaz
                                                                      foobarbazbaaz=1
      end program], [my_check_f77_extend=yes], [my_check_f77_extend=no])
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
# for default, gfortran, ifort, xlf, Fujitsu f90
for flag in "" "-fdefault-real-8" "-autodouble" "-qrealsize=8" "-A d"; do
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

AC_ARG_ENABLE([defaults],
	AS_HELP_STRING([--disable-defaults],
		[Disable deafult flags]))

AS_IF([test x$enable_defaults != xno],
CFLAGS="$CFLAGS -g"
FFLAGS="$FFLAGS -g -cpp -Dtrjctry"
FCFLAGS="$FCFLAGS -g -cpp -Dtrjctry")

AC_ARG_ENABLE([perf],
	AS_HELP_STRING([--enable-perf],
		[(internal use only)]),
	AS_IF([test x$enable_perf = xyes],
		     [FFLAGS="$FFLAGS -fno-omit-frame-pointer -DPERF"
	              FCFLAGS="$FCFLAGS -fno-omit-frame-pointer -DPERF"
                      CCFLAGS="$CCFLAGS -fno-omit-frame-pointer"
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

AM_CONDITIONAL([USE_FFTW], [test x$with_fft = xfftw])
AM_CONDITIONAL([USE_MKL], [test x$with_fft = xmkl])

AS_IF([test "x$lapack_exist" != xyes],
	    [AC_MSG_ERROR(
		[dnl
LAPACK library was not found.
Try specifying --with-lapack=<lib> options or set FCFLAGS/FFLAGS/FLIBS environment variables to fix the problem.])])

lapack_result="$LAPACK_LIBS $BLAS_LIBS"


AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_RESULT([dnl
LAPACK: $lapack_result
FFT: $with_fft
MPI: $mpi_result
FFLAGS=$FFLAGS
FCFLAGS=$FCFLAGS

Configuration complete; type \"make\" to start compiling])

