# Process this file with autoconf to produce a configure script.
AC_INIT(hdtg, 0.1, msuchard@ucla.edu)

# find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

# platform check
RHOST=`"${R_HOME}/bin/Rscript" -e 'cat(Sys.info()[["sysname"]])'`

# get RcppParallel directory
RCPPPAR_LIB=`"${R_HOME}/bin/Rscript" -e "RcppParallel::RcppParallelLibs()"`

SSE_FLAGS=`"${R_HOME}/bin/Rscript" -e 'if(RcppXsimd::supportsSSE()) cat(RcppXsimd::getSSEFlags(),"-DUSE_SSE")'`

AVX_FLAGS=`"${R_HOME}/bin/Rscript" -e 'if(RcppXsimd::supportsAVX()) cat(RcppXsimd::getAVXFlags(),"-DUSE_AVX","-mfma")'`

SIMD_FLAGS="-DUSE_SIMD"

AC_SUBST(SIMD_FLAGS)
AC_SUBST(RCPPPAR_LIB)
AC_SUBST(SSE_FLAGS)
AC_SUBST(AVX_FLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
