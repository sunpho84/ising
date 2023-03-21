AC_INIT([ising],[1.0],[fr.sanfilippo@gmail.com])

#include all configure flags
CONFIGURE_FLAGS="$*"
AC_DEFINE_UNQUOTED([CONFIGURE_FLAGS],["$CONFIGURE_FLAGS"],[Flags passed to configure])

#silent automake
AM_SILENT_RULES([yes])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIGURE_TIME=`date`
AC_DEFINE_UNQUOTED([CONFIGURE_TIME], "$CONFIGURE_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_TIME])

#parameters definition
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.14 no-define foreign subdir-objects])
AC_CONFIG_HEADERS([src/config.hpp])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#checks for programs
AC_PROG_CC
AC_PROG_FC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG([C++])
#AX_CXX_COMPILE_STDCXX(17,noext,mandatory)

LIBRARY_RESULT="---"
SUMMARY_RESULT="---"

AX_ASSEMBLY_REPORT

#check openmp
AX_OPENMP

#AX_CXXFLAGS_WARN_ALL
AC_MSG_CHECKING(whether compiler understands -Wall)
OLD_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Wall"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no)
  CXXFLAGS="$OLD_CXXFLAGS")


AC_ARG_ENABLE(measurement-cache,AS_HELP_STRING([--enable-measurement-cache],[Enable cache for measurement]),enable_measurement_cache="${enableval}",enable_measurement_cache="no")
AC_MSG_RESULT([enabling measurement cache... $enable_measurement_cache])
AS_IF([test "x$enable_measurement_cache" = "xyes"],
      [AC_DEFINE([useMeasurementCache],[1],[Enable cached measurement])],
      [AC_DEFINE([useMeasurementCache],[0],[Disable cached measurement])])

AC_ARG_ENABLE(local-energy-change,AS_HELP_STRING([--enable-local-energy-change],[Enable local energy change measurement]),enable_local_energy_change="${enableval}",enable_local_energy_change="no")
AC_MSG_RESULT([enabling local energy change measurement... $enable_local_energy_change])
AS_IF([test "x$enable_local_energy_change" = "xyes"],
      [AC_DEFINE([useLocalEnergyChange],[1],[Enable local energy change measurement])],
      [AC_DEFINE([useLocalEnergyChange],[0],[Disable local energy change measurement])])


AC_ARG_ENABLE(lookup-table-accept-reject,AS_HELP_STRING([--enable-lookup-table-accept-reject],[Enable lookup table for accept/reject step]),enable_lookup_table_accept_reject="${enableval}",enable_lookup_table_accept_reject="no")
AC_MSG_RESULT([enabling lookup table for accept/reject step... $enable_lookup_table_accept_reject])
AS_IF([test "x$enable_lookup_table_accept_reject" = "xyes"],
      [AC_DEFINE([useLookupTableAcceptReject],[1],[Enable lookup table for accept/reject step])],
      [AC_DEFINE([useLookupTableAcceptReject],[0],[Disable lookup table for accept/reject step])])

#result of conf to screen
AC_MSG_RESULT([
--------------------------------------------------
Configuration summary for $PACKAGE_NAME ($VERSION)
--------------------------------------------------
CC                  : $CC
CPPFLAGS            : $CPPFLAGS
CXXFLAGS            : $CXXFLAGS
CXX                 : $CXX
---
LDFLAGS             : $LDFLAGS
LIBS                : $LIBS
---
Installation prefix : $prefix
$LIBRARY_RESULT
$SUMMARY_RESULT
--------------------------------------------------
Configuration OK
--------------------------------------------------
])

#files to configure
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
