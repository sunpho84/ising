# usage: AX_GET_ENABLE(switch_name,default,additional_info)
AC_DEFUN([AX_GET_ENABLE], [

AC_ARG_ENABLE($1,
	AS_HELP_STRING([--enable-$1],[Enable $1 $3]),
	enable_$1="${enableval}",
	enable_$1="$2")
AC_MSG_RESULT([enabling $1 ... ${enable_$1}])
SUMMARY_RESULT="$SUMMARY_RESULT
$1 enabled        : $enable_$1"

if test "$enable_$1" == "yes"
then
	if test "$$1_found" == "no"
	then
		AC_MSG_ERROR(["Cannot enable $1, library/functionality not found, please provide/improve hint using --with-$1 flag (if available)"])
	fi
fi
])
