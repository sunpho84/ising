#include m4 dir
ACLOCAL_AMFLAGS=-I m4

#init to empty
include_HEADERS=
bin_PROGRAMS=
BUILT_SOURCES=
CLEANFILES=
assembly_reports=

#rule to create assembly output and get it filtered
%.s: %.cpp %.o $(top_srcdir)/config/asm_inspect.sh
	$(CXXCOMPILE) -S $< -c -o -|$(top_srcdir)/config/asm_inspect.sh $@

bin_PROGRAMS+= \
        $(top_builddir)/bin/ising

__top_builddir__bin_ising_SOURCES= \
	%D%/src/ising.cpp

#generate assembly report if asked
if ASSEMBLY_REPORT
BUILT_SOURCES+=$(assembly_reports)
CLEANFILES+=$(assembly_reports)
endif
