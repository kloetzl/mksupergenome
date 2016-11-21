SUBDIRS = libs .
DIST_SUBDIRS = $(SUBDIRS)

bin_PROGRAMS = mksupergenome
mksupergenome_SOURCES = mksupergenome.cxx io.cxx io.h sequence.cxx sequence.h esa.cxx esa.h process.cxx process.h global.h
mksupergenome_CPPFLAGS = -std=c++14 -Wall -Wextra -I$(top_srcdir)/libs
mksupergenome_CXXFLAGS = -O3 -ggdb
mksupergenome_LDADD = $(top_builddir)/libs/libpfasta.a $(top_builddir)/libs/libcompat.a

.PHONY: format debug hack

hack:
	$(MAKE) -C libs all

debug: CXXFLAGS=-O0 -fno-omit-frame-pointer -ggdb
debug: hack mksupergenome

format:
	clang-format -i *.h *.cxx
