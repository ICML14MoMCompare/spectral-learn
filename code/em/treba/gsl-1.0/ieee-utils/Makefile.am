noinst_LTLIBRARIES = libgslieeeutils.la 

pkginclude_HEADERS = gsl_ieee_utils.h

libgslieeeutils_la_SOURCES =  print.c make_rep.c gsl_ieee_utils.h env.c fp.c read.c

noinst_HEADERS = fp-aix.c fp-darwin.c fp-hpux.c fp-hpux11.c fp-irix.c fp-m68klinux.c fp-ppclinux.c fp-solaris.c fp-sparclinux.c fp-sunos4.c fp-tru64.c fp-unknown.c fp-x86linux.c fp-freebsd.c fp-os2emx.c fp-netbsd.c fp-openbsd.c endian.c standardize.c

INCLUDES= -I$(top_builddir)

TESTS = test
check_PROGRAMS = test
test_LDADD = libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la
test_SOURCES = test.c
