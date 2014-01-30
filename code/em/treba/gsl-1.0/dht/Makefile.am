noinst_LTLIBRARIES = libgsldht.la 

pkginclude_HEADERS = gsl_dht.h

INCLUDES= -I$(top_builddir)

TESTS = test

check_PROGRAMS = test

test_LDADD = libgsldht.la ../specfunc/libgslspecfunc.la ../sys/libgslsys.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

test_SOURCES = test.c

libgsldht_la_SOURCES = dht.c
