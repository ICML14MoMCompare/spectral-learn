noinst_LTLIBRARIES = libgslhistogram.la 

pkginclude_HEADERS = gsl_histogram.h gsl_histogram2d.h

INCLUDES= -I$(top_builddir)


check_PROGRAMS = test test2d test_resample test2d_resample test_trap test2d_trap

TESTS = test test2d test_resample test2d_resample test_trap test2d_trap
EXTRA_DIST = urand.c

noinst_HEADERS = urand.c find.c find2d.c

test_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la
test2d_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la
test_resample_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la
test2d_resample_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la
test_trap_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la
test2d_trap_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la

test_SOURCES = test.c
test2d_SOURCES = test2d.c
test_resample_SOURCES = test_resample.c
test2d_resample_SOURCES = test2d_resample.c
test_trap_SOURCES = test_trap.c
test2d_trap_SOURCES = test2d_trap.c

CLEANFILES = test.txt test.dat

libgslhistogram_la_SOURCES = add.c  get.c init.c params.c reset.c file.c pdf.c gsl_histogram.h add2d.c get2d.c init2d.c params2d.c reset2d.c file2d.c pdf2d.c gsl_histogram2d.h calloc_range.c calloc_range2d.c copy.c copy2d.c maxval.c maxval2d.c oper.c oper2d.c stat.c


