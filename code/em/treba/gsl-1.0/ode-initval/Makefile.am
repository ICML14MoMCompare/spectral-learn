noinst_LTLIBRARIES = libgslodeiv.la 

pkginclude_HEADERS = gsl_odeiv.h

INCLUDES= -I$(top_builddir)

libgslodeiv_la_SOURCES = control.c cstd.c evolve.c step.c rk2.c rk2imp.c rk4.c rk4imp.c rkf45.c rk8pd.c rkck.c bsimp.c gear1.c gear2.c

noinst_HEADERS = odeiv_util.h

check_PROGRAMS = test

TESTS = test

test_LDADD = libgslodeiv.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../test/libgsltest.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../permutation/libgslpermutation.la ../vector/libgslvector.la ../block/libgslblock.la ../utils/libutils.la

test_SOURCES = test.c

