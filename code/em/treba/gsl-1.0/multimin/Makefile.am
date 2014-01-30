noinst_LTLIBRARIES = libgslmultimin.la 

pkginclude_HEADERS = gsl_multimin.h

noinst_HEADERS = 

INCLUDES= -I$(top_builddir)

libgslmultimin_la_SOURCES = fdfminimizer.c steepest_descent.c conjugate_fr.c conjugate_pr.c convergence.c diff.c vector_bfgs.c

noinst_HEADERS = directional_minimize.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultimin.la ../min/libgslmin.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../linalg/libgsllinalg.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

