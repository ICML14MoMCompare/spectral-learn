# -*-makefile-*-

noinst_LTLIBRARIES = libgslmultiroots.la 

pkginclude_HEADERS = gsl_multiroots.h

noinst_HEADERS = enorm.c dogleg.c

INCLUDES= -I$(top_builddir)

libgslmultiroots_la_SOURCES = fdjac.c fsolver.c fdfsolver.c convergence.c  newton.c gnewton.c dnewton.c broyden.c hybrid.c hybridj.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test_funcs.c test_funcs.h
test_LDADD = libgslmultiroots.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../sys/libgslsys.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

