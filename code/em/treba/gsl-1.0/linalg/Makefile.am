noinst_LTLIBRARIES = libgsllinalg.la 

pkginclude_HEADERS = gsl_linalg.h

INCLUDES= -I$(top_builddir)

libgsllinalg_la_SOURCES = multiply.c tridiag.c tridiag.h lu.c luc.c hh.c qr.c qrpt.c svd.c householder.c householdercomplex.c cholesky.c symmtd.c hermtd.c bidiag.c balance.c

noinst_HEADERS =  givens.c apply_givens.c svdstep.c tridiag.h 

TESTS = test

check_PROGRAMS = test

test_LDADD = libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../sys/libgslsys.la ../test/libgsltest.la ../utils/libutils.la

test_SOURCES = test.c


