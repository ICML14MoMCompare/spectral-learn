noinst_LTLIBRARIES = libgslrng.la 

pkginclude_HEADERS = gsl_rng.h

INCLUDES= -I$(top_builddir)

libgslrng_la_SOURCES = rng.c types.c default.c cmrg.c gfsr4.c slatec.c minstd.c mrg.c mt.c r250.c ran0.c ran1.c ran2.c ran3.c rand.c random.c rand48.c randu.c ranf.c ranlux.c ranlxs.c ranlxd.c ranmar.c taus.c transputer.c tt.c uni.c uni32.c vax.c zuf.c 

test_SOURCES = test.c
test_LDADD = libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

TESTS = test
check_PROGRAMS = test benchmark rng-dump

benchmark_SOURCES = benchmark.c 
benchmark_LDADD = libgslrng.la ../err/libgslerr.la ../utils/libutils.la

rng_dump_SOURCES = rng-dump.c 
rng_dump_LDADD = libgslrng.la ../err/libgslerr.la ../utils/libutils.la
