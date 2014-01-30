CC = g++
CPPFLAGS = -O2 -fPIC 
LDFLAGS = -Wl,-rpath,/usr/lib

all: objects libADMM.so

inc = includes/projections.h includes/admm.h includes/admmpythonwrapper.h 
obj = projections.o accel_admm.o admmpythonwrapper.o util.o redsvd_util.o

objects: $(obj)

main: $(obj)

libADMM.so: $(obj)
	$(CC) -shared $(obj) -o libADMM.so

clean:
	rm -f *.o
	rm -f libADMM.so
