CPPFLAGS = -O2 -Iinclude
CXX = g++

all: objects main

#inc = includes/projections.h includes/admm.h includes/test.h includes/fista.h includes/frank_wolfe.h
obj = util.o tensor_decomp.o test.o tensor_whiten.o learn_Otilde.o projections.o tensor_symmetrize.o
objects: $(obj)
main: $(obj)

clean:
	-rm *.o main

