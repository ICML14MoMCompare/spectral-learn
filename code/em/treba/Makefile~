############################################
# Treba makefile                           #
# To compile with NVIDIA CUDA support run: #
# make CUDA=1                              #
############################################

PREFIX = /usr/local
CC = gcc -Wall
BINPREFIX = $(PREFIX)/bin/
MANPREFIX = $(PREFIX)/share/man/man1/
RM = /bin/rm -f
CFLAGS = -O3 -ffast-math

ifeq ($(CUDA),1)
	CUDA_INSTALL_PATH ?= /usr/local/cuda
	CFLAGS += -DUSE_CUDA
	LFLAGS = -lm -lpthread -L$(CUDA_INSTALL_PATH)/lib -lcudart -lgsl -lgslcblas
	TREBADEPS = treba.o dffa.o gibbs.o observations.o io.o treba.h treba_cuda.o fastlogexp.h
	TREBACMD = $(CC) $(CFLAGS) -DUSE_CUDA -o treba treba_cuda.o treba.o dffa.o gibbs.o observations.o io.o $(LFLAGS)
else
	LFLAGS = -lm -lpthread -lgsl -lgslcblas
	TREBADEPS = treba.o dffa.o gibbs.o observations.o io.o treba.h fastlogexp.h
	TREBACMD = $(CC) $(CFLAGS) -o treba treba.o dffa.o gibbs.o observations.o io.o $(LFLAGS)
endif


all:	treba

treba:	$(TREBADEPS)
	$(TREBACMD)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

treba_cuda.o: treba_cuda.cu treba.h fastlogexp.h
	nvcc -m64 -I$(CUDA_INSTALL_PATH)/include -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -o treba_cuda.o -c treba_cuda.cu

clean:
	$(RM) treba treba.o dffa.o gibbs.o observations.o io.o treba_cuda.o

install: treba treba.1
	-@if [ ! -d $(BINPREFIX) ]; then mkdir -p $(BINPREFIX); fi
	-@if [ ! -d $(MANPREFIX) ]; then mkdir -p $(MANPREFIX); fi
	cp treba $(BINPREFIX)
	cp ./man/treba.1 $(MANPREFIX)
