CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3

default: all

v4_seq:
	$(CC) $(CFLAGS) -o v4_seq v4_seq.c mmio.c

v4_omp: 
	$(CC) $(CFLAGS) -o v4_omp v4_omp.c mmio.c -fopenmp

v4_threads:
	$(CC) $(CFLAGS) -o v4_threads v4_threads.c mmio.c -lpthread

v4_cilk:
	$(CILKCC) $(CFLAGS) -o v4_cilk v4_cilk.c -fcilkplus

.PHONY: clean

all: v4_seq v4_omp v4_cilk v4_threads

clean:
	rm -f v4_seq v4_omp v4_cilk v4_threads
