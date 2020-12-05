CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3

default: all

v3_seq:
	$(CC) $(CFLAGS) -o v3_seq v3_seq.c mmio.c

v3_omp: 
	$(CC) $(CFLAGS) -o v3_omp v3_omp.c mmio.c -fopenmp

v3_cilk:
	$(CILKCC) $(CFLAGS) -o v3_cilk v3_cilk.c mmio.c -fcilkplus

.PHONY: clean

all: v3_seq v3_omp v3_cilk


clean:
	rm -f v3_seq v3_omp v3_cilk
