CC := mpicc
CFLAGS := -Wall -Werror -O2 -std=c99

prog1 = mmm_mpi.c matrix_checksum.c
prog2 = mandelbrot_mpi.c matrix_checksum.c 

exes = mmm_mpi mandelbrot_mpi 

all: $(exes)

mmm_mpi: $(prog1)
	$(CC) $(CFLAGS) $(prog1) -o $@ 

mandelbrot_mpi: $(prog2)
	$(CC) $(CFLAGS) $(prog2) -o $@ -lm

clean:
	rm -f $(exes)

