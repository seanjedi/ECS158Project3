all: mmm_mpi.c
	gcc -Wall -Werror -O2 -mpicc matrix_checksum.c mmm_mpi.c -o mmm_mpi

clean: 
	rm -f mmm_mpi
