all: mmm_mpi.c matrix_checksum.c
	mpicc -Wall -Werror -O2 matrix_checksum.c mmm_mpi.c -o mmm_mpi

clean: 
	rm -f mmm_mpi
