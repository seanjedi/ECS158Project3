#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>
#include<string.h>

unsigned int matrix_checksum(int N, void *M, unsigned int size);
/////////////////////
//Multiply Function//
/////////////////////
void multiply_mpi(double* A, double* B, double* C, int N, int myChunk){
	// int start = chunk_size * rank;
	// int stop = start + myChunk;
	for(int i = 0 ; i < myChunk; i++){
        	for(int k = 0; k < N; k++){
	                double temp = A[i*N + k];
                        for(int j = 0; j < N; j++){
          	              C[i*N + j] += temp * B[k*N + j];
                        }
                }
    }
	return;
}

void com0(int com_rank, int com_size, int argc, char **argv){
	if(argc != 2){//if wrong number of inputs
		fprintf(stderr, "Usage: %s N\n", *argv);
		// MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}

	int N, chunk_size;

	if(!(N = atoi(argv[1]))){
		fprintf(stderr, "Error: wrong matrix order (N > 0)\n");
		exit(1);
		// MPI_Abort(MPI_COMM_WORLD, 1);

	}
	if(N <= 0){
		fprintf(stderr, "Error: wrong matrix order (N > 0)\n");
		exit(1);
		// MPI_Abort(MPI_COMM_WORLD, 1);
	}
	struct timespec before, after;
	clock_gettime(CLOCK_MONOTONIC, &before);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	chunk_size = N/(com_size);
	double* A = malloc(sizeof(double) * N * N);
	double* B = malloc(sizeof(double) * N * N);
	double* C = malloc(sizeof(double) * N * N);
	for(int i = 0; i < N * N; i++){ // Matrix initialization
		A[i] = i/N + i%N;
		B[i] = i/N + (i%N)*2;
		C[i] = 0;
	}

	if(N < com_size)
		multiply_mpi(A, B, C, N, N);
	else{
		int last_chunk = N - (chunk_size * (com_size - 1 ));
		for(int i = 1; i < com_size; i++){
			if(i == com_size - 1){
				MPI_Send((A+(N*i*chunk_size)), N*last_chunk, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send((B+(N*i*chunk_size)), N*last_chunk, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send((C+(N*i*chunk_size)), N*last_chunk, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}else{
				MPI_Send((A+(N*i*chunk_size)), N*chunk_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send((B+(N*i*chunk_size)), N*chunk_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send((C+(N*i*chunk_size)), N*chunk_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		
		multiply_mpi(A, B, C, N, chunk_size);

		// indexing for our C
		int offset = com_rank+chunk_size*N;
		
		for(int i = 1; i < com_size; i++){
			if(i  == com_size - 1) {
				double *buffer = (double*) calloc(N*last_chunk, sizeof(double));

				// using buffer instead of temp because MPI_Recv kept on receiving data from index 0 
				// regardless of how I tried to index the location, it would always receive the data at 0
				// I may be wrong about this, but hey its working.
				MPI_Recv(buffer, N*last_chunk, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
					
				int count = 0;
				for(int i = offset; i < offset+N*last_chunk; i++) {
					C[i] += buffer[count];
					count++;
				}
				offset += N*last_chunk;
				free(buffer);
			} else {
				double *buffer = (double*) calloc(N*chunk_size, sizeof(double));
					
				MPI_Recv(buffer, N*chunk_size, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);

				int count = 0;
				for(int i = offset; i < offset+N*chunk_size; i++) {
					C[i] += buffer[count];
					count++;
				}
				offset += N*chunk_size;
				free(buffer);
			}
		}

		clock_gettime(CLOCK_MONOTONIC, &after);
		unsigned long elapsed_ns = (after.tv_sec - before.tv_sec)*(1E9) + after.tv_nsec - before.tv_nsec;
		double seconds = elapsed_ns / (1E9);
		printf("Running time: %f secs\n", seconds);

		printf("A: %u\n", matrix_checksum(N, A, sizeof(double)));
		printf("B: %u\n", matrix_checksum(N, B, sizeof(double)));
		printf("C: %u\n", matrix_checksum(N, C, sizeof(double)));
		free(A);
    	free(B);
    	free(C);


	}
}

void coms(int com_rank, int com_size){
	int N, my_chunk, chunk_size;
	
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(N > com_size){

		chunk_size = N/(com_size);
		if(com_rank == com_size - 1)
			my_chunk = N - (chunk_size * (com_size - 1 ));
		else
			my_chunk = N/(com_size);
			
		double* A = (double*) calloc(N*my_chunk, sizeof(double));
		double* B = (double*) calloc(N*my_chunk, sizeof(double));
		double* C = (double*) calloc(N*my_chunk, sizeof(double));
		MPI_Recv(A, N*my_chunk, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
		MPI_Recv(B, N*my_chunk, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
		MPI_Recv(C, N*my_chunk, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);

		multiply_mpi(A, B, C, N, my_chunk);
		MPI_Send(C, N*my_chunk, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		free(A);
    	free(B);
    	free(C);
	}
}


int main(int argc, char **argv)
{
	//MPI Intialization//
	MPI_Init(&argc, &argv);
	int com_size, com_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &com_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &com_rank);

	if(com_rank == 0)
        com0(com_rank, com_size, argc, argv);
	else
        coms(com_rank, com_size);
		
	MPI_Finalize();
	
	return 0;
}