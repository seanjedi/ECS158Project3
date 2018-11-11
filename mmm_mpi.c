#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>

unsigned int matrix_checksum(int N, void *M, unsigned int size);

/////////////////////
//Multiply Function//
/////////////////////
void multiply_mpi(double* A, double* B, double* C, int N, int chunk, int myChunk, int rank){
	int start = chunk * rank;
	int stop = start + myChunk;

	for(int i = start ; i < stop; i++){
        	for(int k = 0; k < N; k++){
	                double temp = A[i*N + k];
                        for(int j = 0; j < N; j++){
          	              C[i*N + j] += temp * B[k*N + j];
                        }
                }
    }

	return;
}


int main(int argc, char **argv)
{
	//MPI Intialization//
	MPI_Init(&argc, &argv);
	int com_size, chunk, com_rank, chunk_size, last_chunk, my_chunk;
	MPI_Comm_size(MPI_COMM_WORLD, &com_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &com_rank);

	if(argc != 2){//if wrong number of inputs
		fprintf(stderr, "Usage: %s N\n", *argv);
		exit(1);
	}

	int N;

	if(!(N = atoi(argv[1]))){
		printf("Error: wrong matrix order (0 < N <= 2000)\n");
		exit(1);
	}
	if(N <= 0 || N > 2000){
                printf("Error: wrong matrix order (0 < N <= 2000)\n");
                exit(1);
	}
	
		chunk_size = N/(com_size);
		chunk = chunk_size;
		last_chunk = N - (chunk_size * (com_size - 1 ));

		if(com_rank == com_size - 1)
			my_chunk = last_chunk;
		else
			my_chunk = chunk_size;

	// double* A = malloc(sizeof(double) * chunk * chunk);
	// double* B = malloc(sizeof(double) * chunk * chunk);
	// double* C = malloc(sizeof(double) * chunk * chunk);
	// for(int i = 0; i < chunk * chunk; i++){ // Matrix initialization
	// 	A[i] = i/(chunk * com_rank) + i % (chunk * com_rank);
	// 	B[i] = i/(chunk * com_rank) + (i % (chunk * com_rank)) *2;
	// }
	double* A = malloc(sizeof(double) * N * N);
	double* B = malloc(sizeof(double) * N * N);
	double* C = malloc(sizeof(double) * N * N);
	for(int i = 0; i < N * N; i++){ // Matrix initialization
		A[i] = i/N + i%N;
		B[i] = i/N + (i%N)*2;
	}

	struct timespec before, after;

	if(com_rank == 0){
		// double* C = malloc(sizeof(double) * N * N);
		clock_gettime(CLOCK_MONOTONIC, &before);
	}
	
	//Implement MPI code//
	if(N < com_size){
		if(com_rank == 0)
			multiply_mpi(A, B, C, N, N, N, 0);
	}else{
		multiply_mpi(A, B, C, N, chunk, my_chunk, com_rank);

		//Receive other data
		if(com_rank == 0){
			for(int i = 1; i < com_size; i++){
				if(i  == com_size - 1)
					MPI_Recv(C, last_chunk, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
				else
					MPI_Recv(C, chunk, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
				
			}
		}else{
			//Send other data
			MPI_Send(C, chunk, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}

	//End MPI Code 
	if(com_rank == 0){
		clock_gettime(CLOCK_MONOTONIC, &after);

		unsigned long elapsed_ns = (after.tv_sec - before.tv_sec)*(1E9) + after.tv_nsec - before.tv_nsec;
		double seconds = elapsed_ns / (1E9);

		printf("Running time: %f secs\n", seconds);

		printf("A: %u\n", matrix_checksum(N, A, sizeof(double)));
		printf("B: %u\n", matrix_checksum(N, B, sizeof(double)));
		printf("C: %u\n", matrix_checksum(N, C, sizeof(double)));
	}

	MPI_Finalize();
    free(A);
    free(B);
    free(C);	
	
	return 0;
}
