#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>
#include<string.h>

unsigned int matrix_checksum(int N, void *M, unsigned int size);
void print_matrix(double* M, int matrix_size, int n);
/////////////////////
//Multiply Function//
/////////////////////
void multiply_mpi(double* A, double* B, double* C, int N, int chunk_size, int myChunk, int rank){
	int start = chunk_size * rank;
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
	int com_size, com_rank, chunk_size, last_chunk, my_chunk;
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
	if(N <= 0){
	    printf("Error: wrong matrix order (0 < N <= 2000)\n");
	    exit(1);
	}
	
		chunk_size = N/(com_size);
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
	double* temp = malloc(sizeof(double) * N * N);
	for(int i = 0; i < N * N; i++){ // Matrix initialization
		A[i] = i/N + i%N;
		B[i] = i/N + (i%N)*2;
		temp[i] = 0;
		C[i] = 0;
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
		
		multiply_mpi(A, B, temp, N, chunk_size, my_chunk, com_rank);	
		//Receive matrix data
		if(com_rank == 0){
			// printf("==========temp0============\n");
			// print_matrix(temp, N*N, N);
			for(int i = 0; i < N*N; i++)
				C[i] += temp[i];

			// indexing for our C
			int offset = com_rank+chunk_size*N;

			for(int i = 1; i < com_size; i++){
				if(i  == com_size - 1) {
					double *buffer = (double*) malloc(sizeof(double)*N*last_chunk);
					
					for(int i = 0; i < N*last_chunk; i++) {
						buffer[i] = 0.0;
					}

					// using buffer instead of temp because MPI_Recv kept on receiving data from index 0 
					// regardless of how I tried to index the location, it would always receive the data at 0
					// I may be wrong about this, but hey its working.
					MPI_Recv(&buffer[0], N*last_chunk, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
					
					int count = 0;
					for(int i = offset; i < offset+N*last_chunk; i++) {
						// printf("buffer[%d] = %lf\n", i, buffer[count]);
						C[i] += buffer[count];
						count++;
					}
					offset += N*last_chunk;

					free(buffer);
				} else {
					double *buffer = (double*) malloc(sizeof(double)*N*chunk_size);
					for(int i = 0; i < N*chunk_size; i++) {
						buffer[i] = 0.0;
					}
					
					MPI_Recv(&buffer[0], N*chunk_size, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
					
					int count = 0;
					for(int i = offset; i < offset+N*chunk_size; i++) {
						// printf("buffer[%d] = %lf\n", i, buffer[count]);
						C[i] += buffer[count];
						count++;
					}
					offset += N*chunk_size;

					free(buffer);
				}
			}
		}else{	
			//Send matrix data
			MPI_Send((temp+(N*com_rank*chunk_size)), my_chunk*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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

	free(A);
    free(B);
    free(C);	
	MPI_Finalize();
	
	return 0;
}

void print_matrix(double* M, int matrix_size, int n) {
    printf("-----Printing Matrix-----\n");
    int count = 0;
    for(int i = 0; i < matrix_size; i++) {

        printf("%f , ", M[i]);
        count++;
        if(count%n == 0) {
        	printf("\n");
        	count = 0;
        }

    }
    printf("\n");
}