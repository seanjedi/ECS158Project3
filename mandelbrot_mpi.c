#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>

unsigned int matrix_checksum(int N, double *M);

void input_validation(char **argv);
int validate_int_boundary(char *arg, int start, int end, char *output);
int validate_double_boundary(char *arg, double start, double end, char *output);

void com0(int com_rank, int com_size, int argc, char **argv);
void coms(int com_rank, int com_size);
int compute_mandelbrot(double xcenter, double ycenter, int my_chunk, int cutoff, double increment, int* matrix, int start);

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int com_size, com_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &com_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &com_rank);

    if(com_rank == 0)
        com0(com_rank, com_size, argc, argv);
    else
        coms(com_rank, com_size);
        
    MPI_Finalize();
}

/*
* Boundary check:
* Validate integer and boundary
*/
int validate_int_boundary(char *arg, int start, int end, char *output) {
    int input_length_1 = strlen(arg);
    for(int i = 0; i < input_length_1; i++)
        if(!isdigit(arg[i])){
            fprintf(stderr, "Error: wrong %s (%d <= N <= %d)\n", output, start, end);
            return 0;
        }

    int size = 0;
    size = atoi(arg);
    if (size < start || size > end) {
        fprintf(stderr, "Error: wrong %s (%d <= N <= %d)\n", output, start, end);
        return 0;
    }
    return 1;
}

/*
* Boundary check:
* Validate type double and boundary
*/
int validate_double_boundary(char *arg, double start, double end, char *output) {
    int input_length_1 = strlen(arg);
    int dot_flag = 0;

    for(int i = 0; i < input_length_1; i++) {
        if (arg[i] == '.') {
            dot_flag++;
        }

        if ((!isdigit(arg[i]) && (arg[i] != '.')) | (dot_flag > 1)){
            fprintf(stderr, "Error: wrong %s (%.6f  <= N <= %.6f)\n", output, start, end);
            return 0;
        }

    }

    double size = 0;
    size = atof(arg);
    if (size < start || size > end) {
        fprintf(stderr, "Error: wrong %s (%.6f <= N <= %.6f)\n", output, start, end);
        return 0;
    }
    return 1;
}

/*
* Call correct validation function per user argument
*/
void input_validation(char **argv) {
    if (!(validate_double_boundary(argv[1], -10.000000, 10.000000, "x-center") && 
        validate_double_boundary(argv[2], -10.000000, 10.000000, "y-center") && 
        validate_int_boundary(argv[3], 0, 100, "zoom") && 
        validate_int_boundary(argv[4], 0, 1000, "cutoff")))
        exit(1);
}

int compute_mandelbrot(double xcenter, double ycenter, int my_chunk, int cutoff, double increment, int* matrix, int start) {
    double zcurr = 0, znext = 0;
    for(int y = start; y < my_chunk; y++){
        for(int x = 0; x < 1024; x++){
            int xvalue = x * increment;
            int yvalue = y * increment;
            double complex C = xvalue + yvalue * I;
            for(int iteration = 0; iteration < cutoff; iteration++){
                znext = pow(zcurr, zcurr) + cabs(C);
                matrix[((y - start)*1024) + x] += 1;
                if(abs(znext - zcurr) > 2)
                    break;
            }
        }
    }
    return 0;
}

void com0(int com_rank, int com_size, int argc, char **argv) {
    int zoom, cutoff, chunk_size, last_chunk, start;
    double xcenter, ycenter, increment; 
    if (argc == 5) {
        input_validation(argv);
    } else {
        fprintf(stderr, 
        "Usage: ./mandelbrot_mpi_reference xcenter ycenter zoom cutoff\n");
        exit(1);
    }

    chunk_size = 1024/(com_size);
    last_chunk = 1024 - (chunk_size * (com_size - 1 ));
    xcenter = atof(argv[1]);
    ycenter = atof(argv[2]);
    zoom = atoi(argv[3]);
    cutoff = atoi(argv[4]);
    increment = pow(2,(double) -zoom);

    for(int i = 1; i < com_size; i ++){
        if(i == com_size - 1)
            MPI_Send(&last_chunk, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        else
            MPI_Send(&chunk_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    for(int i = 1; i < com_size; i ++){
        start = i * chunk_size;
        MPI_Send(&start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    MPI_Bcast(&xcenter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ycenter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cutoff, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&increment, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int *matrix = (int*) calloc(1024*1024, sizeof(int));

    compute_mandelbrot(xcenter, ycenter, chunk_size, cutoff, increment, matrix, 0);
    
    int offset = com_rank+chunk_size*1024;
    // int rowNumber = 0;

    for(int i = 1; i < com_size; i ++){
        if(i == com_size - 1) {
            int* buffer = calloc(1024*last_chunk, sizeof(double));
            MPI_Recv(buffer, 1024*last_chunk, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);

            int count = 0;
            for(int i = offset; i < offset+1024*last_chunk; i++) {
                matrix[i] += buffer[count];
                count++;
            }

            free(buffer);
        } else {
            int* buffer = calloc(1024*chunk_size, sizeof(double));
            MPI_Recv(buffer, 1024*chunk_size, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);

            int count = 0;
            for(int i = offset; i < offset+1024*chunk_size; i++) {
                matrix[i] += buffer[count];
                count++;
            }

            free(buffer);
        }

    }



    free(matrix);

      
} 

void coms(int com_rank, int com_size) {
    int cutoff, my_chunk, start;
    double xcenter, ycenter, increment;

    MPI_Recv(&my_chunk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
    MPI_Recv(&start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
    MPI_Bcast(&xcenter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ycenter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cutoff, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&increment, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int *matrix = (int*) calloc(1024*my_chunk, sizeof(int));

    compute_mandelbrot(xcenter, ycenter, my_chunk, cutoff, increment, matrix, start);

    MPI_Send(matrix, my_chunk * 1024, MPI_INT, 0, 0, MPI_COMM_WORLD);
    free(matrix);
}