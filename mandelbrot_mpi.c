#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>

unsigned int matrix_checksum(int N, void *M, unsigned int size);

void input_validation(char **argv);
int validate_int_boundary(char *arg, int start, int end, char *output);
int validate_double_boundary(char *arg, double start, double end, char *output);
void print_matrix(int* M, int matrix_size, int n);

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

        if ((!isdigit(arg[i]) && (arg[i] != '.') && (arg[i] != '-')) | (dot_flag > 1)){
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
    
    double complex zcurr = 0, znext = 0;
    for(int y = start; y < start + my_chunk; y++){
        for(int x = 0; x < 1024; x++){
            double xvalue = xcenter + (512 * increment) - x * increment;
            double yvalue = ycenter - (512 * increment) + y * increment ;  // Need to rotate it 90
         
            double complex C = xvalue + yvalue * I;
            for(int iteration = 0; iteration < cutoff; iteration++){
                znext = zcurr * zcurr + C;
             
                matrix[((y - start)*1024) + x] += 1;
                if(abs(znext) > 2)
                    break;
                
                zcurr = znext;
            }
            zcurr = 0;
            znext = 0;
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

    struct timespec before, after;
    clock_gettime(CLOCK_MONOTONIC, &before);

    chunk_size = 1024/(com_size);
    last_chunk = 1024 - (chunk_size * (com_size - 1 ));
    xcenter = atof(argv[1]);
    ycenter = atof(argv[2]);
    zoom = atoi(argv[3]);
    cutoff = atoi(argv[4]);
    increment = pow(2,-zoom);
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
    
    int offset = chunk_size*1024;

    for(int i = 1; i < com_size; i ++){
        if(i == com_size - 1) {

            int* buffer = (int*) calloc(1024*last_chunk, sizeof(int));
            MPI_Recv(buffer, 1024*last_chunk, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);

            int count = 0;
            for(int i = offset; i < offset+1024*last_chunk; i++) {
                matrix[i] += buffer[count];
                count++;
            }
            offset += 1024*last_chunk;
            free(buffer);
        } else {
            int* buffer = (int*) calloc(1024*chunk_size, sizeof(int));
            MPI_Recv(buffer, 1024*chunk_size, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, 0);
            int count = 0;
            for(int i = offset; i < offset+1024*chunk_size; i++) {
                matrix[i] += buffer[count];
                // printf("buffer[%d] = %d\n", count, buffer[count]);
                count++;
            }
            offset += 1024*chunk_size;
            free(buffer);
        }

    }


    clock_gettime(CLOCK_MONOTONIC, &after);
    unsigned long elapsed_ns = (after.tv_sec - before.tv_sec)*(1E9) + after.tv_nsec - before.tv_nsec;
    double seconds = elapsed_ns / (1E9);
    // printf("==========matrix============\n");
    // print_matrix(matrix, 1024*1024, 1024);     
    
    printf("Running time: %f secs\n", seconds);
    printf("M: %u\n", matrix_checksum(1024, matrix, sizeof(int)));


    char name[255];
    sprintf(name, "mandel_%.6f_%.6f_%d_%d.pgm", xcenter, ycenter, zoom, cutoff);
    
    FILE *fp;
    fp = fopen(name, "w+");
    fprintf(fp, "P2\n");
    fprintf(fp, "1024 1024\n");
    fprintf(fp, "%d\n", cutoff);
    int count = 0;
    for(int i = 0; i < 1024*1024; i++) {
        fprintf(fp, "%d ", matrix[i]);
        count++;
        if(count == 1024) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fclose(fp);

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

    // print_matrix(matrix, 1024*my_chunk, 1024);

    compute_mandelbrot(xcenter, ycenter, my_chunk, cutoff, increment, matrix, start);
    // print_matrix(matrix, 1024*my_chunk, 1024);
    MPI_Send(matrix, my_chunk * 1024, MPI_INT, 0, 0, MPI_COMM_WORLD);
    free(matrix);
}

void print_matrix(int* M, int matrix_size, int n) {
    printf("-----Printing Matrix-----\n");
    int count = 0;
    for(int i = 0; i < matrix_size; i++) {

        printf("%d , ", M[i]);
        count++;
        if(count%n == 0) {
            printf("\n");
            count = 0;
        }

        if(i > 1024) {
            break;
        }

    }
    printf("\n");
}