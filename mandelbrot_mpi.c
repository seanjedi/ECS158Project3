#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include<mpi.h>

unsigned int matrix_checksum(int N, double *M);

void input_validation(char **argv);
int validate_int_boundary(char *arg, int start, int end, char *output);
int validate_double_boundary(char *arg, double start, double end, char *output);

void com0(int com_rank, int com_size, int argc, char **argv);
void coms(int com_rank, int com_size);
int compute_mandelbrot(char **argv);

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
        validate_int_boundary(argv[4], 50, 1000, "cutoff")))
        exit(1);
}

int compute_mandelbrot(char **argv) {

    printf("Hello\n");
    return 0;
}

void com0(int com_rank, int com_size, int argc, char **argv) {
    if (argc == 5) {
        input_validation(argv);
    } else {
        fprintf(stderr, 
        "Usage: ./mandelbrot_mpi_reference xcenter ycenter zoom cutoff\n");
        exit(1);
    }
    
    compute_mandelbrot(argv);
      
}

void coms(int com_rank, int com_size) {


}