#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>        // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>     // GSL sparse matrix library
#include <gsl/gsl_spblas.h>       // GSL sparse BLAS functions
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 16384
#define DEFAULT_DENSITY 0.1

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
    unsigned int nnz = 0;

    srand(seed);

    for (unsigned int i = 0; i < n * n; i++) {
        if ((rand() % 100) / 100.0 < density) {
            mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
            nnz++;
        } else {
            mat[i] = 0;
        }
    }

    return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
    srand(seed);

    for (unsigned int i = 0; i < size; i++) {
        vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }

    return size;
}

int is_nearly_equal(double x, double y)
{
    const double epsilon = 1e-5;
    return fabs(x - y) <= epsilon * fabs(x);
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
    for(unsigned int i = 0; i < size; i++) {
        if (!is_nearly_equal(ref[i], result[i]))
            return 0;
    }
    return 1;
}

int main(int argc, char *argv[])
{
    int size;        // Number of rows and columns (size x size matrix)
    double density;  // Ratio of non-zero values

    // Parsing command-line arguments
    if (argc < 2) {
        size = DEFAULT_SIZE;
        density = DEFAULT_DENSITY;
    } else if (argc < 3) {
        size = atoi(argv[1]);
        density = DEFAULT_DENSITY;
    } else {
        size = atoi(argv[1]);
        density = atof(argv[2]);
    }

    // Declare timing variables
    timeinfo start, now;

    // Allocate memory for matrices and vectors
    double *mat = (double *) malloc(size * size * sizeof(double));
    double *vec = (double *) malloc(size * sizeof(double));
    double *refsol = (double *) malloc(size * sizeof(double));
    double *mysol = (double *) malloc(size * sizeof(double));

    // Check for successful memory allocation
    if (!mat || !vec || !refsol || !mysol) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix and vector
    unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
    populate_vector(vec, size, 2);

    printf("Matrix size: %d x %d (%d elements)\n", size, size, size * size);
    printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size * size) * 100.0);

    //
    // Dense computation using CBLAS (e.g., GSL's CBLAS implementation)
    //
    printf("Dense computation\n----------------\n");

    timestamp(&start);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);
    timestamp(&now);

    printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

    //
    // Using your own dense implementation
    //
    timestamp(&start);
    my_dense(size, mat, vec, mysol);
    timestamp(&now);

    printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

    if (check_result(refsol, mysol, size))
        printf("Result is ok!\n");
    else
        printf("Result is wrong!\n");

    //
    // Sparse computation using GSL
    //
    printf("\nSparse computation using GSL\n----------------\n");

    // Convert the dense matrix to GSL sparse format (CSR)
    gsl_spmatrix *gsl_sparse = convert_to_gsl_sparse(mat, size);
    gsl_vector_view gsl_vec = gsl_vector_view_array(vec, size);
    gsl_vector *gsl_result = gsl_vector_alloc(size);

    timestamp(&start);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, gsl_sparse, &gsl_vec.vector, 0.0, gsl_result);
    timestamp(&now);

    printf("Time taken by GSL sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));

    // Check the GSL result
    for (unsigned int i = 0; i < size; i++) {
        mysol[i] = gsl_vector_get(gsl_result, i);
    }

    if (check_result(refsol, mysol, size))
        printf("GSL sparse result is ok!\n");
    else
        printf("GSL sparse result is wrong!\n");

    //
    // COO SpMV Computation
    //
    printf("\nCOO SpMV Computation\n-------------------\n");

    // Declare variables for COO
    unsigned int coo_nnz;
    int *coo_row_indices, *coo_col_indices;
    double *coo_values;

    // Convert to COO format
    dense_to_coo(mat, size, &coo_nnz, &coo_row_indices, &coo_col_indices, &coo_values);

    timestamp(&start);
    coo_spmv(size, coo_nnz, coo_row_indices, coo_col_indices, coo_values, vec, mysol);
    timestamp(&now);

    printf("Time taken by COO SpMV: %ld ms\n", diff_milli(&start, &now));

    if (check_result(refsol, mysol, size))
        printf("COO SpMV result is ok!\n");
    else
        printf("COO SpMV result is wrong!\n");

    //
    // CSC SpMV Computation
    //
    printf("\nCSC SpMV Computation\n-------------------\n");

    // Declare variables for CSC
    unsigned int csc_nnz;
    int *csc_col_ptr;
    int *csc_row_indices;
    double *csc_values;

    // Convert to CSC format
    dense_to_csc(mat, size, &csc_nnz, &csc_col_ptr, &csc_row_indices, &csc_values);

    timestamp(&start);
    csc_spmv(size, csc_col_ptr, csc_row_indices, csc_values, vec, mysol);
    timestamp(&now);

    printf("Time taken by CSC SpMV: %ld ms\n", diff_milli(&start, &now));

    if (check_result(refsol, mysol, size))
        printf("CSC SpMV result is ok!\n");
    else
        printf("CSC SpMV result is wrong!\n");

    //
    // Free allocated resources
    //

    // Free COO resources
    free(coo_row_indices);
    free(coo_col_indices);
    free(coo_values);

    // Free CSC resources
    free(csc_col_ptr);
    free(csc_row_indices);
    free(csc_values);

    // Free GSL resources
    gsl_spmatrix_free(gsl_sparse);
    gsl_vector_free(gsl_result);

    // Free other resources
    free(mat);
    free(vec);
    free(refsol);
    free(mysol);

    return 0;
}
