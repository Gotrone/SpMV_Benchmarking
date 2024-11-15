#include "spmv.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>        // Needed for size_t
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

// Convert dense matrix to GSL sparse matrix in CSR format
gsl_spmatrix* convert_to_gsl_sparse(double *mat, unsigned int n)
{
    gsl_spmatrix *sparse = gsl_spmatrix_alloc(n, n);  // Allocating sparse matrix with size n x n

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            double value = mat[i * n + j];
            if (value != 0.0) {
                gsl_spmatrix_set(sparse, i, j, value);  // Insert non-zero values
            }
        }
    }
    // Removed the non-existent gsl_spmatrix_compact() function, as it's not needed
    return sparse;
}

// Free memory for CSR matrix
void free_gsl_sparse(gsl_spmatrix *sparse) {
    gsl_spmatrix_free(sparse);
}

// Function to convert dense matrix to COO format
void dense_to_coo(double *mat, unsigned int n, unsigned int *nnz, int **row_indices, int **col_indices, double **values) {
    unsigned int count = 0;

    // First pass to count non-zero elements
    for (unsigned int i = 0; i < n * n; i++) {
        if (mat[i] != 0.0) {
            count++;
        }
    }
    *nnz = count;

    // Allocate memory for COO arrays
    *row_indices = (int *)malloc(count * sizeof(int));
    *col_indices = (int *)malloc(count * sizeof(int));
    *values = (double *)malloc(count * sizeof(double));

    // Second pass to fill arrays
    unsigned int idx = 0;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            double val = mat[i * n + j];
            if (val != 0.0) {
                (*row_indices)[idx] = i;
                (*col_indices)[idx] = j;
                (*values)[idx] = val;
                idx++;
            }
        }
    }
}

// COO SpMV function
void coo_spmv(unsigned int n, unsigned int nnz, int *row_indices, int *col_indices, double *values, double *vec, double *result) {
    // Initialize result vector
    for (unsigned int i = 0; i < n; i++) {
        result[i] = 0.0;
    }

    // Perform multiplication
    for (unsigned int k = 0; k < nnz; k++) {
        result[row_indices[k]] += values[k] * vec[col_indices[k]];
    }
}

// Function to convert dense matrix to CSC format
void dense_to_csc(double *mat, unsigned int n, unsigned int *nnz, int **col_ptr, int **row_indices, double **values) {
    unsigned int count = 0;

    // First pass to count non-zero elements
    for (unsigned int i = 0; i < n * n; i++) {
        if (mat[i] != 0.0) {
            count++;
        }
    }
    *nnz = count;

    // Allocate arrays
    *col_ptr = (int *)malloc((n + 1) * sizeof(int));
    *row_indices = (int *)malloc(count * sizeof(int));
    *values = (double *)malloc(count * sizeof(double));

    // Initialize col_ptr
    for (unsigned int i = 0; i <= n; i++) {
        (*col_ptr)[i] = 0;
    }

    // Count non-zero elements per column
    for (unsigned int j = 0; j < n; j++) {
        for (unsigned int i = 0; i < n; i++) {
            if (mat[i * n + j] != 0.0) {
                (*col_ptr)[j + 1]++;
            }
        }
    }

    // Cumulative sum to get col_ptr
    for (unsigned int j = 0; j < n; j++) {
        (*col_ptr)[j + 1] += (*col_ptr)[j];
    }

    // Fill row_indices and values
    unsigned int *col_counts = (unsigned int *)calloc(n, sizeof(unsigned int));
    for (unsigned int j = 0; j < n; j++) {
        for (unsigned int i = 0; i < n; i++) {
            double val = mat[i * n + j];
            if (val != 0.0) {
                unsigned int idx = (*col_ptr)[j] + col_counts[j];
                (*row_indices)[idx] = i;
                (*values)[idx] = val;
                col_counts[j]++;
            }
        }
    }
    free(col_counts);
}

// CSC SpMV function
void csc_spmv(unsigned int n, int *col_ptr, int *row_indices, double *values, double *vec, double *result) {
    // Initialize result vector
    for (unsigned int i = 0; i < n; i++) {
        result[i] = 0.0;
    }

    // Perform multiplication
    for (unsigned int j = 0; j < n; j++) {
        for (unsigned int idx = col_ptr[j]; idx < col_ptr[j + 1]; idx++) {
            result[row_indices[idx]] += values[idx] * vec[j];
        }
    }
}



