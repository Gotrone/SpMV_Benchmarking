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
