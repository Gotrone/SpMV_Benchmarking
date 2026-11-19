#ifndef SPMV_H
#define SPMV_H

#include <stddef.h>        // Include stddef.h for size_t
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

// Dense matrix-vector multiplication
int my_dense(const unsigned int n, const double mat[], double vec[], double result[]);

// Convert dense matrix to GSL sparse matrix format (CSR)
gsl_spmatrix* convert_to_gsl_sparse(double *mat, unsigned int n);

// Free memory for GSL sparse matrix
void free_gsl_sparse(gsl_spmatrix *sparse);

#endif
