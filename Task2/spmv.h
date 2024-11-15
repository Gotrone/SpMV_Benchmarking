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

// COO functions
void dense_to_coo(double *mat, unsigned int n, unsigned int *nnz, int **row_indices, int **col_indices, double **values);
void coo_spmv(unsigned int n, unsigned int nnz, int *row_indices, int *col_indices, double *values, double *vec, double *result);

// CSC functions
void dense_to_csc(double *mat, unsigned int n, unsigned int *nnz, int **col_ptr, int **row_indices, double **values);
void csc_spmv(unsigned int n, int *col_ptr, int *row_indices, double *values, double *vec, double *result);

#endif
