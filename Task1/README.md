# SpMV: Sparse Matrix-Vector Product

This code is based on the use of GSL (GNU Scientific Library) for the
implementation of the baseline operations used for comparison:
- dense matrix-vector product: `cblas_dgemv()`, you need to link against *libgslcblas*
- sparse matrix-vector product: `gsl_spblas_dgemv()`, you need to link against *libgsl*

The dense product, cblas_dgemv(), can be found in other CBLAS
implementation. You just need to change the library to be linked,
eg. `-lopenblas` instead of `-lgslcblas`

The basetype in GSL for working with sparse matrices is `gsl_spmatrix`.
GSL also provides functions to help convert you dense matrices into a
sparse format.


## Features

- **Dense Matrix-Vector Product** using the CBLAS library (GSL or OpenBLAS).
- **Sparse Matrix-Vector Product** using GSL’s `gsl_spblas_dgemv()` for sparse matrices in CSR format.
- Custom implementations for dense and sparse matrix-vector multiplication.
- Performance comparison between dense and sparse matrix operations.

## Contents
- [Usage](#usage)
- [Performance](#performance)
- [Matrix Formats](#matrix-formats)
- [Project Structure](#project-structure)

## Usage

1. **Running the Program**:
   You can run the program with default settings (matrix size 1024x1024, 25% non-zero elements for sparse matrix):

   `./spmv`

    You can also specify custom matrix sizes and densities for the sparse matrix (e.g., matrix size 512 and 50% density):

   `./spmv 512 0.5`


2. **Sample Output**:
   ```bash
    Matrix size: 1024 x 1024 (1048576 elements)
    261996 non-zero elements (24.99%)

    Dense computation
    ----------------
    Time taken by CBLAS dense computation: 0 ms
    Time taken by my dense matrix-vector product: 3 ms
    Result is ok!

    Sparse computation using GSL
    ----------------
    Time taken by GSL sparse matrix-vector product: 1 ms
    GSL sparse result is ok!

## Performance

This project compares the performance of the following operations:

- **Dense Matrix-Vector Product**:
    - Using **CBLAS** (provided by GSL or OpenBLAS).
    - Custom dense matrix-vector multiplication implementation.
    
- **Sparse Matrix-Vector Product**:
    - Using GSL's `gsl_spblas_dgemv()` for sparse matrix-vector multiplication in **CSR (Compressed Sparse Row)** format.

The execution times of these operations are printed after each computation, allowing you to observe how the performance varies with matrix size and density.

## Matrix Formats
### Dense Matrix Format

A dense matrix is stored in row-major order, meaning that all elements of a row are stored contiguously in memory. This format is inefficient when most of the elements are zero, as it uses extra memory and computational resources.
### Sparse Matrix Format (CSR)

Sparse matrices are stored using the Compressed Sparse Row (CSR) format. This format stores only the non-zero elements and their indices, resulting in significant memory savings for large, sparse matrices.

In CSR format, the matrix is represented using three arrays:

    values[]: Contains all non-zero values.
    col_indices[]: Contains the column index for each non-zero value.
    row_ptr[]: Points to the index in values[] where each row starts.

## Project Structure

- **`spmv.c`**: The main driver file for the program. It initializes the matrices, performs computations, and compares performance between dense and sparse matrix-vector multiplication.
- **`my_dense.c`**: Contains the custom dense matrix-vector multiplication function (`my_dense()`).
- **`my_sparse.c`**: Contains the function to convert a dense matrix to GSL’s sparse matrix format (CSR) and the memory management functions for sparse matrices.
- **`spmv.h`**: Header file that contains the function declarations and includes necessary GSL headers.
- **`timer.c`**: Contains utility functions to measure time with nanosecond precision.
- **`Makefile`**: Automates the build process.
