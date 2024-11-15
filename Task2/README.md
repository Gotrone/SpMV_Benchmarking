
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


## Result 
With : 
- Matrix size: 16384 x 16384 (268435456 elements)
- 26837519 non-zero elements (10.00%)

I don't understand what need to be in the "Ref" column

### Result for GCC

|              | O0 | O2-novec | O3-vec | Ofast-vec | Ref |
| :----------- | :------: | :------: | :------: | :------: | ----: |
| my_dense     |   764   | 355 | 366 | 782 | ? |
| my_coo       |   97   | 87 | 72 | 72 | ? |   
| my_csr       |  85   | 74 | 131 | 74 | ? |
| my_csc       |  105   | 33 | 30 | 31 | ? |

### Result for ICC

|              | O0 | O2-novec | O3-vec | Ofast-vec | Ref |
| :----------- | :------: | :------: | :------: | :------: | ----: |
| my_dense     |   749   | 357 | 157 | 158 | ? |
| my_coo       |   97   | 75 | 73 | 73 | ? |   
| my_csr       |  77   | 72 | 74 | 183 | ? |
| my_csc       |  92   | 28 | 31 | 30 | ? |
