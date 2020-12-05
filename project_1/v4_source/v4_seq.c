#include <stdio.h>
#include <stdint.h>
#include <stdlib.h> 
#include <time.h>
#include "mmio.h"

uint32_t commonElements(uint32_t* A, uint32_t* B, uint32_t size_A, uint32_t size_B);

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

int main(int argc, char* argv[]){

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t i, *I, *J;
    double *val;

    if (argc < 2){
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else{
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) ){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    I = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
    val = (double *) malloc(2*nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    int offset = 0;

    if (!mm_is_pattern(matcode)){
        for (i=0; i<nz; i++){
            fscanf(f, "%d %d %lg\n", &I[i+offset], &J[i], &val[i+offset]);
            I[i+offset]--;  /* adjust from 1-based to 0-based */
            J[i+offset]--;
            I[i+offset+1] = J[i+offset];
            J[i+offset+1] = I[i+offset];

            offset++;
        }
    }
    else{
        for (i=0; i<nz; i++){
            fscanf(f, "%d %d\n", &I[i+offset], &J[i+offset]);
            val[i]=1;
            I[i+offset]--;  /* adjust from 1-based to 0-based */
            J[i+offset]--;
            I[i+offset+1] = J[i+offset];
            J[i+offset+1] = I[i+offset];

            offset++;
        }
    }


    if (f !=stdin) fclose(f);

    // Checks if the matrix is square
    if (M!=N){
        printf("Matrix is not square\n");
        exit(1);
    }


    const uint32_t nnz = 2*nz;
    const uint32_t n   = M;

    // starting matrices
    uint32_t * row = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
    uint32_t * col = (uint32_t *)malloc((n + 1) * sizeof(uint32_t));

    uint32_t isOneBased = 0;

    // Transform the coo I and J matrices to CSC form
    coo2csc(row, col, I, J, nnz, n, isOneBased);


    // Stores the element values of the product A.*(A*A) in CSC form
    uint32_t *values = (uint32_t*)calloc(nnz, sizeof(uint32_t));

    // Array that stores the # of triangles a node belongs
    uint32_t *c3 = (uint32_t*)calloc(n, sizeof(uint32_t));

    struct timespec begin, end;

    // Starting the clock
    clock_gettime(CLOCK_REALTIME, &begin);
    
    for (uint32_t i=0; i<n; i++){

        // Array that contains the rows in column i that exist
        uint32_t col_i_size = col[i+1]-col[i];
        uint32_t *col_i = (uint32_t*)calloc(col_i_size, sizeof(uint32_t));
        
        for (uint32_t l=0; l<col_i_size; l++){
            col_i[l] = row[col[i]+l];
        }


        /*
            Runs for every element in column i that exists.
            In this way I mask the elements in the A*A matrix that have
            to become zero due to the Hadamard product.
            Since the matrix A is symmetrical the row j == col j
        */
        for (uint32_t j=0; j<col_i_size; j++){

            // Array that contain the rows in column j that exist
            uint32_t col_j_size = col[col_i[j]+1]-col[col_i[j]];
            uint32_t *col_j = (uint32_t*)calloc(col_j_size,sizeof(uint32_t));

            for (uint32_t l=0; l<col_j_size; l++){
                col_j[l] = row[col[col_i[j]]+l];
            }

            // Find the number of common elements in column i and column j
            values[col[i]+j] = commonElements(col_i, col_j, col_i_size, col_j_size);
            
            free ( col_j );
        }

        free ( col_i );
    }

    // This is the multiplication with the vector of size n and elements 1
    for (uint32_t i=0; i<n; i++){
        for (uint32_t k=col[i]; k<col[i+1]; k++){
            c3[i] += values[k];
        }
        c3[i] = c3[i]/2;
    }

    // Stopping the clock
    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    uint32_t sum = 0;
    for (uint32_t l=0; l<n; l++){        
        sum += c3[l];
    }
    printf("# of trianges: %d\n", sum/3);

    printf("Time elapsed: %.5f seconds.\n", elapsed);

    /* cleanup variables */
    free( row );
    free( col );
    free ( I );
    free ( J );
    free ( val );
    free ( values );
    free ( c3 );

    return 0;
}


/* ..... The functions that are used ..... */

// Function that returns the number of common elements in array A and B
// ..... Arrays MUST BE SORTED
uint32_t commonElements(uint32_t* A, uint32_t* B, uint32_t size_A, uint32_t size_B){
    // 2 iterators
    uint32_t i = 0; // Points to an index in array A
    uint32_t j = 0; // Points to an index in array B

    uint32_t common = 0;

    while (i<size_A && j<size_B){
        if (A[i] == B[j]){
            common++;
            i++;
            j++;
        }            
        else if (A[i] < B[j])
            i++;
        else
            j++;
    }
    return common;
}


/**
 *  \brief COO to CSC conversion
 *
 *  Converts a square matrix from COO to CSC format.
 *
 *  Note: The routine assumes the input COO and the output CSC matrix
 *  to be square.
 *
 */
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}