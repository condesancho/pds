#include <stdio.h>
#include <stdint.h>
#include <stdlib.h> 
#include <time.h>
#include "mmio.h"

#include <omp.h>

// My functions
int findValue(uint32_t *arr, uint32_t size, uint32_t toSearch);

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);


/* ..... The main ..... */

int main(int argc, char* argv[]){

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t i, *I, *J;
    double *val;

    if (argc < 3){
		fprintf(stderr, "Usage: %s [martix-market-filename] [number-of-threads]\n", argv[0]);
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

    I = (uint32_t *) malloc(nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode)){
        for (i=0; i<nz; i++){
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
        }
    }
    else{
        for (i=0; i<nz; i++){
        fscanf(f, "%d %d\n", &I[i], &J[i]);
        val[i]=1;
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
        }
    }

    if (f !=stdin) fclose(f);

    // Checks if the matrix is square
    if (M!=N){
        printf("Matrix is not square\n");
        exit(1);
    }


    const uint32_t nnz = nz;
    const uint32_t n   = M;

    // starting matrices
    uint32_t * row = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
    uint32_t * col = (uint32_t *)malloc((n + 1) * sizeof(uint32_t));

    uint32_t isOneBased = 0;

    // Transform the coo I and J matrices to csc form
    coo2csc(row, col, I, J, nnz, n, isOneBased);


    // Array that stores the # of triangles a node belongs
    uint32_t c3[n];
    for (uint32_t l = 0; l < n; l++) c3[l] = 0;

    struct timespec begin, end;

    // Starting the clock
    clock_gettime(CLOCK_REALTIME, &begin);

    int num_of_threads = atoi(argv[2]);
    omp_set_num_threads(num_of_threads);

    #pragma omp parallel
    {
        for (uint32_t i=0; i<n-2; i++){

            // Array that contains the adjacent nodes to i
            uint32_t adj_i_size = col[i+1]-col[i];
            uint32_t adj_i[adj_i_size];
            for (uint32_t l = 0; l < adj_i_size; l++) adj_i[l] = 0;
            
            for (uint32_t l=0; l<adj_i_size; l++){
                adj_i[l] = row[col[i]+l];
            }
            #pragma omp for nowait
            // Runs for every j node that is adjacent to i
            for (uint32_t j=0; j<adj_i_size; j++){

                // Array that contains the adjacent nodes to j
                uint32_t adj_j_size = col[adj_i[j]+1]-col[adj_i[j]];
                uint32_t adj_j[adj_j_size];

                for (uint32_t l = 0; l < adj_j_size; l++) adj_j[l] = 0;

                for (uint32_t l=0; l<adj_j_size; l++){
                    adj_j[l] = row[col[adj_i[j]]+l];
                }


                // Runs for every k node that is adjacent to j        
                for (uint32_t k=0; k<adj_j_size; k++){
                    if (i<adj_i[j] && adj_i[j]<adj_j[k] && findValue(adj_i, adj_i_size, adj_j[k])){
                        #pragma omp critical
                        {
                            c3[i]++;
                            c3[adj_i[j]]++;
                            c3[adj_j[k]]++;
                        }
                    }
                }
            }
        }
    }

    uint32_t sum = 0;
    for (uint32_t l=0; l<n; l++){        
        sum += c3[l];
    }
    printf("# of trianges: %d\n", sum/3);

    // Stopping the clock
    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    printf("Time elapsed: %.5f seconds.\n", elapsed);

    /* cleanup variables */
    free( row );
    free( col );
    free ( I );
    free ( J );
    free ( val );

    return 0;
}


/* ..... The functions that are used ..... */

// Returns 1 if k is a value of adj[]
int findValue(uint32_t *arr, uint32_t size, uint32_t toSearch){
    int found = 0;
    
    for(uint32_t i = 0; i < size; i++){
        /* 
         * If element is found in array then raise found flag
         * and terminate from loop.
         */
        if(arr[i] == toSearch){
            found = 1;
            return found;
        }
    }
    return found;
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
