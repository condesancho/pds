#include <stdio.h>
#include <stdint.h>
#include <stdlib.h> 


uint32_t commonElements(uint32_t* A, uint32_t* B, uint32_t size_A, uint32_t size_B);



int main(int argc, char* argv[]){

    uint32_t nnz=12;
    uint32_t n=5;

    // starting matrices
    uint32_t row[] = {1, 2, 4, 0, 3, 4, 0, 4, 1, 0, 1, 2};
    uint32_t col[] = {0, 3, 6, 8, 9, 12};



    // Stores the element values of the product A.*(A*A) in CSC form
    uint32_t *values = (uint32_t*)calloc(nnz, sizeof(uint32_t));

    // Array that stores the # of triangles a node belongs
    uint32_t *c3 = (uint32_t*)calloc(n, sizeof(uint32_t));

    
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

            printf("\nfor i = %d and j = %d\n", i,j);

            for (uint32_t l=0; l<col_j_size; l++){
                col_j[l] = row[col[col_i[j]]+l];
            }

            // Find the number of common elements in column i and column j
            values[col[i]+j] = commonElements(col_i, col_j, col_i_size, col_j_size);
            printf("%d ", values[col[i]+j]);
        }
    }

    // This is the multiplication with the vector of size n and elements 1
    for (uint32_t i=0; i<n; i++){
        for (uint32_t k=col[i]; k<col[i+1]; k++){
            c3[i] += values[k];
        }
        c3[i] = c3[i]/2;
    }

    for (uint32_t i=0; i<n; i++){
        printf("%d ", c3[i]);
    }
    printf("\n");

    uint32_t sum = 0;
    for (uint32_t l=0; l<n; l++){        
        sum += c3[l];
    }
    printf("# of trianges: %d\n", sum/3);

    return 0;

}


/* ..... The functions that are used ..... */

// Function that returns the number of common elements in array A and B
// ...... Arrays MUST BE SORTED
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