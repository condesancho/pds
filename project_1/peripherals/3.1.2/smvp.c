#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

uint32_t* smvp(uint32_t *values, uint32_t *rowIndex, uint32_t *colIndex, uint32_t *v, uint32_t size);

int main (void){
    uint32_t n=4;
    uint32_t *C = (uint32_t*)calloc(n, sizeof(uint32_t));
    uint32_t row[] = { 2, 3, 0, 3, 0, 2 };
    int col[] = { 0, 2, 2, 4, 6 };
    uint32_t val[] = {1, 1, 1, 1, 1, 1};

    uint32_t v[] = {1, 1, 1, 1};

    C = smvp(val, row, col, v, n);

    for(uint32_t i=0; i<4; i++)
        printf("%d ", C[i]);
    printf("\n");

    return 0;
}

// Sparse Matrix-Vector Product
// We assume that the matrix in CSC form is square
uint32_t* smvp(uint32_t *values, uint32_t *rowIndex, uint32_t *colIndex, uint32_t *v, uint32_t size){

    uint32_t* C = (uint32_t*) malloc(size * sizeof(int32_t));

    for (uint32_t i=0; i<size; i++){
        for (uint32_t k=colIndex[i]; k<colIndex[i+1]; k++){
            C[i] += values[k]*v[rowIndex[k]];
        }
    }
    return C;
}