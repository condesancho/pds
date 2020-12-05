#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 

int main(int argc, char* argv[]){
    struct timespec begin, end;

    // Starting the clock
    clock_gettime(CLOCK_REALTIME, &begin);

    int n=0;

    if (argc!=2){
        printf("Usage: %s n\nwhere n is size of matrix.\n", argv[0]);
        return 1;
    }

    n = atoi(argv[1]);

    if (n<1){
        printf("Must enter positive number.\n");
        return 1;
    }

    int adj[n][n];
    srand(time(0));

    // Fill random symmetric matrix with zeros and ones and set the diagonal zero
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (j>i)
                adj[i][j] = rand() % 2;
            else if (i==j)
                adj[i][j] = 0;
            else
                adj[i][j] = adj[j][i];
        }
    }

    // Print the matrix I made
    /*
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            printf("%d\t", adj[i][j]);
        }
        printf("\n");
    }
    */
    
    
    int c3[n];
    for (int i=0; i<n; i++){
        c3[i] = 0;
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                if (adj[i][j] == adj[j][k] && adj[j][k] == adj[k][i] && adj[i][j] == 1){
                    c3[i]++;
                    c3[j]++;
                    c3[k]++;
                }
            }
        }
    }

    /*
    printf("\n");

    // Prints the number of triangles adjacent to a node x6
    for (int i=0; i<n; i++){
        printf("%d\t", c3[i]);
    }
    printf("\n");
    */

    // Stopping the clock
    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    printf("Time elapsed: %.5f seconds.\n", elapsed);

    return 0;
} 