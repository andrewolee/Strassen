#include <stdio.h>
#include <stdlib.h>

#define ADJACENCY 0
#define RANDOM_1 1

int main (int argc, char* argv[]) {
    int mode = atoi(argv[1]);
    int n = atoi(argv[2]);

    FILE* fp = fopen("test.txt", "w");

    switch (mode) {
        case ADJACENCY: {
            float p = atof(argv[3]);
            int A[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    if (i == j) {
                        A[i][j] = 0;
                        continue;
                    }
                    if (rand() < p * RAND_MAX) {
                        A[i][j] = A[j][i] = 1;
                    } else {
                        A[i][j] = A[j][i] = 0;
                    }
                    
                }
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    fprintf(fp, "%i\n", A[i][j]);
                }
            }
            break;
        }
        case RANDOM_1: {
            for (int i = 0; i < n * n * 2; i++) {
                if (rand() < RAND_MAX / 2) {
                    fprintf(fp, "0\n");
                } else {
                    fprintf(fp, "1\n");
                }
            }
            break;
        }
    }

    fclose(fp);
}