#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

#define N_0 512 // Crossover point
#define BUFFER_SIZE 40

// Print diagonal of matrix
void print (int n, int A[][n]) {
    for (int i = 0; i < n; i++) {
        printf("%i\n", A[i][i]);
    }
}

// Print full matrix
void print_full (int n, int A[][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%i ", A[i][j]);
        }
        printf("\n");
    }
}

// Standard matrix multiplication algorithm
void standard (int n, int A[][n], int B[][n], int C[][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Helpers for Strassen's
void add (int n, int A[][n], int B[][n], int C[][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void sub (int n, int A[][n], int B[][n], int C[][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// Copy part of larger square matrix to smaller square matrix
void split_cpy (int n_a, int A[][n_a], int i_a, int j_a, int n_b, int B[][n_b]) {
    for (int i = 0; i < n_b; i++) {
        for (int j = 0; j < n_b; j++) {
            if (i + i_a < n_a && j + j_a < n_a) {
                B[i][j] = A[i + i_a][j + j_a];
            } else {
                B[i][j] = 0;
            }
        }
    }
}

// Copy smaller square marix to part of larger square matrix
void stitch_cpy (int n_a, int A[][n_a], int n_b, int B[][n_b], int i_b, int j_b) {
    for (int i = 0; i < n_a; i++) {
        for (int j = 0; j < n_a; j++) {
            if (i + i_b < n_b && j + j_b < n_b) {
                B[i + i_b][j + j_b] = A[i][j];
            }
        }
    }
}

void strassen (int n, int A[][n], int B[][n], int C[][n]) {
    if (n <= N_0) {
        standard(n, A, B, C);
        return;
    }

    int m = (n + 1) / 2;

    int (*M)[m][m] = malloc(sizeof(int) * 8 * m * m);
    
    split_cpy(n, A, 0, 0, m, M[0]);
    split_cpy(n, A, 0, m, m, M[1]);
    split_cpy(n, A, m, 0, m, M[2]);
    split_cpy(n, A, m, m, m, M[3]);

    split_cpy(n, B, 0, 0, m, M[4]);
    split_cpy(n, B, 0, m, m, M[5]);
    split_cpy(n, B, m, 0, m, M[6]);
    split_cpy(n, B, m, m, m, M[7]);


    int (*P)[m][m] = malloc(sizeof(int) * 7 * m * m);
    int temp[2][m][m];

    // P1 = A(F - H)
    sub(m, M[5], M[7], temp[0]);
    strassen(m, M[0], temp[0], P[0]);

    // P2 = (A + B)H
    add(m, M[0], M[1], temp[0]);
    strassen(m, temp[0], M[7], P[1]);

    // P3 = (C + D)E
    add(m, M[2], M[3], temp[0]);
    strassen(m, temp[0], M[4], P[2]);

    // P4 = D(G - E)
    sub(m, M[6], M[4], temp[0]);
    strassen(m, M[3], temp[0], P[3]);

    // P5 = (A + D)(E + H)
    add(m, M[0], M[3], temp[0]);
    add(m, M[4], M[7], temp[1]);
    strassen(m, temp[0], temp[1], P[4]);

    // P6 = (B - D)(G + H)
    sub(m, M[1], M[3], temp[0]);
    add(m, M[6], M[7], temp[1]);
    strassen(m, temp[0], temp[1], P[5]);

    // P7 = (A - C)(E + F)
    sub(m, M[0], M[2], temp[0]);
    add(m, M[4], M[5], temp[1]);
    strassen(m, temp[0], temp[1], P[6]);

    free(M);

    // AE + BG
    add(m, P[4], P[3], temp[0]);
    sub(m, temp[0], P[1], temp[1]);
    add(m, temp[1], P[5], temp[0]);
    stitch_cpy(m, temp[0], n, C, 0, 0);

    // AF + BH
    add(m, P[0], P[1], temp[0]);
    stitch_cpy(m, temp[0], n, C, 0, m);

    // CE + DG
    add(m, P[2], P[3], temp[0]);
    stitch_cpy(m, temp[0], n, C, m, 0);

    // CF + DH
    add(m, P[4], P[0], temp[0]);
    sub(m, temp[0], P[2], temp[1]);
    sub(m, temp[1], P[6], temp[0]);
    stitch_cpy(m, temp[0], n, C, m, m);

    free(P);
}

// Timing functions
void start_time(struct timeval *start) {
    gettimeofday(start, NULL); 
}

void stop_time(struct timeval *stop) {
    gettimeofday(stop, NULL);
}

void print_time(struct timeval start, struct timeval stop) {
    double secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + 
        (double)(stop.tv_sec - start.tv_sec); 
    printf("%f seconds\n", secs);
}

int main (int argc, char* argv[]) {
    int mode = atoi(argv[1]);

    if (argc != 4) {
        return -1;
    }
    
    int n = atoi(argv[2]);
    char* file = argv[3];

    int (*A)[n] = malloc(sizeof(int) * n * n);
    int (*B)[n] = malloc(sizeof(int) * n * n);
    int (*C)[n] = malloc(sizeof(int) * n * n);

    // Read in file, and fill matrices A and B
    FILE* fp = fopen(file, "r");
    char buf[BUFFER_SIZE];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            assert(fgets(buf, BUFFER_SIZE, fp) != NULL);
            A[i][j] = atoi(buf);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            assert(fgets(buf, BUFFER_SIZE, fp) != NULL);
            B[i][j] = atoi(buf);
        }
    }
    fclose(fp);

    switch (mode) {
        case 0: {
            strassen(n, A, B, C);
            print(n, C);
            break;
        }
        case 1: {
            struct timeval start, stop;
            /*
            start_time(&start);
            standard(n, A, B, C);
            stop_time(&stop);
            printf("Standard Algorithm: ");
            print_time(start, stop);*/

            start_time(&start);
            strassen(n, A, B, C);
            stop_time(&stop);
            printf("n_0 %i: ", N_0);
            print_time(start, stop);
            break;
        }
        case 2: {
            strassen(n, A, B, C);
            strassen(n, B, C, A);
            int triangles = 0;
            for (int i = 0; i < n; i++) {
                triangles += A[i][i];
            }
            triangles /= 6;
            printf("%i Triangles\n", triangles);
        }
    }
    

    free(A);
    free(B);
    free(C);
}
