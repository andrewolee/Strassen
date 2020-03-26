#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BUFFER_SIZE 40

void print (int n, int A[][n]) {
    for (int i = 0; i < n; i++) {
        printf("%i\n", A[i][i]);
    }
}

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
    if (n < 30) {
        standard(n, A, B, C);
        return;
    }

    int m = (n + 1) / 2;

    int a[m][m], b[m][m], c[m][m], d[m][m], e[m][m], f[m][m], g[m][m], h[m][m];

    split_cpy(n, A, 0, 0, m, a);
    split_cpy(n, A, 0, m, m, b);
    split_cpy(n, A, m, 0, m, c);
    split_cpy(n, A, m, m, m, d);

    split_cpy(n, B, 0, 0, m, e);
    split_cpy(n, B, 0, m, m, f);
    split_cpy(n, B, m, 0, m, g);
    split_cpy(n, B, m, m, m, h);


    int P[7][m][m];
    int temp[2][m][m];

    // P1 = A(F - H)
    sub(m, f, h, temp[0]);
    strassen(m, a, temp[0], P[0]);

    // P2 = (A + B)H
    add(m, a, b, temp[0]);
    strassen(m, temp[0], h, P[1]);

    // P3 = (C + D)E
    add(m, c, d, temp[0]);
    strassen(m, temp[0], e, P[2]);

    // P4 = D(G - E)
    sub(m, g, e, temp[0]);
    strassen(m, d, temp[0], P[3]);

    // P5 = (A + D)(E + H)
    add(m, a, d, temp[0]);
    add(m, e, h, temp[1]);
    strassen(m, temp[0], temp[1], P[4]);

    // P6 = (B - D)(G + H)
    sub(m, b, d, temp[0]);
    add(m, g, h, temp[1]);
    strassen(m, temp[0], temp[1], P[5]);

    // P7 = (A - C)(E + F)
    sub(m, a, c, temp[0]);
    add(m, e, f, temp[1]);
    strassen(m, temp[0], temp[1], P[6]);

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
}

int main (int argc, char* argv[]) {
    
    int n = atoi(argv[2]);
    char* file = argv[3];

    int A[n][n], B[n][n], C[n][n];

    // Read in file, and fill matrices A and B
    FILE* fp = fopen(file, "r");
    char buf[BUFFER_SIZE];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fgets(buf, BUFFER_SIZE, fp);
            A[i][j] = atoi(buf);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fgets(buf, BUFFER_SIZE, fp);
            B[i][j] = atoi(buf);
        }
    }
    fclose(fp);

    strassen(n, A, B, C);
    print(n, C);
    /*
    int n = 3;
    int A[3][3] = {{1,2,3},{5,6,7},{9,10,11}};
    int B[3][3] = {{1,2,2},{3,4,7},{2,4,5}};
    int C[n][n];

    print_full(n, A);
    print_full(n, B);

    strassen(n, A, B, C);
    printf("Strassen\n");
    print_full(n, C);

    standard(n, A, B, C);
    printf("Standard\n");
    print_full(n, C);*/
}
