/* Amrith Deepak
 * amrithd@andrew.cmu.edu
 * February 26, 2014
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"
#include "contracts.h"




int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that takes
 * the matrix and transposes it in a way that minimizes misses and evictions.
 */

char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    REQUIRES(M > 0);
    REQUIRES(N > 0);
    if(M == 32) {
        int diagonal = 0, diagvalue, i, j, k, l;
        int block_size = 8;
        for (i = 0; i < N; i += block_size) {
            for (j = 0; j < M; j += block_size) {
                for (k = i; k < (i + block_size); k++) {
                    for (l = j; l < (j + block_size); l++) {
                        if(k != l) {
                            B[l][k] = A[k][l];
                        }
                        else {
                            diagvalue = A[k][l];
                            diagonal = 1;
                        }
                    }
                    if(diagonal == 1) {
                        B[k][k] = diagvalue;
                        diagonal = 0;
                    }
                }
            }
        }
    }
    else if(M==64) {
        int block_size = 8;
        int n;
        int a, b, c, d;
        int i, j, k, l, m; // 11 local variables
        for (i = 0; i < N; i += block_size) {
            for (j = 0; j < M; j += block_size) {
                n = j + 4; //diagonal
                for(k = i + 4; k < i + 8; k++) {
                    if(k == n){ //in-cache trans
                        a = A[k][n];
                        b = A[k][n+1];
                        c = A[k][n+2];
                        d = A[k][n+3];
                        B[k][n] = a;
                        B[k][n+1] = b;
                        B[k][n+2] = c;
                        B[k][n+3] = d;
                        
                        a = A[k+1][n];
                        b = A[k+1][n+1];
                        c = A[k+1][n+2];
                        d = A[k+1][n+3];
                        B[k+1][n] = a;
                        B[k+1][n+1] = b;
                        B[k+1][n+2] = c;
                        B[k+1][n+3] = d;
                        
                        a = A[k+2][n];
                        b = A[k+2][n+1];
                        c = A[k+2][n+2];
                        d = A[k+2][n+3];
                        B[k+2][n] = a;
                        B[k+2][n+1] = b;
                        B[k+2][n+2] = c;
                        B[k+2][n+3] = d;
                        
                        a = A[k+3][n];
                        b = A[k+3][n+1];
                        c = A[k+3][n+2];
                        d = A[k+3][n+3];
                        B[k+3][n] = a;
                        B[k+3][n+1] = b;
                        B[k+3][n+2] = c;
                        B[k+3][n+3] = d;
                        
                        
                        //in-cache-trans
                        for(l = k; l < k + 4; l++){
                            for(m = n; m <n + 4; m++){
                                if(l == m)
                                    break;
                                else {
                                    int temp = B[m][l];
                                    B[m][l]=B[l][m];
                                    B[l][m]=temp;
                                }
                            }
                        }
                        
                        break;
                        
                    }
                    else{
                        a = A[k][n];
                        b = A[k][n+1];
                        c = A[k][n+2];
                        d = A[k][n+3];
                        B[n][k] = a;
                        B[n+1][k] = b;
                        B[n+2][k] = c;
                        B[n+3][k] = d;
                    }
                }
                n = j+4;
                for(k = i; k < i + 4; k++) {
                    a = A[k][n];
                    b = A[k][n+1];
                    c = A[k][n+2];
                    d = A[k][n+3];
                    B[n][k] = a;
                    B[n+1][k] = b;
                    B[n+2][k] = c;
                    B[n+3][k] = d;
                }
                n = j;//diagonal
                for(k = i; k < i + 4; k++){
                    if(k == n){//in-cache trans
                        a = A[k][n];
                        b = A[k][n+1];
                        c = A[k][n+2];
                        d = A[k][n+3];
                        B[k][n] = a;
                        B[k][n+1] = b;
                        B[k][n+2] = c;
                        B[k][n+3] = d;
                        
                        a = A[k+1][n];
                        b = A[k+1][n+1];
                        c = A[k+1][n+2];
                        d = A[k+1][n+3];
                        B[k+1][n] = a;
                        B[k+1][n+1] = b;
                        B[k+1][n+2] = c;
                        B[k+1][n+3] = d;
                        
                        a = A[k+2][n];
                        b = A[k+2][n+1];
                        c = A[k+2][n+2];
                        d = A[k+2][n+3];
                        B[k+2][n] = a;
                        B[k+2][n+1] = b;
                        B[k+2][n+2] = c;
                        B[k+2][n+3] = d;
                        
                        a = A[k+3][n];
                        b = A[k+3][n+1];
                        c = A[k+3][n+2];
                        d = A[k+3][n+3];
                        B[k+3][n] = a;
                        B[k+3][n+1] = b;
                        B[k+3][n+2] = c;
                        B[k+3][n+3] = d;
                        
                        
                        //in-cache-trans
                        for(l = k; l < k + 4; l++){
                            for(m = n; m < n + 4; m++){
                                if(l == m)
                                    break;
                                else {
                                    int temp = B[m][l];
                                    B[m][l]=B[l][m];
                                    B[l][m]=temp;
                                }
                            }
                        }
                        
                        break;
                        
                    }
                    else{
                        a = A[k][n];
                        b = A[k][n+1];
                        c = A[k][n+2];
                        d = A[k][n+3];
                        B[n][k] = a;
                        B[n+1][k] = b;
                        B[n+2][k] = c;
                        B[n+3][k] = d;
                    }
                }
                n = j;
                for(k = i+4; k < i + 8; k++) {
                    a = A[k][n];
                    b = A[k][n+1];
                    c = A[k][n+2];
                    d = A[k][n+3];
                    B[n][k] = a;
                    B[n+1][k] = b;
                    B[n+2][k] = c;
                    B[n+3][k] = d;
                }
            }
        }
    }
    else {
        int i, j, k, l;
        int block_size = 20;
        for (i = 0; i < N; i += block_size) {
            for (j = 0; j < M; j += block_size) {
                for (k = i; k < (i + block_size) && k < N; k++) {
                    for (l = j; l < (j + block_size) && l < M; l++) {
                        B[l][k] = A[k][l];
                    }
                }
            }
        }
    }
    ENSURES(is_transpose(M, N, A, B));
}


/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;
    
    REQUIRES(M > 0);
    REQUIRES(N > 0);
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }
    
    ENSURES(is_transpose(M, N, A, B));
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);
    
    /* Register any additional transpose functions */
    //registerTransFunction(trans, trans_desc);
}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}



