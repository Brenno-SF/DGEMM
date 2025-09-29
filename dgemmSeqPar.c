#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define blockSize 64

//type this in linux terminal "gcc -o teste -Wall -O3 -fopenmp -march=native -mfma dgemmSeqParGPT.c -lm" to compile

void randomMatrix(double *aa, int nn){
    int i;
    for (i = 0; i < nn; i++)
        aa[i] = (double)(rand()) / RAND_MAX;
}

double calculateDiff(double *seq, double *par, int z) {
    double diff = 0.0;
    for (int i = 0; i < z; i++) {
        diff += fabs(seq[i] - par[i]);
    }
    return diff;
}

void transposeSeq(int n, double *restrict matrix, double *restrict result) {
    int i, j, ii, jj;
    //block tiling to improve cache performance
    for (ii = 0; ii < n; ii+=blockSize){
        for(jj =0; jj < n; jj+=blockSize){

            for (i = ii; i < ii + blockSize && i < n ; i++) {
                for (j = jj; j < jj + blockSize && j < n; j++) {
                    result[j * n + i] = matrix[i * n + j];
                }

            }
        }
    }
}

void transposePar(int n, double *restrict matrix, double *restrict result, int numThreads) {
    int i, j, ii, jj;
    //block tiling to improve cache performance
    //colapse(2) to parallelize both loops at the same time
    #pragma omp parallel for collapse(2) private(i,j, ii, jj) num_threads(numThreads)
    for (ii = 0; ii < n; ii+=blockSize){
        for(jj =0; jj < n; jj+=blockSize){

            for (i = ii; i < ii + blockSize && i < n ; i++) {
                for (j = jj; j < jj + blockSize && j < n; j++) {
                    result[j * n + i] = matrix[i * n + j];
                }
            }
        }
    }
}

void dgemmSeq(int n, double alpha, double *restrict a, double *restrict b, double *restrict bT,double beta, double *restrict c){

    int i, j, k, ii, jj, kk;  

    transposeSeq(n, b, bT); // Transpose matrix B to improve performance

    //block tiling to improve cache performance
    for (ii = 0; ii < n; ii += blockSize) {
        for (jj = 0; jj < n; jj += blockSize) {
            for (kk = 0; kk < n; kk += blockSize) {   // <<< tiling no eixo k
                for (i = ii; i < ii + blockSize && i < n; i += 2) {
                    for (j = jj; j < jj + blockSize && j < n; j += 2) {
                        //register blocking to improve performance
                        double c00 = 0.0, c01 = 0.0, c10 = 0.0, c11 = 0.0;

                        for (k = kk; k < kk + blockSize && k < n; k++) {
                            double a0k = a[i * n + k];
                            double a1k = (i+1 < n) ? a[(i+1) * n + k] : 0.0;
                            double b0k = bT[j * n + k];
                            double b1k = (j+1 < n) ? bT[(j+1) * n + k] : 0.0;

                            c00 += a0k * b0k;
                            c01 += a0k * b1k;
                            c10 += a1k * b0k;
                            c11 += a1k * b1k;
                        }

                        c[i * n + j]     += alpha * c00;
                        if (j+1 < n) c[i * n + j+1] += alpha * c01;
                        if (i+1 < n) c[(i+1) * n + j] += alpha * c10;
                        if (i+1 < n && j+1 < n) c[(i+1) * n + j+1] += alpha * c11;
                    }
                }
            }
        }
    }

}

void dgemmPar(int n, double alpha, double *restrict a, double *restrict b, double *restrict bT,double beta, double *restrict c, int numThreads){
    int i, j, k, ii, jj, kk;       

    transposePar(n, b, bT, numThreads); // Transpose matrix B to improve  performance
    
    //block tiling to improve cache performance
    //colapse(2) to parallelize both loops at the same time
    #pragma omp parallel for collapse(2) private(i,j,k,ii,jj,kk) num_threads(numThreads)
    for (ii = 0; ii < n; ii += blockSize) {
        for (jj = 0; jj < n; jj += blockSize) {
            for (kk = 0; kk < n; kk += blockSize) {   // <<< tiling no eixo k
                for (i = ii; i < ii + blockSize && i < n; i += 2) {
                    for (j = jj; j < jj + blockSize && j < n; j += 2) {
                        //register blocking to improve performance
                        double c00 = 0.0, c01 = 0.0, c10 = 0.0, c11 = 0.0;

                        for (k = kk; k < kk + blockSize && k < n; k++) {
                            double a0k = a[i * n + k];
                            double a1k = (i+1 < n) ? a[(i+1) * n + k] : 0.0;
                            double b0k = bT[j * n + k];
                            double b1k = (j+1 < n) ? bT[(j+1) * n + k] : 0.0;

                            c00 += a0k * b0k;
                            c01 += a0k * b1k;
                            c10 += a1k * b0k;
                            c11 += a1k * b1k;
                        }

                        c[i * n + j]     += alpha * c00;
                        if (j+1 < n) c[i * n + j+1] += alpha * c01;
                        if (i+1 < n) c[(i+1) * n + j] += alpha * c10;
                        if (i+1 < n && j+1 < n) c[(i+1) * n + j+1] += alpha * c11;
                    }
                }
            }
        }
    }

}

int main(){
    int n, z, numThreads = 2;
    double alpha = 2.0, beta = 0.0; // arbitrary values
    double *a, *b, *cSeq, *cPar2, *cPar4, *cPar8, *cPar12;

    printf("Matrix size:  ");

    if (scanf("%d", &n) != 1) {
        printf("Error: invalid input.\n");
        return 1; 
    }

    z = n * n;

    a = (double *)calloc(z, sizeof(double));
    b = (double *)calloc(z, sizeof(double));

    //Calloc to intialize c matrix equal to zero 
    cSeq = (double *)calloc(z, sizeof(double)); 
    cPar2 = (double *)calloc(z, sizeof(double));
    cPar4 = (double *)calloc(z, sizeof(double));
    cPar8 = (double *)calloc(z, sizeof(double));
    cPar12 = (double *)calloc(z, sizeof(double));

    double *bT = (double *)calloc(z, sizeof(double)); // Transposed matrix B
    printf("Generate random matrix 1 (%d x %d)\n", n, n);
    randomMatrix(&a[0], z);
    printf("Generate random matrix 2 (%d x %d)\n", n, n);
    randomMatrix(&b[0], z);

    double t0 = omp_get_wtime();
    dgemmSeq(n, alpha, a, b, bT, beta, cSeq);
    double t1 = omp_get_wtime();
    double timeDgemmSeq = t1 - t0;
    printf("dgemmSeq time: %f seconds\n", timeDgemmSeq);

    t0 = omp_get_wtime();
    dgemmPar(n, alpha, a, b, bT, beta, cPar2, numThreads);
    t1 = omp_get_wtime();
    double timeDgemmPar2 = t1 - t0;
    printf("dgemmPar with 2 threads time: %f seconds\n", timeDgemmPar2);
    
    numThreads = 4;
    t0 = omp_get_wtime();
    dgemmPar(n, alpha, a, b, bT, beta, cPar4, numThreads);
    t1 = omp_get_wtime();
    double timeDgemmPar4 = t1 - t0;
    printf("dgemmPar with 4 threads time: %f seconds\n", timeDgemmPar4);
 
    numThreads = 8;
    t0 = omp_get_wtime();
    dgemmPar(n, alpha, a, b, bT, beta, cPar8, numThreads);
    t1 = omp_get_wtime();
    double timeDgemmPar8 = t1 - t0;
    printf("dgemmPar with 8 threads time: %f seconds\n", timeDgemmPar8);

    numThreads = 12;
    t0 = omp_get_wtime();
    dgemmPar(n, alpha, a, b, bT, beta, cPar12, numThreads);
    t1 = omp_get_wtime();
    double timeDgemmPar12 = t1 - t0;
    printf("dgemmPar with 12 threads time: %f seconds\n", timeDgemmPar12);
    
    // Calculate speedup
    double speedup2 = timeDgemmSeq / timeDgemmPar2;
    printf("\nSpeedup with 2 threads: %f\n", speedup2);

    double speedup4 = timeDgemmSeq / timeDgemmPar4;
    printf("Speedup with 4 threads: %f\n", speedup4);
    
    double speedup8 = timeDgemmSeq / timeDgemmPar8;
    printf("Speedup with 8 threads: %f\n", speedup8);

    double speedup12 = timeDgemmSeq / timeDgemmPar12;
    printf("Speedup with 12 threads: %f\n", speedup12);

    // Calculate efficiency
    double efficiency2 = speedup2 / 2.0;
    printf("\nEfficiency with 2 threads: %f\n", efficiency2);

    double efficiency4 = speedup4 / 4.0;
    printf("Efficiency with 4 threads: %f\n", efficiency4);

    double efficiency8 = speedup8 / 8.0;
    printf("Efficiency with 8 threads: %f\n", efficiency8);

    double efficiency12 = speedup12 / 12.0;
    printf("Efficiency with 12 threads: %f\n", efficiency12);

    // Verify correctness by calculating the difference between sequential and parallel results
    printf("\nCalculating differences between sequential and parallel results...\n");
    printf("Difference with 2 threads: %e\n", calculateDiff(cSeq, cPar2, z));
    printf("Difference with 4 threads: %e\n", calculateDiff(cSeq, cPar4, z));
    printf("Difference with 8 threads: %e\n", calculateDiff(cSeq, cPar8, z));
    printf("Difference with 12 threads: %e\n", calculateDiff(cSeq, cPar12, z));

    // Free allocated memory
    free(a);
    free(b);
    free(cSeq);
    free(cPar2);
    free(cPar4);
    free(cPar8);
    free(cPar12);
    free(bT);

    return 0;
}