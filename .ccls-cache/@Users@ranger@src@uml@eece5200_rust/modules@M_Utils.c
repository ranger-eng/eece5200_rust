#define MAX 10
#include <stdio.h>
#include <stdlib.h>

// gess, svd, and eigen from port3 library, and svdc.f
void eigen_(int *, int *, float *, float *, float *, float *);
void svd_(int *, int *, int *, float *, float *, int *, float *, int *, float *, int *, float *);

//void M_Print(int, int, float [][]); 
void M_Print(int M_n, int M_m, float M[][M_m]) {
    printf("[");
    for(int i=0; i<M_n; i++) {
        for(int j=0; j<M_m; j++) {
           printf("[%f],",M[i][j]); 
        }
        printf("\n");
    }
    printf("]\n");
    }

//void M_Clear(float**, int, int); 
void M_Clear(int M_n, int M_m, float M[][M_m]) {
    for(int i=0; i<M_n; i++) {
        for(int j=0; j<M_m; j++) {
           M[i][j] = 0.0;
        }
    }
    
}

//void M_Transpose(float**, float**, int, int);
void M_Transpose(int M_n, int M_m, float M[][M_m], float Mt[][M_n]) {
    float value;
    for(int i=0; i<M_n; i++) {
        for(int j=0; j<M_m; j++) {
            value = M[i][j];
            Mt[i][j] = value;
        }
    }

}

//void M_Multiply(float**, float**, float**, int, int, int);
void M_Multiply(int A_n, int A_m_B_n, int B_m, float A[][A_m_B_n], float B[][B_m], float C[][B_m]) {
    for(int i=0; i<A_n; i++) {
        for(int j=0; j<A_m_B_n; j++) {
            for(int k=0; k<B_m; k++) {
                C[i][k] = C[i][k] + A[i][j]*B[j][k];            
            }
        }
    }
}

void gess_( int *, float *, int *, float *, int * , int *, float * );
//float GESS(float **, float **, float **, int, int);
float GESS(int A_n, int B_m, float A[][A_n], float X[][B_m], float B[][B_m]) {
    float cond = 0.0;
    // malloc an 2-d array that stores fortran style memory
    float* A_4tran = (float*)malloc(A_n * A_n * sizeof(float));
    float* X_4tran = (float*)malloc(A_n * B_m * sizeof(float));
    float* B_4tran = (float*)malloc(A_n * B_m * sizeof(float));

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<A_n; j++) {
            *(A_4tran+(i)+(j)*A_n) = A[i][j];
        }
    }

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<B_m; j++) {
            *(X_4tran+(i)+(j)*A_n) = X[i][j];
        }
    }

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<B_m; j++) {
            *(B_4tran+(i)+(j)*A_n) = B[i][j];
        }
    }

    gess_((int *)&A_n, (float *)A_4tran, (int *)&A_n, (float *)B_4tran, (int *)&A_n, (int *)&B_m, (float *)&cond);
    

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<B_m; j++) {
            X[i][j] = *(B_4tran+(i)+(j)*A_n);
        }
    }


    free(B_4tran);
    free(X_4tran);
    free(A_4tran);

    return cond;
}

void eigen_(int *, int *, float *, float *, float * , float *);
float EIGEN(int A_n, float A[][A_n], float EVALreal[][1], float EVALimaginary[][1],float EVEC[][A_n*2]) {
    // malloc an 2-d array that stores fortran style memory
    float* A_4tran = (float*)malloc(A_n * A_n * sizeof(float));
    float* EVALreal_4tran = (float*)malloc(A_n * 1 * sizeof(float));
    float* EVALimaginary_4tran = (float*)malloc(A_n * 1 * sizeof(float));
    float* EVEC_4tran = (float*)malloc(A_n * A_n*2 * sizeof(float));

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<A_n; j++) {
            *(A_4tran+(i)+(j)*A_n) = A[i][j];
        }
    }

    // do something
    eigen_((int *)&A_n, (int *)&A_n, (float *)A_4tran, (float *)EVALreal_4tran, (float *)EVALimaginary_4tran , (float *)EVEC_4tran);

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<1; j++) {
            EVALreal[i][j] = *(EVALreal_4tran+(i)+(j)*A_n);
        }
    }

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<1; j++) {
            EVALimaginary[i][j] = *(EVALimaginary_4tran+(i)+(j)*A_n);
        }
    }

    for (int i=0; i<A_n; i++) {
        for (int j=0; j<A_n*2; j++) {
            EVEC[i][j] = *(EVEC_4tran+(i)+(j)*A_n);
        }
    }

    for (int i=0; i<A_n*A_n*2; i++) {
        printf("%f\n",*(EVEC_4tran+i));
    }

    free(A_4tran);
    free(EVALreal_4tran);
    free(EVALimaginary_4tran);
    free(EVEC_4tran);
}

//TODO
//void sdv_(int *, int *, int *, float *, float *, int *, float *, int *, float *, int *, float *);
//float SVD(int A_n, int A_m) {
//}
//
void Hey(void);

void Hey() {
    printf("Hello World");
}
