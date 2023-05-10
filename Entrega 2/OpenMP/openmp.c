#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>

#define DEBUG 0

//Declaración de funciones 

int getRandomNumber();
double dwalltime();
void blkmulWithEscalar(double *ablk, double *bblk, double *cblk, int n, int bs, double scalar);
void blkmulwithIntMat(double *ablk, int *bblk, double *cblk, int n, int bs);

int main(int argc, char*argv[]){
    //Validar parametros
    if (argc < 4){
        fprintf(stderr, "usage %s matrix_size block_size num_threads\n",argv[0]);
        exit(1);
    }

    //Tomar argumentos del main 
    int N = atoi(argv[1]);
    int bs = atoi(argv[2]);
    int num_threads = atoi(argv[3]);

    #if DEBUG == 1
        printf("N = %i\n",N);
        printf("bs = %i\n",bs);
        printf("NUM_THREADS = %i\n", num_threads);
    #endif


    //Asignar matrices (Se utiliza alternativa 4 mostrada en la teoria)
    double *A,*B,*C,*R,*CxD;
    int *D;

    //declaración de variables
    int i, j, k, offI, offJ;
    double timetick, endtime;
    double item;
    int size = N*N;

    A=(double*)malloc(sizeof(double)*size);
    B=(double*)malloc(sizeof(double)*size);
    C=(double*)malloc(sizeof(double)*size);
    R=(double*)malloc(sizeof(double)*size);
    D=(int*)malloc(sizeof(int)*size);
    CxD=(double*)malloc(sizeof(double)*size);

    double maxA, minA, promA = 0.0;
    double maxB, minB, promB = 0.0;

    //Inicializar matrices 
    srand(time(NULL));

    for (i=0; i<size; i++){
        A[i] = 1.0;
        B[i] = 1.0;
        C[i] = 1.0;
        R[i] = 0.0;
        CxD[i] = 0.0;
        D[i] = getRandomNumber();
    } 

    //OpenMP
    omp_set_num_threads(num_threads);

    //Empieza a contar el tiempo
    timetick = dwalltime();
    
    //Operar matrices
    minA = maxA = A[0];
    minB = maxB = B[0];
    #pragma omp parallel private(i, j, offI, item)
    {
        #pragma omp for reduction(max:maxA) reduction(min:minA) reduction(+:promA) nowait   
        for(i=0;i<N;i++) {
            int tid = omp_get_thread_num();
            for(j=0;j<N;j++) {
	            item = A[offI+j];
                if (item < minA) {
                    minA = item;
                }
                if (item > maxA) {
                    maxA = item;
                }
                promA += item;
            }           
        }
        
        #pragma omp for reduction (max:maxB) reduction(min:minB) reduction(+:promB) 
        for(j=0;j<N;j++) {
            offJ = j * N;
            for(i=0;i<N;i++) {
	            item = B[i+offJ];
                if (item < minB) {
                    minB = item;
                }
                if (item > maxB) {
                    maxB = item;
                }
                promB += item;
            }
        } 
    }
    //printf("valor preliminar de promA %f\n",promA);
    //printf("valor preliminar de promB %f\n",promB);
    promA = promA/(size);
    promB = promB/(size);

    double scalar = (maxA * maxB - minA * minB) / (promA * promB);

    #pragma omp parallel private(i, j, k, offI, offJ)
    {   
        //Multiplicación AxB
        //matmulblksWithEscalar(A, B, R, N, bs, escalar);
        #pragma omp for nowait
        for (i = 0; i < N; i += bs){
            offI = i * N;
            for (j = 0; j < N; j += bs){
                offJ = j * N;
                for (k = 0; k < N; k += bs){
                    blkmulWithEscalar(&A[offI + k], &B[offJ + k], &R[offI + j], N, bs, scalar);
                }
            }
        }

        //Pot2(D)
        #pragma omp for
        for(j=0;j<N;j++) {
            offJ = j*N;
            for(i=0;i<N;i++) {
                int v = D[i+offJ];
                D[i+offJ] = v * v;
            }
        }

        //Multiplicación CxD
        #pragma omp for 
        for (i = 0; i < N; i += bs){
            offI = i * N;
            for (j = 0; j < N; j += bs){
                offJ = j * N;
                for (k = 0; k < N; k += bs){
                    blkmulwithIntMat(&C[offI + k], &D[offJ + k], &CxD[offI + j], N, bs);
                }
            }
        }

        //Suma entre (escalar*AxB) + CxD
        #pragma omp for 
        for (i=0; i<N; i++) {
            offI = i * N;
            for (j=0; j<N; j++) {
                R[offI+j] += CxD[offI+j];
            }
        }
    }

    //Detener el tiempo
    endtime = dwalltime();

    printf("Tiempo empleado en segundos %f\n", endtime - timetick);
        
    #if DEBUG == 1
        printf("valor de maxA: %f\n",maxA);
        printf("valor del minA: %f\n",minA);
        printf("valor del maxB: %f\n",maxB);
        printf("valor del minB: %f\n",minB);
        printf("valor del promA: %f\n",promA);
        printf("valor del promB: %f\n",promB);
        printf("valor del escalar: %f\n",scalar);
    #endif
    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(CxD);

    return 0;
}

/*****************************************************************/

int getRandomNumber(){
    return rand() % 41 + 1;
}

/*****************************************************************/

//Multiplicación de matrices y escalar
void blkmulWithEscalar(double *ablk, double *bblk, double *cblk, int n, int bs, double scalar){
    int i, j, k, offI, offJ;    

    for (i = 0; i < bs; i++){
        int offI = i * n;
        for (j = 0; j < bs; j++){
            int offJ = j * n;
            for (k = 0; k < bs; k++){
                cblk[offI + j] += ablk[offI + k] * bblk[offJ + k];
            }
            cblk[offI + j] *= scalar;
        }
    }
}

void blkmulwithIntMat(double *ablk, int *bblk, double *cblk, int n, int bs){
    int i, j, k, offI, offJ;  

    for (i = 0; i < bs; i++){
        offI = i * n;
        for (j = 0; j < bs; j++){
            offJ = j * n;
            for (k = 0; k < bs; k++){
                cblk[offI + j] += ablk[offI + k] * bblk[offJ + k];
            }
        }
    }
}

/*****************************************************************/

/* Funcion para depurar el tiempo de ejecución */
double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}