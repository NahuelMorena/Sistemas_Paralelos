#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define DEBUG 0

//Declaración de funciones 

int getRandomNumber();

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

    //indices
    int i, j, k, offI, offJ;

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

    //Operar matrices
    minA = maxA = A[0];
    minB = maxB = B[0];
    #pragma omp parallel 
    {
        #pragma omp for reduction(max:maxA) reduction(min:minA) reduction(+:promA) nowait   
        for(i=0;i<N;i++) {
            offI = i * N;
            for(j=0;j<N;j++) {
	            double valor = A[offI+j];
                if (valor < minA) {
                    minA = valor;
                }
                if (valor > maxA) {
                    maxA = valor;
                }
                promA += valor;
            }           
        }
        
        #pragma omp for reduction (max:maxB) reduction(min:minB) reduction(+:promB) nowait 
        for(j=0;j<N;j++) {
            offJ = j * N;
            for(i=0;i<N;i++) {
	            double valor = B[i+offJ];
                if (valor < minB) {
                    minB = valor;
                }
                if (valor > maxB) {
                    maxB = valor;
                }
                promB += valor;
            }
        } 
    }

    promA = promA/(size);
    promB = promB/(size);

    double scalar = (maxA * maxB - minA * minB) / (promA * promB);

    printf("valor del escalar: %f",scalar);
}

/*****************************************************************/

int getRandomNumber(){
    return rand() % 41 + 1;
}

/*****************************************************************/