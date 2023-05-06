#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define DEBUG 0

//Declaración de funciones

int getRandomNumber();

void* calculate_R(void* id);

//mutex
pthread_mutex_t promA_lock, promB_lock;

//barrier
pthread_barrier_t barrier;

//Variables compartidas
int block_size_by_threads, N, size;
double promA, promB, maxA, maxB, scalar;
double minA, minB;
double *A,*B,*C,*R,*CxD;
int *D;

int main(int argc, char*argv[]){
    //Validar parametros
    if (argc < 4){
        fprintf(stderr, "usage %s matrix_size block_size num_threads\n",argv[0]);
        exit(1);
    }

    //Tomar argumentos del main 
    N = atoi(argv[1]);
    #if DEBUG == 1
    printf("N = %i\n",N);
    #endif

    int bs = atoi(argv[2]);
    #if DEBUG == 1
    printf("bs = %i\n",bs);
    #endif

    int num_threads = atoi(argv[3]);
    #if DEBUG == 1
    printf("NUM_THREADS = %i\n", num_threads);
    #endif

    //indices
    int i;

    size = N*N;

    A=(double*)malloc(sizeof(double)*size);
    B=(double*)malloc(sizeof(double)*size);
    C=(double*)malloc(sizeof(double)*size);
    R=(double*)malloc(sizeof(double)*size);
    D=(int*)malloc(sizeof(int)*size);
    CxD=(double*)malloc(sizeof(double)*size);

    promA = promB = maxA = maxB = scalar = 0.0;
    minA = minB = 100.0;

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

    //Pthreads

    pthread_t threads[num_threads];
    int ids[num_threads];

    pthread_mutex_init(&promA_lock,NULL);
    pthread_mutex_init(&promB_lock,NULL);

    pthread_barrier_init(&barrier,NULL,num_threads);

    block_size_by_threads = N/num_threads;
    
    //Crear hilos
    for (int i = 0; i < num_threads; i++){
        ids[i] = i;
        pthread_create(&threads[i], NULL, calculate_R, &ids[i]);
    }

    //Espera a que los hilos terminen
    for (int i = 0; i < num_threads; i++){
        pthread_join(threads[i], NULL);
    }

    printf("valor de maxA %f\n", maxA);
    printf("valor de maxB %f\n", maxB);
    printf("valor de minA %f\n", minA);
    printf("valor de minB %f\n", minB);
    printf("valor de promA %f\n", promA);
    printf("valor de promB %f\n", promB);
    printf("valor del escalar %f\n", scalar);

    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(CxD);

    pthread_mutex_destroy(&promA_lock);
    pthread_mutex_destroy(&promB_lock);
    pthread_barrier_destroy(&barrier);

    return 0;
}

/*****************************************************************/

void* calculate_R(void* ptr){
    int* p, id;
    p = (int*) ptr;
    id = *p;

    printf("empieza ejecución hilo %i\n",id);

    int start_block = block_size_by_threads * id;
    int end_block = block_size_by_threads * (id + 1);

    int i, j, offI, offJ;


    //obteniendo minA, maxA y promA
    double localPromA = 0.0;
    double localMaxA, localMinA = A[start_block];
    for(i=start_block ; i<end_block ;i++) {
        offI = i * N;
        for(j=0;j<N;j++) {
	        double valor = A[offI+j];

            if (valor < localMinA) {
                localMinA = valor;
            }

            if (valor > localMaxA) {
                localMaxA = valor;
            }

            localPromA += valor;
        }
    } 

    pthread_mutex_lock(&promA_lock);
    printf("soy hilo %d localMInA= %f, valor de minA= %f\n",id, localMinA, minA);
    if (localMinA < minA){
        minA = localMinA;
    }
    if (localMaxA > maxA){
        maxA = localMaxA;
    }
    promA += localPromA;
    pthread_mutex_unlock(&promA_lock);

    //obteniendo minB, maxB y promB
    double localPromB = 0.0;
    double localMaxB, localMinB = B[start_block];
    for(j=start_block; j<end_block; j++) {
        offJ = j * N;
        for(i=0;i<N;i++) {
	        double valor = B[i+offJ];

            if (valor < localMinB) {
                localMinB = valor;
            }

            if (valor > localMaxB) {
                localMaxB = valor;
            }

            localPromB += valor;
        }
    } 

    printf("soy el hilo %d\n", id);

    pthread_mutex_lock(&promB_lock);
    printf("soy hilo %d localMInB= %f, valor de minB= %f\n",id, localMinB, minB);
    if (localMinB < minB){
        minB = localMinB;
    }
    if (localMaxB > maxB){
        maxB = localMaxB;
    }
    promB += localPromB;
    pthread_mutex_unlock(&promB_lock);

    pthread_barrier_wait(&barrier);
    if (id == 0){
        promA = promA/size;
        promB = promB/size;
        scalar = (maxA * maxB - minA * minB) / (promA * promB);
    }
    pthread_barrier_wait(&barrier);

    //Hasta este punto se calculo es escalar

    pthread_exit(0);
}

/*****************************************************************/

int getRandomNumber(){
    return rand() % 41 + 1;
}

/*****************************************************************/