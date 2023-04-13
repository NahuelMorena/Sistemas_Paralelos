#include <stdio.h>
#include <stdlib.h>
//Utilizado para calcular las potencias
#include <math.h>
//Utilizado para contar tiempo empleado
#include <sys/time.h>
//Utilizado para obtener numeros random
#include <time.h>
//Utilizado para asignar espacio en memoria de string
#include <string.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

//Declaración de funciones 
int getValor(double *matriz,int fila,int columna,int N, int orden);
int getRandomNumber();
void matmulblksWithEscalar(double *a, double *b, double *c, int n, int bs, double escalar);
void blkmulWithEscalar(double *ablk, double *bblk, double *cblk, int n, int bs, double escalar);

void matIntmulblks(double *a, int *b, double *c, int n, int bs);
void blkmulwithIntMat(double *ablk, int *bblk, double *cblk, int n, int bs);
double dwalltime();

void recorrerArreglo(double *matriz, int N, int orden);
void recorrerArreglo2(int *matriz, int N, int orden);

void timelog_start();
void timelog_total();
void timelog(const char *desc);

static double stime, ttime = 0.0;

int main(int argc, char*argv[]){
    //Validar parametros
    if (argc < 3){
        fprintf(stderr, "usage %s matrix_size block_size\n",argv[0]);
        exit(1);
    }

    //Tomar parametro n para el tamaño de las matrices NxN
    int N = atoi(argv[1]);
    printf("N = %i\n",N);

    int bs = atoi(argv[2]);
    printf("bs = %i\n",bs);

    
    //Asignar matrices (Se utiliza alternativa 4 mostrada en la teoria)
    double *A,*B,*C,*R,*CxD;
    int *D;

    //indices
    int i, j, k, offI, offJ;
    int fila, columna;

    //variables varias
    int sum;
    int size = N*N;
    double timetick, endtime;

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

    //Empieza a contar el tiempo
    timetick = dwalltime();


    //Operar matrices

    //obteniendo minA, maxA y promA
    timelog_start();
    minA = maxA = A[0];
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
    promA = promA/(size);
    timelog("minA, maxA, promA");

    //obteniendo minB, maxB y promB
    timelog_start();
    minB = maxB = B[0];
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
    promB = promB/(size);
    timelog("minB, maxB, promB");

    timelog_start();
    double escalar = (maxA * maxB - minA * minB) / (promA * promB);
    timelog("e = (maxA * maxB - minA * minB) / (promA * promB)");
    //printf("valor del escalar: %f\n",escalar);

    //Multiplicación AxB

    timelog_start();
    matmulblksWithEscalar(A, B, R, N, bs, escalar);
    timelog("R = (AxB) * e");

    //printf("MATRIZ AxB\n");
    //recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Pot2(D)

    //printf("MATRIZ D\n");
    //recorrerArreglo2(D,N,ORDENXCOLUMNAS);
    
    timelog_start();
    for(int j=0;j<N;j++) {
        offJ = j*N;
        for(int i=0;i<N;i++) {
            int v = D[i+offJ];
            D[i+j*N] = v * v;
        }
    }
    timelog("D = Pot2(D)");

    //printf("MATRIZ D DESPUES DE POTENCIA DE 2\n");
    //recorrerArreglo2(D,N,ORDENXCOLUMNAS);


    //Multiplicación CxD

    timelog_start();
    matIntmulblks(C, D, CxD, N, bs);
    timelog("CxD = CxD");

    //printf("MATRIZ CxD\n");
    //recorrerArreglo(CxD,N,ORDENXCOLUMNAS);

    //Suma entre escalar*AxB + CxD
    timelog_start();
    for (i=0; i<N; i++) {
        offI = i * N;
        for (j=0; j<N; j++) {
            R[offI+j] += CxD[offI+j];
        }
    }
    timelog("R = R + CxD");

    //printf("MATRIZ (escalar *AxB) + CxD\n");
    //recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Finaliza conteo de tiempo
    //endtime = dwalltime() - timetick;
    //printf("Tiempo en segundos %f\n", endtime);

    timelog_total();

    //Obtención del resultado
    /*
    char file[10000];
    bzero(file,10000);
    char temp_str[2000];

    sprintf(temp_str, "Analisis de tiempo obtenidos en matrices de %i. \n",N);
    strcat(file, temp_str);
    sprintf(temp_str, "Tiempo obtenido de %f segundos\n", endtime);
    strcat(file, temp_str);
    sprintf(temp_str, "-----------------------------------------------\n\n");
    strcat(file, temp_str);
    */

    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(CxD);

    return 0;
}

/* Funciones para depurar el tiempo de ejecución */
double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void timelog_start() {
    stime = dwalltime();
}

void timelog_total() {
    printf("Tiempo total: %.02fms\n", ttime * 1000.0);
}

void timelog(const char *desc) {
    double time = dwalltime() - stime;
    ttime += time;

    if (desc)
        printf("Tiempo '%s': %.02fms\n", desc, time * 1000.0);
}

//--------------------------------------------------------------------------------------------------

int getIntValor(int *matriz,int fila,int columna, int N, int orden){
    if(orden==ORDENXFILAS){
        return(matriz[fila*N+columna]);
    }else{
        return(matriz[fila+columna*N]);
    }
}

double getDoubleValor(double *matriz, int fila, int columna, int N, int orden){
    if(orden==ORDENXFILAS){
        return(matriz[fila*N+columna]);
    }else{
        return(matriz[fila+columna*N]);
    }
}

int getRandomNumber(){
    return rand() % 41 + 1;
}

void matmulblksWithEscalar(double *a, double *b, double *c, int n, int bs, double escalar) {
    int i, j, k, offI, offJ;    /* Guess what... */
  
    for (i = 0; i < n; i += bs)
    {
        offI = i * n;
        for (j = 0; j < n; j += bs)
        {
            offJ = j * n;
            for (k = 0; k < n; k += bs)
            {
                blkmulWithEscalar(&a[offI + k], &b[offJ + k], &c[offI + j], n, bs, escalar);
            }
        }
    }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmulWithEscalar(double *ablk, double *bblk, double *cblk, int n, int bs, double escalar){
    int i, j, k, offI, offJ;    /* Guess what... again... */

    for (i = 0; i < bs; i++)
    {
        int offI = i * n;
        for (j = 0; j < bs; j++)
        {
            int offJ = i * n;
            for (k = 0; k < bs; k++)
            {
                cblk[offI + j] += ablk[offI + k] * bblk[offJ + k];
            }
            cblk[offI + j] *= escalar;
        }
    }
}

void matIntmulblks(double *a, int *b, double *c, int n, int bs){
    int i, j, k, offI, offJ;    /* Guess what... */

    for (i = 0; i < n; i += bs)
    {
        offI = i * n;
        for (j = 0; j < n; j += bs)
        {
            offJ = j * n;
            for (k = 0; k < n; k += bs)
            {
                blkmulwithIntMat(&a[offI + k], &b[offJ + k], &c[offI + j], n, bs);
            }
        }
    }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmulwithIntMat(double *ablk, int *bblk, double *cblk, int n, int bs){
    int i, j, k, offI, offJ;    /* Guess what... again... */

    for (i = 0; i < bs; i++)
    {
        offI = i * n;
        for (j = 0; j < bs; j++)
        {
            offJ = j * n;
            for (k = 0; k < bs; k++)
            {
                cblk[offI + j] += ablk[offI + k] * bblk[offJ + k];
            }
        }
    }
}

//Pruebas

void recorrerArreglo(double *matriz, int N, int orden){
    for(int i=0;i<N;i++){
        int cont = 0;
        for(int j=0;j<N;j++){
	        printf("%f ||| ",getDoubleValor(matriz,i,j,N,orden));
            cont++;
            if (cont == N){
                printf("\n");
                cont = 0;
            }
        }
    } 
}

void recorrerArreglo2(int *matriz, int N, int orden){
    for(int i=0;i<N;i++){
        int cont = 0;
        for(int j=0;j<N;j++){
	        printf("%d ||| ",getIntValor(matriz,i,j,N,orden));
            cont++;
            if (cont == N){
                printf("\n");
                cont = 0;
            }
        }
    } 
}