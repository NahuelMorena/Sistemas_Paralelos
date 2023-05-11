#include <stdio.h>
#include <stdlib.h>
//Utilizado para contar tiempo empleado
#include <sys/time.h>
//Utilizado para obtener numeros random
#include <time.h>

#define DEBUG 0

static double start_time, total_time = 0.0;

//Declaración de funciones 
double dwalltime();
void timelog_start();
void timelog_total();
void timelog(const char *desc);

static void matReadRow(FILE *f, void *mat, int size, int n) 
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            void *p = mat + ((i * n) + j) * size;
            fread(p, size, 1, f);
        }
    }
}

static void matReadCol(FILE *f, void *mat, int size, int n) 
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            void *p = mat + ((j * n) + i) * size;
            fread(p, size, 1, f);
        }
    }
}

static void blkMultDD(double *r, double *a, double *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            for (k = 0; k < blksize; k++)
            {
                r[i*size+j] += a[i*size+k] * b[j*size+k];
            }
        }
    }
}

static void blkMultDS(double *r, double *a, double s, int size, int blksize)
{
    int i, j;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            r[i*size+j] = a[i*size+j] * s;
        }
    }
}

static void matMultDDS(double *r, double *a, double *b, double s, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < size; i += blksize)
    {
        for (j = 0; j < size; j += blksize)
        {
            r[i*size+j] = 0.0;
            for (k = 0; k < size; k += blksize)
            {
                blkMultDD(&r[i*size+j], &a[i*size+k], &b[j*size+k], size, blksize);
            }
            blkMultDS(&r[i*size+j], &r[i*size+j], s, size, blksize);
        }
    }
}

static void blkMultDI(double *r, double *a, int *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            for (k = 0; k < blksize; k++)
            {
                r[i*size+j] += a[i*size+k] * b[j*size+k];
            }
        }
    }
}

static void matMultDI(double *r, double *a, int *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < size; i += blksize)
    {
        for (j = 0; j < size; j += blksize)
        {
            r[i*size+j] = 0.0;
            for (k = 0; k < size; k += blksize)
            {
                blkMultDI(&r[i*size+j], &a[i*size+k], &b[j*size+k], size, blksize);
            }
        }
    }
}

int main(int argc, char*argv[]){
    //Validar parametros
    if (argc < 3){
        fprintf(stderr, "usage %s matrix_size block_size [input_file output_file]\n",argv[0]);
        exit(1);
    }

    FILE *ifile, *ofile;
    if (argc >= 5) {
        ifile = fopen(argv[3], "rb");
        ofile = fopen(argv[4], "wb");
        if (!ifile)
            fprintf(stderr, "no se pudo abrir el archivo de entrada '%s'\n", argv[3]);
        if (!ofile)
            fprintf(stderr, "no se pudo abrir el archivo de salida '%s'\n", argv[3]);
    } else {
        ifile = ofile = NULL;
    }

    //Tomar parametro n para el tamaño de las matrices NxN
    int N = atoi(argv[1]);
    #if DEBUG == 1
    printf("N = %i\n",N);
    #endif

    int bs = atoi(argv[2]);
    #if DEBUG == 1
    printf("bs = %i\n",bs);
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

    srand(time(NULL));

    if (ifile && ofile) {
        matReadRow(ifile, A, sizeof(double), N);
        matReadCol(ifile, B, sizeof(double), N);
        matReadRow(ifile, C, sizeof(double), N);
        matReadCol(ifile, D, sizeof(int), N);
    } else {
        for (i = 0; i < size; i++) {
            A[i] = 1.0;
            B[i] = 1.0;
            C[i] = 1.0;
            D[i] = rand() % 41 + 1;
        }
    }

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

    //Multiplicación AxB
    timelog_start();
    matMultDDS(R, A, B, escalar, N, bs);
    timelog("R = (AxB) * e");

    //Pot2(D)  
    timelog_start();
    for(int j=0;j<N;j++) {
        offJ = j*N;
        for(int i=0;i<N;i++) {
            int v = D[i+offJ];
            D[i+offJ] = v * v;
        }
    }
    timelog("D = Pot2(D)");

    //Multiplicación CxD
    timelog_start();
    matMultDI(CxD, C, D, N, bs);
    timelog("CxD = CxD");

    //Suma entre (escalar*AxB) + CxD
    timelog_start();
    for (i=0; i<N; i++) {
        offI = i * N;
        for (j=0; j<N; j++) {
            R[offI+j] += CxD[offI+j];
        }
    }
    timelog("R = R + CxD");

    timelog_total();

    if (ofile)
        fwrite(R, sizeof(double), size, ofile);

    if (ifile)
        fclose(ifile);

    if (ofile)
        fclose(ofile);

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

/* Funciones para depurar el tiempo de ejecución */
double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void timelog_start() {
    start_time = dwalltime();
}

void timelog_total() {
    #if DEBUG == 1
    printf("Tiempo total: %.02fms\n\n\n", total_time * 1000.0);
    #else
    printf("%f\n", total_time);
    #endif
}

void timelog(const char *desc) {
    double time = dwalltime() - start_time;
    total_time += time;

    #if DEBUG == 1
    if (desc)
        printf("Tiempo '%s': %.02fms\n", desc, time * 1000.0);
    #endif
}