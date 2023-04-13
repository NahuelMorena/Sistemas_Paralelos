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
void matmulblks(double *a, double *b, double *c, int n, int bs);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);

void matmulblks2(double *a, int *b, double *c, int n, int bs);
void blkmul2(double *ablk, int *bblk, double *cblk, int n, int bs);
double dwalltime();

void recorrerArreglo(double *matriz, int N, int orden);
void recorrerArreglo2(int *matriz, int N, int orden);

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
    int i,j,k;
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
    minA = maxA = A[0];
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
	        double valor = A[i*N+j];
            if (valor < minA){
                minA = valor;
            }

            if (valor > maxA){
                maxA = valor;
            }

            promA += valor;
        }
    } 
    promA = promA/(size);

    //obteniendo minB, maxB y promB
    minB = maxB = B[0];
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
	        double valor = B[i+j*N];
            if (valor < minB){
                minB = valor;
            }

            if (valor > maxB){
                maxB = valor;
            }
            promB += valor;
        }
    } 
    promB = promB/(size);

    double escalar = (maxA * maxB - minA * minB)/(promA * promB);
    //printf("valor del escalar: %f\n",escalar);

    //Multiplicación AxB

    matmulblks(A, B, R, N, bs);

    //for (i=0; i<N; i++){
    //    fila = i*N;
    //    for (j=0; j<N; j++){
    //        columna = j*N;
    //        sum = 0;
    //        for (k=0; k<N; k++){
    //            sum += A[fila+k] * B[k+columna]; 
    //        }
    //        R[i+columna] = sum * escalar;
    //    }
    //}
    //printf("MATRIZ AxB\n");
    //recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Multiplicación de AxB por escalar
    //for (i=0; i<N; i++){
    //    for (j=0; j<N; j++){
    //        R[i+j*N] *= escalar;
    //    }
    //}

    //printf("MATRIZ (AxB) * escalar\n");
    //recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Pot2(D)

    //printf("MATRIZ D\n");
    //recorrerArreglo2(D,N,ORDENXCOLUMNAS);
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            D[i+j*N] = (int) pow(D[i+j*N],2);
        }
    } 

    //printf("MATRIZ D DESPUES DE POTENCIA DE 2\n");
    //recorrerArreglo2(D,N,ORDENXCOLUMNAS);


    //Multiplicación CxD

    matmulblks2(C, D, CxD, N, bs);

    //for (i=0; i<N; i++){
    //    fila = i*N;
    //    for (j=0; j<N; j++){
    //        columna = j*N;
    //        sum = 0;
    //        for (k=0; k<N; k++){
    //            sum += C[fila+k] * D[k+columna]; 
    //        }
    //        CxD[i+columna] = sum;
    //    }
    //}

    //printf("MATRIZ CxD\n");
    //recorrerArreglo(CxD,N,ORDENXCOLUMNAS);

    //Suma entre escalar*AxB + CxD
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            R[i+j*N] += CxD[i+j*N];
        }
    }

    //printf("MATRIZ (escalar *AxB) + CxD\n");
    //recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Finaliza conteo de tiempo
    endtime = dwalltime() - timetick;
    //printf("Tiempo en segundos %f\n", endtime);

    //Obtención del resultado
    char file[10000];
    bzero(file,10000);
    char temp_str[2000];

    sprintf(temp_str, "Analisis de tiempo obtenidos en matrices de %i. \n",N);
    strcat(file, temp_str);
    sprintf(temp_str, "Tiempo obtenido de %f segundos\n", endtime);
    strcat(file, temp_str);
    sprintf(temp_str, "-----------------------------------------------\n\n");
    strcat(file, temp_str);

    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(CxD);

    printf("%s", file);
    return(0);
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

void matmulblks(double *a, double *b, double *c, int n, int bs){
  int i, j, k;    /* Guess what... */

  /* Init matrix c, just in case */  
  //initvalmat(c, n, 0.0, 0);
  
  for (i = 0; i < n; i += bs)
  {
    for (j = 0; j < n; j += bs)
    {
      for  (k = 0; k < n; k += bs)
      {
        blkmul(&a[i*n + k], &b[j*n + k], &c[i*n + j], n, bs);
      }
    }
  }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs){
  int i, j, k;    /* Guess what... again... */

  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      for  (k = 0; k < bs; k++)
      {
        cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
      }
    }
  }
}

void matmulblks2(double *a, int *b, double *c, int n, int bs){
  int i, j, k;    /* Guess what... */

  /* Init matrix c, just in case */  
  //initvalmat(c, n, 0.0, 0);
  
  for (i = 0; i < n; i += bs)
  {
    for (j = 0; j < n; j += bs)
    {
      for  (k = 0; k < n; k += bs)
      {
        blkmul2(&a[i*n + k], &b[j*n + k], &c[i*n + j], n, bs);
      }
    }
  }
}

/*****************************************************************/

/* Multiply (block)submatrices */
void blkmul2(double *ablk, int *bblk, double *cblk, int n, int bs){
  int i, j, k;    /* Guess what... again... */

  for (i = 0; i < bs; i++)
  {
    for (j = 0; j < bs; j++)
    {
      for  (k = 0; k < bs; k++)
      {
        cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
      }
    }
  }
}

//Para calcular tiempo
double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
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