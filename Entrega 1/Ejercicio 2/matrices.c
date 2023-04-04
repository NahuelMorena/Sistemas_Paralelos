#include <stdio.h>
#include <stdlib.h>
//Utilizado para calcular las potencias
#include <math.h>
//Utilizado para contar tiempo empleado
#include <sys/time.h>
#include <time.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

//Declaración de funciones 
int getValor(double *matriz,int fila,int columna,int N, int orden);
int getRandomNumber();
double dwalltime();

void recorrerArreglo(double *matriz, int N, int orden);
void recorrerArreglo2(int *matriz, int N, int orden);

int main(int argc, char*argv[]){
    //Validar parametros
    if (argc < 2){
        fprintf(stderr, "usage %s matrix_size\n",argv[0]);
        exit(1);
    }

    //Tomar parametro n para el tamaño de las matrices NxN
    int N = atoi(argv[1]);
    printf("N = %i\n",N);

    
    //Asignar matrices (Se utiliza alternativa 4 mostrada en la teoria)
    double *A,*B,*C,*R,*CxD;
    int *D;

    //indices
    int i,j,k;

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
    printf("valor del escalar: %f\n",escalar);

    //Multiplicación AxB
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            sum = 0;
            for (k=0; k<N; k++){
                sum += A[i*N+k] * B[k+j*N]; 
            }
            R[i+j*N] = sum;
        }
    }
    printf("MATRIZ AxB\n");
    recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Multiplicación de AxB por escalar
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            R[i+j*N] *= escalar;
        }
    }

    printf("MATRIZ (AxB) * escalar\n");
    recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Pot2(D)

    printf("MATRIZ D\n");
    recorrerArreglo2(D,N,ORDENXCOLUMNAS);
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            D[i+j*N] = (int) pow(D[i+j*N],2);
        }
    } 

    printf("MATRIZ D DESPUES DE POTENCIA DE 2\n");
    recorrerArreglo2(D,N,ORDENXCOLUMNAS);


    //Multiplicación CxD
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            sum = 0;
            for (k=0; k<N; k++){
                sum += C[i*N+k] * D[k+j*N]; 
            }
            CxD[i+j*N] = sum;
        }
    }

    printf("MATRIZ CxD\n");
    recorrerArreglo(CxD,N,ORDENXCOLUMNAS);

    //Suma entre escalar*AxB + CxD
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            R[i+j*N] += CxD[i+j*N];
        }
    }

    printf("MATRIZ (escalar *AxB) + CxD\n");
    recorrerArreglo(R,N,ORDENXCOLUMNAS);

    //Finaliza conteo de tiempo
    endtime = dwalltime() - timetick;
    printf("Tiempo en segundos %f\n", endtime);

    //Obtención del resultado


    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(CxD);

    printf("finish!\n");

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