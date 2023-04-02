#include <stdio.h>
#include <stdlib.h>
//Utilizado para calcular las potencias
#include <math.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

//Declaración de funciones 
void setDoubleValor(double *matriz, int fila, int columna, int orden, int N, double valor);
void setIntValor(int *matriz, int fila, int columna, int orden, int N, int valor);
int getValor(double *matriz,int fila,int columna,int N, int orden);
int getRandomNumber();
void getDataFromArray(double *matriz, int N, int orden, double *min, double *max, double *prom);
void applyPow(int pow, int *matriz, int N);

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
    double *A,*B,*C,*R,*AxB,*CxD;
    int *D;

    //indices
    int i,j,k;
    int sum;

    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    AxB=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    R=(double*)malloc(sizeof(double)*N*N);

    D=(int*)malloc(sizeof(int)*N*N);

    CxD=(double*)malloc(sizeof(double)*N*N);

    double maxA, minA, promA = 0.0;
    double maxB, minB, promB = 0.0;

    //Inicializar matrices
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
	        setDoubleValor(A,i,j,ORDENXFILAS,N,1.0);
	        setDoubleValor(B,i,j,ORDENXCOLUMNAS,N,1.0);
            setDoubleValor(C,i,j,ORDENXFILAS,N,1.0);
            setIntValor(D,i,j,ORDENXCOLUMNAS,N,getRandomNumber());
        }
    }   

    //Empieza a contar el tiempo

    //Operar matrices

    getDataFromArray(A, N, ORDENXFILAS, &maxA, &minA, &promA);
    getDataFromArray(B, N, ORDENXCOLUMNAS, &maxB, &minB, &promB);
    //applyPow(2,D,N);

    printf("MATRIZ A\n");
    recorrerArreglo(A,N,ORDENXFILAS);
    printf("maxA: ");
    printf("%f",maxA);
    printf("\n");
    printf("minA: ");
    printf("%f",minA);
    printf("\n");
    printf("promA: ");
    printf("%f",promA);
    printf("\n");
    printf("MATRIZ B\n");
    recorrerArreglo(B,N,ORDENXFILAS);
    printf("maxB: ");
    printf("%f",maxB);
    printf("\n");
    printf("minB: ");
    printf("%f",minB);
    printf("\n");
    printf("promB: ");
    printf("%f",promB);
    printf("\n");
    printf("MATRIZ C\n");
    recorrerArreglo(C,N,ORDENXFILAS);

    printf("MATRIZ D\n");
    recorrerArreglo2(D,N,ORDENXCOLUMNAS);
    applyPow(2,D,N);
    printf("MATRIZ D DESPUES DE POTENCIA DE 2\n");
    recorrerArreglo2(D,N,ORDENXCOLUMNAS);

    double escalar = (maxA * maxB - minA * minB)/(promA * promB);
    printf("valor del escalar: %f\n",escalar);

    //Multiplicación AxB
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            sum = 0;
            for (k=0; k<N; k++){
                sum += A[i*N+k] * B[k+j*N]; 
            }
            AxB[i+j*N] = sum;
        }
    }
    printf("MATRIZ AxB\n");
    recorrerArreglo(AxB,N,ORDENXCOLUMNAS);

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

    //Suma entre AxB + CxD
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            AxB[i+j*N] += CxD[i+j*N];
        }
    }

    printf("MATRIZ AxB + CxD\n");
    recorrerArreglo(AxB,N,ORDENXCOLUMNAS);

    //Multiplicación de matriz por escalar
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            AxB[i+j*N] *= escalar;
        }
    }

    printf("MATRIZ (AxB + CxD) * escalar\n");
    recorrerArreglo(AxB,N,ORDENXCOLUMNAS);

    //Finaliza conteo de tiempo

    //Obtención del resultado


    //Liberando memoria
    free(A);
    free(B);
    free(C);
    free(R);
    free(D);
    free(AxB);
    free(CxD);

    printf("finish!\n");

    return(0);
}

//--------------------------------------------------------------------------------------------------

void setDoubleValor(double *matriz, int fila, int columna, int orden, int N, double valor){
    if(orden==ORDENXFILAS){
        matriz[fila*N+columna]=valor;
    }else{
        matriz[fila+columna*N]=valor;
    }
}

void setIntValor(int *matriz, int fila, int columna, int orden, int N, int valor){
    if(orden==ORDENXFILAS){
        matriz[fila*N+columna]=valor;
    }else{
        matriz[fila+columna*N]=valor;
    }
}

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

void applyPow(int p, int *matriz, int N){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
	        setIntValor(matriz,i,j,ORDENXCOLUMNAS,N,(int)pow(getIntValor(matriz,i,j,N,ORDENXCOLUMNAS),p));
        }
    } 
}

void getDataFromArray(double *matriz, int N, int orden, double *min, double *max, double *prom){
    *min = *max = getDoubleValor(matriz,0,0,N,orden);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
	        double valor = getDoubleValor(matriz,i,j,N,orden);
            if (valor < *min){
                *min = valor;
            }

            if (valor > *max){
                *max = valor;
            }

            *prom += valor;
        }
    } 
    printf("valor de prom: %f\n",*prom);
    *prom = *prom/(N*N);
    printf("valor de prom despues de dividir N*N: %f\n",*prom);
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