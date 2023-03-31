#include <stdio.h>
#include <stdlib.h>

//int N = 10;

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
    double *A,*B,*C,*R;
    int *D;

    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    R=(double*)malloc(sizeof(double)*N*N);

    D=(int*)malloc(sizeof(int)*N*N);


    //Inicializar matrices

    //Operar matrices

    //Obtención del resultado


    printf("finish!\n");
    
}