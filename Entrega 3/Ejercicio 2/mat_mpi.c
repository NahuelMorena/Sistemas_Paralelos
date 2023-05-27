#include <stdio.h>
#include <mpi.h>
#include "matlib.h"

int main(int argc, char *argv[]){
    mats_t m;
    mat_init(argc, argv, &m);

    int N = m.N;
    int S = m.S;
    int B = m.B;
    int T = m.T;
    int P = m.P;

    double *MA = m.MA;
    double *MB = m.MB;
    double *MC = m.MC;
    double *MT = m.MT;
    double *MR = m.MR;
    int *MD = m.MD;

    double min_a, max_a, sum_a, avg_a;
    double local_min_a, local_max_a, local_sum_a;
    double min_b, max_b, sum_b, avg_b;
    double local_min_b, local_max_b, local_sum_b;
    double e;

    int i, j, k;

    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int strip_size = N / numProcs;
    //int cellAmount = local_size * N;
    int COORDINATOR = 0;

    if (rank == COORDINATOR) {
        //double st = dwalltime();
    }

    // Reparto de las matrices A y B
    MPI_Scatter(MA, strip_size, 
                MPI_DOUBLE, MA, strip_size,
                MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    MPI_Scatter(MB, strip_size,
                MPI_DOUBLE, MB, strip_size,
                MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    /* MaxA, MinA, AvgA */
    local_max_a = local_min_a = MA[0];
    local_sum_a = 0.0;
    for (i = 0; i < strip_size; i++){
        vecd_calc_min_max_sum(&MA[i*N], N, &local_min_a, &local_max_a, &local_sum_a);
    }

    /* MaxB, MinB, AvgB */
    local_max_b = local_min_b = MB[0];
    local_sum_b = 0.0;
    for (i = 0; i < strip_size; i++){
        vecd_calc_min_max_sum(&MB[i*N], N, &local_min_b, &local_max_b, &local_sum_b);
    }


    MPI_Allreduce(&local_min_a, &min_a, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_a, &max_a, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_sum_a, &sum_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
   
    MPI_Allreduce(&local_min_b, &min_b, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_b, &max_b, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_sum_b, &sum_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == COORDINATOR){
        avg_a = sum_a / (double)S;
        avg_b = sum_b / (double)S;
        e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b);
    }

    
    printf("valor del max_a %f\n", max_a);
    printf("valor del min_a %f\n", min_a);
    printf("valor del avg_a %f\n", avg_a);
    printf("valor del max_b %f\n", max_b);
    printf("valor del min_b %f\n", min_b);
    printf("valor del avg_b %f\n", avg_b);
    printf("valor del escalar %f\n", e);
    
    MPI_Finalize();
    mat_free(&m);
}