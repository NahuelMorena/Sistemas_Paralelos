#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include "matlib.h"

#define COORDINATOR 0

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

    int provided;
    int rank, numProcs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int strip_size = N / numProcs;

    if (rank == COORDINATOR) {
        //double st = dwalltime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Reparto parcial de las matrices A, B, C, D
    MPI_Scatter(MA, strip_size, MPI_DOUBLE, MA, strip_size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MB, strip_size, MPI_DOUBLE, MB, strip_size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MC, strip_size, MPI_DOUBLE, MB, strip_size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MD, strip_size, MPI_INT, MD, strip_size, MPI_INT, COORDINATOR, MPI_COMM_WORLD);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp single
        {
            local_min_a = local_max_a = MA[0];
            local_min_b = local_max_b = MB[0];
            local_sum_a = local_sum_b = 0.0;
        }

        /* MaxA, MinA, PromA */
        #pragma omp for reduction(max:local_max_a) reduction(min:local_min_a) reduction(+:local_sum_a) nowait
        for (i = 0; i < strip_size; i++){
            vecd_calc_min_max_sum(&MA[i*N], N, &local_min_a, &local_max_a, &local_sum_a);
        }

        /* MaxB, MinB, PromB */
        #pragma omp for reduction(max:local_max_b) reduction(min:local_min_b) reduction(+:local_sum_b) nowait
        for (i = 0; i < strip_size; i++){
            vecd_calc_min_max_sum(&MB[i*N], N, &local_min_b, &local_max_b, &local_sum_b);
        }

        /* D = Pot2(D) */
        #pragma omp for
        for (i = 0; i < strip_size; i++){
            int *v = &MD[i*N];
            veci_mult_elem(v, v, v, N);
        }

        MPI_Allreduce(&local_min_a, &min_a, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_max_a, &max_a, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_sum_a, &sum_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
        MPI_Allreduce(&local_min_b, &min_b, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_max_b, &max_b, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_sum_b, &sum_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        #pragma omp single
        {
            avg_a = sum_a / (double)S;
            avg_b = sum_b / (double)S;
            e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b);  
        }

        //Envia la matriz B completa a todos los nodos
        MPI_Bcast(MB, N, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

        /* R = e * (A x B) */
        #pragma omp for nowait
        for (i = 0; i < strip_size; i += B){
            for (j = 0; j < N; j += B){
                MR[i*N+j] = 0.0;
                for (k = 0; k < N; k += B){
                    blk_matd_mult(&MR[i*N+j], &MA[i*N+k], &MB[j*N+k], N, B);
                }
                blk_matd_mult_d(&MR[i*N+j], &MR[i*N+j], e, N, B);
            }
        }

        //Envia la matriz D completa a todos los nodos
        MPI_Bcast(MD, N, MPI_INT, COORDINATOR, MPI_COMM_WORLD);

        /* T = C x D */
        #pragma omp for
        for (i = 0; i < strip_size; i += B){
            for (j = 0; j < N; j += B){
                MT[i*N+j] = 0.0;
                for (k = 0; k < N; k += B){
                    blk_matd_mult_mati(&MT[i*N+j], &MC[i*N+k], &MD[j*N+k], N, B);
                }
            }
        }

        /* R = R + T */
        #pragma omp for
        for (i = 0; i < strip_size; i++){
            vecd_sum(&MR[i*N], &MR[i*N], &MT[i*N], N);
        }
    }
    
    MPI_Gather(MR, strip_size, MPI_DOUBLE, MR, strip_size, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);

    printf("valor del max_a %f\n", max_a);
    printf("valor del min_a %f\n", min_a);
    printf("valor del avg_a %f\n", avg_a);
    printf("valor del max_b %f\n", max_b);
    printf("valor del min_b %f\n", min_b);
    printf("valor del avg_b %f\n", avg_b);        
    printf("valor del escalar %f\n", e);

    mat_free(&m);

    MPI_Finalize();

    return 0;
}