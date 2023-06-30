#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include "matlib.h"

#define COORDINATOR 0
#define DEBUG 0

int main(int argc, char *argv[]) {
    int rank, numProcs, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    mats_t m;
    mat_init(argc, argv, rank == COORDINATOR, numProcs, &m);

    int N = m.N;
    int S = m.S;
    int B = m.B;
    int T = m.T;

    double *MA = m.MA;
    double *MB = m.MB;
    double *MC = m.MC;
    double *MT = m.MT;
    double *MR = m.MR;
    int *MD = m.MD;

    double lmin[2], min[2];
    double lmax[2], max[2];
    double lsum[2], sum[2];
    double avg[2];
    double e;

    int i, j, k;

    int strip_size = N / numProcs;
    int so = rank * strip_size;

    omp_set_num_threads(T);

    //Times
    double st, t, comst, comt = 0.0, comtotal = 0.0;

    if (rank == COORDINATOR) {
        st = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Envio las matrices A, B, C, D
    comst = MPI_Wtime();
    MPI_Bcast(MB, S, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MA, strip_size * N, MPI_DOUBLE, MA, strip_size * N, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MC, strip_size * N, MPI_DOUBLE, MC, strip_size * N, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Scatter(MD, strip_size * N, MPI_INT, MD, strip_size * N, MPI_INT, COORDINATOR, MPI_COMM_WORLD);
    comt += MPI_Wtime() - comst;

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp single
        {
            for (i = 0; i < 2; i++) {
                lmin[i] = lmax[i] = MA[i];
                lsum[i] = 0.0;
            }
        }

        // MaxA, MinA, AvgA (MA: Scatter)
        #pragma omp for reduction(max:lmax[0]) reduction(min:lmin[0]) reduction(+:lsum[0]) nowait
        for (i = 0; i < strip_size; i++) {
            vecd_calc_min_max_sum(&MA[i*N], N, &lmin[0], &lmax[0], &lsum[0]);
        }

        // MaxB, MinB, AvgB (MB: Bcast)
        #pragma omp for reduction(max:lmax[1]) reduction(min:lmin[1]) reduction(+:lsum[1]) nowait
        for (i = 0; i < strip_size; i++) {
            vecd_calc_min_max_sum(&MB[(i+so)*N], N, &lmin[1], &lmax[1], &lsum[1]);
        }

        /* D = Pot2(D) */
        #pragma omp for
        for (i = 0; i < strip_size; i++) {
            int *v = &MD[i*N];
            veci_mult_elem(v, v, v, N);
        }

        #pragma omp single
        {
            if (rank == COORDINATOR) {
                for (i = 0; i < 2; i++) {
                    min[i] = lmin[i];
                    max[i] = lmax[i];
                    sum[i] = 0.0;
                }
            }

            comst = MPI_Wtime();
            MPI_Reduce(lmin, min, 2, MPI_DOUBLE, MPI_MIN, COORDINATOR, MPI_COMM_WORLD);
            MPI_Reduce(lmax, max, 2, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);
            MPI_Reduce(lsum, sum, 2, MPI_DOUBLE, MPI_SUM, COORDINATOR, MPI_COMM_WORLD);
            comt += MPI_Wtime() - comst;

            if (rank == COORDINATOR) {
                for (i = 0; i < 2; i++) {
                    avg[i] = sum[i] / (double)S;
                }
                e = (max[0] * max[1] - min[0] * min[1]) / (avg[0] * avg[1]);
            }

            comst = MPI_Wtime();
            MPI_Bcast(&e, 1, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
            MPI_Allgather(MD, strip_size * N, MPI_INT, MD, strip_size * N, MPI_INT, MPI_COMM_WORLD);
            comt += MPI_Wtime() - comst;
        }

        // R = e * (A x B)
        #pragma omp for
        for (i = 0; i < strip_size; i += B) {
            for (j = 0; j < N; j += B) {
                MR[i*N+j] = 0.0;
                for (k = 0; k < N; k += B) {
                    blk_matd_mult(&MR[i*N+j], &MA[i*N+k], &MB[j*N+k], N, B);
                }
                blk_matd_mult_d(&MR[i*N+j], &MR[i*N+j], e, N, B);
            }
        }

        // T = C x D
        #pragma omp for
        for (i = 0; i < strip_size; i += B) {
            for (j = 0; j < N; j += B) {
                MT[i*N+j] = 0.0;
                for (k = 0; k < N; k += B) {
                    blk_matd_mult_mati(&MT[i*N+j], &MC[i*N+k], &MD[j*N+k], N, B);
                }
            }
        }

        // R = R + T
        #pragma omp for
        for (i = 0; i < strip_size; i++) {
            vecd_sum(&MR[i*N], &MR[i*N], &MT[i*N], N);
        }
    }

    comst = MPI_Wtime();
    MPI_Gather(MR, strip_size * N, MPI_DOUBLE, MR, strip_size * N, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    comt += MPI_Wtime() - comst;

    if (rank == COORDINATOR) {
        t = MPI_Wtime() - st;
        printf("%.04f\n", t);
    }

    MPI_Reduce(&comt, &comtotal, 1, MPI_DOUBLE, MPI_SUM, COORDINATOR, MPI_COMM_WORLD);

    if (rank == COORDINATOR) {
        printf("Tiempo de com.: %.09f\n", comtotal);
    }

    mat_free(&m);

    MPI_Finalize();

    return 0;
}