#include <stdio.h>
#include <omp.h>
#include "../matlib.h"

int main(int argc, char *argv[])
{
    mats_t m;
    mat_init(argc, argv, &m);

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

    double min_a, max_a, sum_a, avg_a;
    double min_b, max_b, sum_b, avg_b;
    double e;

    int i, j, k;
    
    double st = dwalltime();

    omp_set_num_threads(T);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp single
        {
            min_a = max_a = MA[0];
            min_b = max_b = MB[0];
            sum_a = sum_b = 0.0;
        }

        /* MaxA, MinA, PromA */
        #pragma omp for reduction(max:max_a) reduction(min:min_a) reduction(+:sum_a) nowait
        for (i = 0; i < N; i++)
        {
            vecd_calc_min_max_sum(&MA[i*N], N, &min_a, &max_a, &sum_a);
        }

        /* MaxB, MinB, PromB */
        #pragma omp for reduction(max:max_b) reduction(min:min_b) reduction(+:sum_b) nowait
        for (i = 0; i < N; i++)
        {
            vecd_calc_min_max_sum(&MB[i*N], N, &min_b, &max_b, &sum_b);
        }

        /* D = Pot2(D) */
        #pragma omp for
        for (i = 0; i < N; i++)
        {
            int *v = &MD[i*N];
            veci_mult_elem(v, v, v, N);
        }

        #pragma omp single
        {
            avg_a = sum_a / (double)S;
            avg_b = sum_b / (double)S;
            e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b);
        }

        /* R = e * (A x B) */
        #pragma omp for nowait
        for (i = 0; i < N; i += B)
        {
            for (j = 0; j < N; j += B)
            {
                MR[i*N+j] = 0.0;
                for (k = 0; k < N; k += B)
                {
                    blk_matd_mult(&MR[i*N+j], &MA[i*N+k], &MB[j*N+k], N, B);
                }
                blk_matd_mult_d(&MR[i*N+j], &MR[i*N+j], e, N, B);
            }
        }

        /* T = C x D */
        #pragma omp for
        for (i = 0; i < N; i += B)
        {
            for (j = 0; j < N; j += B)
            {
                MT[i*N+j] = 0.0;
                for (k = 0; k < N; k += B)
                {
                    blk_matd_mult_mati(&MT[i*N+j], &MC[i*N+k], &MD[j*N+k], N, B);
                }
            }
        }

        /* R = R + T */
        #pragma omp for
        for (i = 0; i < N; i++)
        {
            vecd_sum(&MR[i*N], &MR[i*N], &MT[i*N], N);
        }
    }

    double t = dwalltime() - st;
    printf("%.04f\n", t);

    mat_free(&m);

    return 0;
}