#include <stdio.h>
#include <pthread.h>
#include "../matlib.h"

int N, S, B, T;
double *MA, *MB, *MC, *MT, *MR;
int *MD;

double min_a, max_a, sum_a, avg_a;
double min_b, max_b, sum_b, avg_b;
double e;

int bs;

pthread_mutex_t avg_mutex;
pthread_barrier_t barrier;

static void *calculate(void *arg)
{
    int id = *((int *)arg);
    
    int sb = bs * id;
    int eb = bs * (id + 1);

    /* MaxA, MinA, PromA */
    double local_min_a, local_max_a, local_sum_a;
    local_min_a = local_max_a = MA[0];
    local_sum_a = 0.0;
    for (int i = sb; i < eb; i++)
    {
        vecd_calc_min_max_sum(&MA[i*N], N, &local_min_a, &local_max_a, &local_sum_a);
    }

    /* MaxB, MinB, PromB */
    double local_min_b, local_max_b, local_sum_b;
    local_min_b = local_max_b = MB[0];
    local_sum_b = 0.0;
    for (int i = sb; i < eb; i++)
    {
        vecd_calc_min_max_sum(&MB[i*N], N, &local_min_b, &local_max_b, &local_sum_b);
    }

    /* D = Pot2(D) */
    for (int i = sb; i < eb; i++)
    {
        int *v = &MD[i*N];
        veci_mult_elem(v, v, v, N);
    }

    pthread_mutex_lock(&avg_mutex);
    if (local_min_a < min_a)
        min_a = local_min_a;
    
    if (local_max_a > max_a)
        max_a = local_max_a;
    
    sum_a += local_sum_a;

    if (local_min_b < min_b)
        min_b = local_min_b;
    
    if (local_max_b > max_b)
        max_b = local_max_b;
    
    sum_b += local_sum_b;
    pthread_mutex_unlock(&avg_mutex);

    pthread_barrier_wait(&barrier);

    if (id == 0)
    {
        avg_a = sum_a / (double)S;
        avg_b = sum_b / (double)S;
        e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b);
    }

    pthread_barrier_wait(&barrier);

    /* R = e * (A x B) */
    for (int i = sb; i < eb; i += B)
    {
        for (int j = 0; j < N; j += B)
        {
            MR[i*N+j] = 0.0;
            for (int k = 0; k < N; k += B)
            {
                blk_matd_mult(&MR[i*N+j], &MA[i*N+k], &MB[j*N+k], N, B);
            }
            blk_matd_mult_d(&MR[i*N+j], &MR[i*N+j], e, N, B);
        }
    }

    /* T = C x D */
    for (int i = sb; i < eb; i += B)
    {
        for (int j = 0; j < N; j += B)
        {
            MT[i*N+j] = 0.0;
            for (int k = 0; k < N; k += B)
            {
                blk_matd_mult_mati(&MT[i*N+j], &MC[i*N+k], &MD[j*N+k], N, B);
            }
        }
    }

    pthread_barrier_wait(&barrier);

    /* R = R + T */
    for (int i = sb; i < eb; i++)
    {
        vecd_sum(&MR[i*N], &MR[i*N], &MT[i*N], N);
    }

    pthread_exit(0);
}

int main(int argc, char *argv[])
{
    mats_t m;
    mat_init(argc, argv, &m);

    N = m.N;
    S = m.S;
    B = m.B;
    T = m.T;

    MA = m.MA;
    MB = m.MB;
    MC = m.MC;
    MT = m.MT;
    MR = m.MR;
    MD = m.MD;

    bs = N / T;

    min_a = max_a = MA[0];
    sum_a = 0.0;
    min_b = max_b = MB[0];
    sum_b = 0.0;

    pthread_t threads[T];
    int thread_ids[T];

    pthread_mutex_init(&avg_mutex, NULL);
    pthread_barrier_init(&barrier, NULL, T);
    
    double st = dwalltime();

    // Crear hilos
    for (int i = 0; i < T; i++)
    {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, calculate, &thread_ids[i]);
    }

    // Espera a que todos los hilos terminen
    for (int i = 0; i < T; i++)
    {
        pthread_join(threads[i], NULL);
    }

    double t = dwalltime() - st;
    printf("%.04f\n", t);

    pthread_mutex_destroy(&avg_mutex);
    pthread_barrier_destroy(&barrier);

    mat_free(&m);

    return 0;
}