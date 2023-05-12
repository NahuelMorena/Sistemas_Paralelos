#include <stdio.h>
#include "matlib.h"

int main(int argc, char *argv[])
{
    mats_t m;
    mat_init(argc, argv, &m);

    int N = m.N;
    int S = m.S;
    int B = m.B;

    double *MA = m.MA;
    double *MB = m.MB;
    double *MC = m.MC;
    double *MT = m.MT;
    double *MR = m.MR;
    int *MD = m.MD;

    double min_a, max_a, sum_a, avg_a;
    double min_b, max_b, sum_b, avg_b;
    double e;
    
    double st = dwalltime();
    
    /* MaxA, MinA, PromA */
    min_a = max_a = MA[0];
    sum_a = 0.0;
    vecd_calc_min_max_sum(MA, S, &min_a, &max_a, &sum_a);
    avg_a = sum_a / (double)S;

    /* MaxB, MinB, PromB */
    min_b = max_b = MB[0];
    sum_b = 0.0;
    vecd_calc_min_max_sum(MB, S, &min_b, &max_b, &sum_b);
    avg_b = sum_b / (double)S;

    e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b);

    /* R = e * (A x B) */
    matd_mult_with_scalar(MR, MA, MB, e, N, B);

    /* D = Pot2(D) */
    veci_mult_elem(MD, MD, MD, S);

    /* T = C x D */
    matd_mult_mati(MT, MC, MD, N, B);

    /* R = R + T */
    vecd_sum(MR, MR, MT, S);

    double t = dwalltime() - st;
    printf("%.04f\n", t);

    mat_free(&m);

    return 0;
}