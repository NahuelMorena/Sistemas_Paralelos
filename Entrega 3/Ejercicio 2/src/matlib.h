#ifndef _MATLIB_H
#define _MATLIB_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

/*
 * Obtiene el tiempo en segundos desde un momento constante
 */

static double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

/*
 * Lectura de matrices por fila o por columna desde archivo
 */

static void matReadRow(FILE *f, void *mat, int size, int n) 
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            void *p = mat + ((i * n) + j) * size;
            fread(p, size, 1, f);
        }
    }
}

static void matReadCol(FILE *f, void *mat, int size, int n) 
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            void *p = mat + ((j * n) + i) * size;
            fread(p, size, 1, f);
        }
    }
}

/*
 * Inicialización y liberación de matrices
 */
typedef struct
{
    /* Parámetro */
    int N; /* Tamaño de las matrices */
    int S; /* Cantidad de elementos de las matrices */
    int B; /* Tamaño de las submatrices */
    int T; /* Número de hilos */
    /* Matrices */
    double *MA, *MB, *MC, *MR, *MT;
    int *MD;
    /* Archivos */
    FILE *ifile, *ofile;
} mats_t;

static void arg_parse(int argc, char *argv[], int open_files, mats_t *m)
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s matrix_size block_size [num_threads = 1] [input_file] [output_file]\n", argv[0]);
        exit(1);
    }

    m->N = atoi(argv[1]);
    m->S = m->N * m->N;
    m->B = atoi(argv[2]);
    
    if (argc >= 4)
        m->T = atoi(argv[3]);
    else
        m->T = 1;
    
    if (open_files && argc >= 5)
    {
        m->ifile = fopen(argv[4], "rb");
        if (!m->ifile)
        {
            fprintf(stderr, "no se pudo abrir el archivo de entrada '%s'\n", argv[4]);
            exit(1);
        }
    }
    else
        m->ifile = NULL;

    if (open_files && argc >= 6)
    {
        m->ofile = fopen(argv[5], "wb");
        if (!m->ofile)
        {
            fprintf(stderr, "no se pudo abrir el archivo de salida '%s'\n", argv[5]);
            exit(1);
        }
    }
    else
        m->ofile = NULL;
}

static void mat_init(int argc, char *argv[], int coord, int procs, mats_t *m)
{
    arg_parse(argc, argv, coord, m);

    int s = coord ? m->S : m->S / procs;

    m->MA = (double *)malloc(sizeof(double) * s);
    m->MB = (double *)malloc(sizeof(double) * m->S);
    m->MC = (double *)malloc(sizeof(double) * s);
    m->MR = (double *)malloc(sizeof(double) * s);
    m->MT = (double *)malloc(sizeof(double) * s);
    m->MD = (int *)malloc(sizeof(int) * m->S);

    if (coord) {
        if (m->ifile)
        {
            matReadRow(m->ifile, m->MA, sizeof(double), m->N);
            matReadCol(m->ifile, m->MB, sizeof(double), m->N);
            matReadRow(m->ifile, m->MC, sizeof(double), m->N);
            matReadCol(m->ifile, m->MD, sizeof(int), m->N);
            fclose(m->ifile);
        }
        else
        {
            srand((int)(dwalltime() * 1000.0));
            for (int i = 0; i < m->S; i++)
            {
                m->MA[i] = 1.0;
                m->MB[i] = 1.0;
                m->MC[i] = 1.0;
                m->MD[i] = rand() % 41 + 1;
            }
        }
    }
}

static void mat_free(const mats_t *m)
{
    if (m->ofile)
    {
        fwrite(m->MR, sizeof(double), m->S, m->ofile);
        fclose(m->ofile);
    }

    free(m->MA);
    free(m->MB);
    free(m->MC);
    free(m->MR);
    free(m->MT);
    free(m->MD);
}

/*
 * Operación de matrices y submatrices de diferentes tipos
 */

static inline void vecd_sum(double *r, double *a, double *b, int size)
{
    int i;
    for (i = 0; i < size; i++)
        r[i] = a[i] + b[i];
}

static inline void vecd_calc_min_max_sum(double *m, int size, double *min, double *max, double *sum)
{
    int i;
    for (i = 0; i < size; i++)
    {
        double v = m[i];

        if (v < *min)
            *min = v;
        
        if (v > *max)
            *max = v;
        
        *sum += v;
    }
}

static inline void veci_mult_elem(int *r, int *a, int *b, int size)
{
    int i;
    for (i = 0; i < size; i++)
        r[i] = a[i] * b[i];
}

static inline void blk_matd_mult_d(double *r, double *a, double b, int size, int blksize)
{
    int i, j;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            r[i*size+j] = a[i*size+j] * b;
        }
    }
}

static inline void blk_matd_mult_mati(double *r, double *a, int *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            for (k = 0; k < blksize; k++)
            {
                r[i*size+j] += a[i*size+k] * b[j*size+k];
            }
        }
    }
}

static inline void blk_matd_mult(double *r, double *a, double *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < blksize; i++)
    {
        for (j = 0; j < blksize; j++)
        {
            for (k = 0; k < blksize; k++)
            {
                r[i*size+j] += a[i*size+k] * b[j*size+k];
            }
        }
    }
}
/*
static inline void matd_mult_mati(double *r, double *a, int *b, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < size; i += blksize)
    {
        for (j = 0; j < size; j += blksize)
        {
            r[i*size+j] = 0.0;
            for (k = 0; k < size; k += blksize)
            {
                blk_matd_mult_mati(&r[i*size+j], &a[i*size+k], &b[j*size+k], size, blksize);
            }
        }
    }
}
*/
/*
static inline void matd_mult_with_scalar(double *r, double *a, double *b, double c, int size, int blksize)
{
    int i, j, k;
    for (i = 0; i < size; i += blksize)
    {
        for (j = 0; j < size; j += blksize)
        {
            r[i*size+j] = 0.0;
            for (k = 0; k < size; k += blksize)
            {
                blk_matd_mult(&r[i*size+j], &a[i*size+k], &b[j*size+k], size, blksize);
            }
            blk_matd_mult_d(&r[i*size+j], &r[i*size+j], c, size, blksize);
        }
    }
}
*/
#endif