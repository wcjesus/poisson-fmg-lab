#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "direto.h"

// FATORACOES UNICAS DA MATRIZ DA MALHA MAIS GROSSA 
static void direto_monta_A(double **A, int nx, int ny)
{
    int nj = ny - 1;
    int n  = (nx - 1) * (ny - 1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = 0.0;
    for (int i = 1; i < nx; i++)
        for (int j = 1; j < ny; j++) {
            int k = (i - 1) * nj + (j - 1);
            A[k][k] = 4.0;
            if (i > 1)      A[k][(i - 2) * nj + (j - 1)] = -1.0;
            if (i < nx - 1) A[k][i * nj + (j - 1)]       = -1.0;
            if (j > 1)      A[k][(i - 1) * nj + (j - 2)] = -1.0;
            if (j < ny - 1) A[k][(i - 1) * nj + j]       = -1.0;
        }
}

static void direto_monta_b(double **u, double **f, int nx, int ny, double h, double *b)
{
    int nj = ny - 1;
    double h2 = h * h;
    for (int i = 1; i < nx; i++)
        for (int j = 1; j < ny; j++) {
            int k = (i - 1) * nj + (j - 1);
            double bi = h2 * f[i][j];
            if (i == 1)      bi += u[0][j];
            if (i == nx - 1) bi += u[nx][j];
            if (j == 1)      bi += u[i][0];
            if (j == ny - 1) bi += u[i][ny];
            b[k] = bi;
        }
}

static void direto_escreve(double **u, double *x, int nx, int ny)
{
    int nj = ny - 1;
    for (int i = 1; i < nx; i++)
        for (int j = 1; j < ny; j++)
            u[i][j] = x[(i - 1) * nj + (j - 1)];
}

static double **LU_L = NULL;
static double **LU_U = NULL;
static int      LU_n = 0;
static double **LDL_L = NULL;
static double  *LDL_D = NULL;
static int      LDL_n = 0;
static double **GJ_Ainv = NULL;
static int      GJ_n    = 0;

static void LU_fatora(int nx, int ny)
{
    int n = (nx - 1) * (ny - 1);
    double **A = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) A[i] = malloc(n * sizeof(double));
    direto_monta_A(A, nx, ny);
    LU_L = malloc(n * sizeof(double *));
    LU_U = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        LU_L[i] = calloc(n, sizeof(double));
        LU_U[i] = calloc(n, sizeof(double));
    }
    for (int k = 0; k < n; k++) {
        LU_L[k][k] = 1.0;
        for (int j = k; j < n; j++) {
            double soma = 0.0;
            for (int s = 0; s < k; s++) soma += LU_L[k][s] * LU_U[s][j];
            LU_U[k][j] = A[k][j] - soma;
        }
        for (int i = k + 1; i < n; i++) {
            double soma = 0.0;
            for (int s = 0; s < k; s++) soma += LU_L[i][s] * LU_U[s][k];
            LU_L[i][k] = (A[i][k] - soma) / LU_U[k][k];
        }
    }
    LU_n = n;
    for (int i = 0; i < n; i++) free(A[i]);
    free(A);
}

static void LU_resolve(double *b, double *x, int n)
{
    double *z = malloc(n * sizeof(double));
    z[0] = b[0];
    for (int i = 1; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < i; j++) soma += LU_L[i][j] * z[j];
        z[i] = b[i] - soma;
    }
    x[n - 1] = z[n - 1] / LU_U[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        double soma = 0.0;
        for (int j = i + 1; j < n; j++) soma += LU_U[i][j] * x[j];
        x[i] = (z[i] - soma) / LU_U[i][i];
    }
    free(z);
}

static void LDL_fatora(int nx, int ny)
{
    int n = (nx - 1) * (ny - 1);
    double **A = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) A[i] = malloc(n * sizeof(double));
    direto_monta_A(A, nx, ny);
    LDL_L = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) LDL_L[i] = calloc(n, sizeof(double));
    LDL_D = malloc(n * sizeof(double));
    for (int j = 0; j < n; j++) {
        LDL_L[j][j] = 1.0;
        double soma_d = 0.0;
        for (int v = 0; v < j; v++) soma_d += LDL_D[v] * LDL_L[j][v] * LDL_L[j][v];
        LDL_D[j] = A[j][j] - soma_d;
        for (int i = j + 1; i < n; i++) {
            double soma_l = 0.0;
            for (int v = 0; v < j; v++) soma_l += LDL_L[i][v] * LDL_D[v] * LDL_L[j][v];
            LDL_L[i][j] = (A[i][j] - soma_l) / LDL_D[j];
        }
    }
    LDL_n = n;
    for (int i = 0; i < n; i++) free(A[i]);
    free(A);
}

static void LDL_resolve(double *b, double *x, int n)
{
    double *z = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < i; j++) soma += LDL_L[i][j] * z[j];
        z[i] = b[i] - soma;
    }
    for (int i = 0; i < n; i++) y[i] = z[i] / LDL_D[i];
    for (int i = n - 1; i >= 0; i--) {
        double soma = 0.0;
        for (int j = i + 1; j < n; j++) soma += LDL_L[j][i] * x[j];
        x[i] = y[i] - soma;
    }
    free(z); free(y);
}

static void GJ_fatora(int nx, int ny)
{
    int n = (nx - 1) * (ny - 1);
    double **M = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) 
    {
        M[i] = calloc(2 * n, sizeof(double));
        M[i][n + i] = 1.0;
    }

    direto_monta_A(M, nx, ny);
    int    *l = malloc(n * sizeof(int));
    double *s = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) 
    {
        l[i] = i; s[i] = 0.0;
        for (int j = 0; j < n; j++)
            if (fabs(M[i][j]) > s[i]) s[i] = fabs(M[i][j]);
    }

    for (int k = 0; k < n; k++) 
    {
        double rmax = 0.0;
        int p = k;
        for (int i = k; i < n; i++) 
	{
            double r = fabs(M[l[i]][k] / s[l[i]]);
            if (r > rmax) { rmax = r; p = i; }
        }

        int temp = l[p]; l[p] = l[k]; l[k] = temp;
        double piv = M[l[k]][k];

        for (int j = 0; j < 2 * n; j++) M[l[k]][j] /= piv;
        for (int i = 0; i < n; i++) 
	{
            if (i != k) 
	    {
                double mult = M[l[i]][k];
                for (int j = 0; j < 2 * n; j++) M[l[i]][j] -= mult * M[l[k]][j];
            }
        }
    }

    GJ_Ainv = malloc(n * sizeof(double *));

    for (int i = 0; i < n; i++) 
    {
        GJ_Ainv[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) GJ_Ainv[i][j] = M[l[i]][n + j];
    }

    GJ_n = n;

    for (int i = 0; i < n; i++) free(M[i]);
    free(M); free(l); free(s);
}

static void GJ_resolve(double *b, double *x, int n)
{
    for (int i = 0; i < n; i++) 
    {
        double soma = 0.0;
        for (int j = 0; j < n; j++) soma += GJ_Ainv[i][j] * b[j];
        x[i] = soma;
    }
}

void direto_fatora_todos(int nx, int ny)
{
    LU_fatora (nx, ny);
    LDL_fatora(nx, ny);
    GJ_fatora (nx, ny);
}

void direto_libera_fatores(void)
{
    if (LU_L) 
    { 
	    for (int i = 0; i < LU_n; i++) free(LU_L[i]); 
	    free(LU_L); LU_L = NULL; 
    }

    if (LU_U) 
    { 
	    for (int i = 0; i < LU_n; i++) free(LU_U[i]); 
	    free(LU_U); 
	    LU_U = NULL; 
    }

    LU_n = 0;

    if (LDL_L) 
    { 
	    for (int i = 0; i < LDL_n; i++) free(LDL_L[i]); 
	    free(LDL_L); 
	    LDL_L = NULL; 
    }

    free(LDL_D); 
    LDL_D = NULL; 
    LDL_n = 0;

    if (GJ_Ainv) 
    { 
	    for (int i = 0; i < GJ_n; i++) free(GJ_Ainv[i]); 
	    free(GJ_Ainv); 
	    GJ_Ainv = NULL; 
    }

    GJ_n = 0;
}
//(B) METODOS DIRETOS (pos-fatorados) 

void LU(double **u, double **f, int nx, int ny, double h)
{
    int n = (nx - 1) * (ny - 1);

    /* Fatora somente se os fatores ainda não existirem. */
    if (LU_L == NULL || LU_n != n)
    {
        if (LU_L != NULL || LU_U != NULL)
            direto_libera_fatores();

        LU_fatora(nx, ny);
    }

    double *b = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));

    if (!b || !x)
    {
        free(b);
        free(x);
        fprintf(stderr, "Erro de alocacao em LU.\n");
        exit(EXIT_FAILURE);
    }

    direto_monta_b(u, f, nx, ny, h, b);
    LU_resolve(b, x, n);
    direto_escreve(u, x, nx, ny);

    free(b);
    free(x);
}

void LDL(double **u, double **f, int nx, int ny, double h)
{
    int n = (nx - 1) * (ny - 1);

    if (LDL_L == NULL || LDL_D == NULL || LDL_n != n)
    {
        if (LDL_L != NULL || LDL_D != NULL)
            direto_libera_fatores();

        LDL_fatora(nx, ny);
    }

    double *b = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));

    if (!b || !x)
    {
        free(b);
        free(x);
        fprintf(stderr, "Erro de alocacao em LDL.\n");
        exit(EXIT_FAILURE);
    }

    direto_monta_b(u, f, nx, ny, h, b);
    LDL_resolve(b, x, n);
    direto_escreve(u, x, nx, ny);

    free(b);
    free(x);
}


void gaussjordan(double **u, double **f, int nx, int ny, double h)
{
    int n = (nx - 1) * (ny - 1);

    if (GJ_Ainv == NULL || GJ_n != n)
    {
        if (GJ_Ainv != NULL)
            direto_libera_fatores();

        GJ_fatora(nx, ny);
    }

    double *b = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));

    if (!b || !x)
    {
        free(b);
        free(x);
        fprintf(stderr, "Erro de alocacao em Gauss-Jordan.\n");
        exit(EXIT_FAILURE);
    }

    direto_monta_b(u, f, nx, ny, h, b);
    GJ_resolve(b, x, n);
    direto_escreve(u, x, nx, ny);

    free(b);
    free(x);
}

