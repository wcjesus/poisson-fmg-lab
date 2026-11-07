#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "matriz.h"
#include "residuo.h"

void calcula_residuo(double **u, double **f, double **r, int nx, int ny, double h)
{
    int i, j;
    zera_matriz(r, nx, ny);
    for (i = 1; i < nx; i++)
        for (j = 1; j < ny; j++) {
            double Au = (4.0 * u[i][j] - u[i-1][j] - u[i+1][j]
                                     - u[i][j-1] - u[i][j+1]) / (h * h);
            r[i][j] = f[i][j] - Au;
        }
}

double norma_residuo(double **r, int nx, int ny, double hx, double hy)
{
    int i, j;
    double soma = 0.0;
    for (i = 1; i < nx; i++)
        for (j = 1; j < ny; j++)
            soma += r[i][j] * r[i][j];
    return sqrt(hx * hy * soma);
}

