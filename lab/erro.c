#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "problema.h"
#include "erro.h"

double norma_erro(double **u, int nx, int ny, double hx, double hy)
{
    int i, j;
    double soma = 0.0, erro;
    for (i = 1; i < nx; i++)
        for (j = 1; j < ny; j++) 
	{
            erro = u[i][j] - uex(i * hx, j * hy);
            soma += erro * erro;
        }
    return sqrt(hx * hy * soma);
}

double norma_infinito_erro(double **u, int nx, int ny, double hx, double hy)
{
    int i, j;
    double max_erro = 0.0, erro;
    for (i = 1; i < nx; i++)
        for (j = 1; j < ny; j++) 
	{
            erro = fabs(u[i][j] - uex(i * hx, j * hy));
            if (erro > max_erro) max_erro = erro;
        }
    return max_erro;
}

