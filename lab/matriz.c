#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "matriz.h"

double **aloca_matriz(int nx, int ny)
{
    int i;
    double **A = malloc((nx + 1) * sizeof(double *));

    if (A == NULL) 
    { 
	    printf("Erro de alocacao.\n"); 
	    return NULL; 
    }

    for (i = 0; i <= nx; i++) 
    {
        A[i] = malloc((ny + 1) * sizeof(double));
        if (A[i] == NULL) 
	{
            int k;
            printf("Erro de alocacao.\n");
            for (k = 0; k < i; k++) free(A[k]);
            free(A);
            return NULL;
        }
    }
    return A;
}

void libera_matriz(double **A, int nx)
{
    int i;
    if (A == NULL) return;
    for (i = 0; i <= nx; i++) free(A[i]);
    free(A);
}

void zera_matriz(double **A, int nx, int ny)
{
    int i, j;
    for (i = 0; i <= nx; i++)
        for (j = 0; j <= ny; j++)
            A[i][j] = 0.0;
}



