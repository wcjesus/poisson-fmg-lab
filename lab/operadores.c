#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "matriz.h"
#include "operadores.h"

void restringe(double **r_fina, double **r_grossa, int nx_fina, int ny_fina)
{
    int I, J;
    int nx_grossa = nx_fina / 2;
    int ny_grossa = ny_fina / 2;
    zera_matriz(r_grossa, nx_grossa, ny_grossa);
    for (I = 1; I < nx_grossa; I++)
        for (J = 1; J < ny_grossa; J++) {
            int i = 2 * I, j = 2 * J;
            r_grossa[I][J] = (
                4.0 * r_fina[i][j]
                + 2.0 * (r_fina[i-1][j] + r_fina[i+1][j] + r_fina[i][j-1] + r_fina[i][j+1])
                + (r_fina[i-1][j-1] + r_fina[i-1][j+1] + r_fina[i+1][j-1] + r_fina[i+1][j+1])
            ) / 16.0;
        }
}

void prolonga(double **e_grossa, double **e_fina, int nx_fina, int ny_fina)
{
    int i, j;
    int nx_grossa = nx_fina / 2;
    int ny_grossa = ny_fina / 2;
    zera_matriz(e_fina, nx_fina, ny_fina);

    for (i = 0; i <= nx_grossa; i++)
        for (j = 0; j <= ny_grossa; j++)
            e_fina[2*i][2*j] = e_grossa[i][j];
    for (i = 0; i < nx_grossa; i++)
        for (j = 0; j <= ny_grossa; j++)
            e_fina[2*i+1][2*j] = 0.5 * (e_grossa[i][j] + e_grossa[i+1][j]);
    for (i = 0; i <= nx_grossa; i++)
        for (j = 0; j < ny_grossa; j++)
            e_fina[2*i][2*j+1] = 0.5 * (e_grossa[i][j] + e_grossa[i][j+1]);
    for (i = 0; i < nx_grossa; i++)
        for (j = 0; j < ny_grossa; j++)
            e_fina[2*i+1][2*j+1] = 0.25 * (e_grossa[i][j] + e_grossa[i+1][j]
                                         + e_grossa[i][j+1] + e_grossa[i+1][j+1]);
}

