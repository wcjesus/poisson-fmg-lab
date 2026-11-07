#ifndef LAB_RESIDUO_H
#define LAB_RESIDUO_H
void calcula_residuo(double **u, double **f, double **r, int nx, int ny, double h);
double norma_residuo(double **r, int nx, int ny, double hx, double hy);
#endif
