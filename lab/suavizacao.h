#ifndef LAB_SUAVIZACAO_H
#define LAB_SUAVIZACAO_H
double gaussseidel(double **u, double **f, int i, int j, double h, double omega);
double jacobi_ponderado(double **u, double **f, int i, int j, double h, double omega);
double **smooth_step(double **u, double **u_novo, double **f,
                     int nx, int ny, double h, int method, double omega);
void suaviza(double **u, double **u_novo, double **f, int nx, int ny, double h);
void suavizacoes(double **u, double **u_novo, double **f,
                 int nx, int ny, double h, int numero);
#endif
