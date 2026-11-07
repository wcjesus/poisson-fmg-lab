#ifndef LAB_DIRETO_H
#define LAB_DIRETO_H
void direto_fatora_todos(int nx, int ny);
void direto_libera_fatores(void);
void LU(double **u, double **f, int nx, int ny, double h);
void LDL(double **u, double **f, int nx, int ny, double h);
void gaussjordan(double **u, double **f, int nx, int ny, double h);
#endif
