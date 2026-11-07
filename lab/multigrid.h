#ifndef LAB_MULTIGRID_H
#define LAB_MULTIGRID_H
#include "resultados.h"
int Multigrid(int n, double custom_tol, double *out_h, double *out_erro_inf,
              ResultEntry *accum);
#endif
