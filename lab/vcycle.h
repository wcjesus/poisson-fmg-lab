#ifndef LAB_VCYCLE_H
#define LAB_VCYCLE_H
void v_cycle(double **u, double **f, int nx, int ny, double h,
             double ***f_niveis, double ***r_niveis,
             double ***e_niveis, double ***u_novo_niveis, int nivel_atual);
#endif
