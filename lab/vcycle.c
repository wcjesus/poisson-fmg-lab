#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"
#include "benchmark.h"
#include "matriz.h"
#include "operadores.h"
#include "residuo.h"
#include "suavizacao.h"
#include "direto.h"
#include "vcycle.h"

void v_cycle(double **u, double **f, int nx, int ny, double h,
             double ***f_niveis, double ***r_niveis,
             double ***e_niveis, double ***u_novo_niveis, int nivel_atual)
{
    int i, j;
    double **u_original = u;
    double **u_novo = u_novo_niveis[nivel_atual];

    if (nx <= 16 || ny <= 16)
    {
        /*
	 *  (B) MALHA GROSSA 
         * - B direto:    executa 1 vez.
         * - B iterativo: executa ate max_grosso ou ate res < tol_B.
         */
        //B 
        int n_iter_max = (cfg_B <= 3) ? 1 : MAX_ITER_GROSSO;
        double **r_B = r_niveis[nivel_atual];

        struct timespec tb0, tb1;
        clock_gettime(CLOCK_MONOTONIC, &tb0);
        int iters_B = 0;

        for (int iter_B = 0; iter_B < n_iter_max; iter_B++)
        {
		//se o metódo for direto, resolve diretamente e sai do loop:
            if (cfg_B == 1) { LU         (u, f, nx, ny, h); iters_B = 1; break; }
            if (cfg_B == 2) { LDL        (u, f, nx, ny, h); iters_B = 1; break; }
            if (cfg_B == 3) { gaussjordan(u, f, nx, ny, h); iters_B = 1; break; }

            //se for iterativo:

            double **res = smooth_step(u, u_novo, f, nx, ny, h, cfg_B-3, cfg_omega_B);
            if (res == u_novo)
	    {
		    double **t = u;
		    u = u_novo;
		    u_novo = t;
	    }

            iters_B++;
            calcula_residuo(u, f, r_B, nx, ny, h);
            if (norma_residuo(r_B, nx, ny, h, h) < cfg_tol_B) break;
        }

        clock_gettime(CLOCK_MONOTONIC, &tb1);
        bench_t_B += ts_diff(tb0, tb1);
        bench_n_B_stages++;
        if (cfg_B >= 4) bench_n_B_iter += iters_B;
        
	if (u!=u_original)
	{
		for (int i=0;i<=nx;i++)
			for (int j=0; j<=ny; j++)
				u_original[i][j]=u[i][j];
	}


        return;
    }

    double **r = r_niveis[nivel_atual];
    int nx_grossa = nx / 2;
    int ny_grossa = ny / 2;
    double **r_grossa = f_niveis[nivel_atual + 1];
    double **e_grossa = e_niveis[nivel_atual + 1];
    

    //Pre-smoothing 
for (i = 0; i < cfg_NU1; i++)
{
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    double **res = smooth_step(u, u_novo, f, nx, ny, h, cfg_A, cfg_omega_A);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    bench_t_pre += ts_diff(ts0, ts1);
    bench_n_pre++;

    if (res == u_novo) { double **t = u; u = u_novo; u_novo = t; }
}

    calcula_residuo(u, f, r, nx, ny, h);
    restringe(r, r_grossa, nx, ny);

    zera_matriz(e_grossa, nx_grossa, ny_grossa);
    v_cycle(e_grossa, r_grossa, nx_grossa, ny_grossa, 2.0 * h,
            f_niveis, r_niveis, e_niveis, u_novo_niveis, nivel_atual + 1);

    prolonga(e_grossa, u_novo, nx, ny);
    for (i = 1; i < nx; i++)
        for (j = 1; j < ny; j++)
            u[i][j] += u_novo[i][j];

    // Pos-smoothing 
    for (i = 0; i < cfg_NU2; i++)
{
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    double **res = smooth_step(u, u_novo, f, nx, ny, h, cfg_A, cfg_omega_A);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    bench_t_pos += ts_diff(ts0, ts1);
    bench_n_pos++;

    if (res == u_novo) { double **t = u; u = u_novo; u_novo = t; }
}

if(u!=u_original)
{
	for (int i=0;i<=nx;i++)
		for (int j=0;j<=ny;j++)
			u_original[i][j]=u[i][j];
}

}

