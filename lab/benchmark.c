#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"
#include "resultados.h"
#include "benchmark.h"

double bench_t_pre      = 0.0;
long   bench_n_pre      = 0;
double bench_t_pos      = 0.0;
long   bench_n_pos      = 0;
double bench_t_B        = 0.0;
long   bench_n_B_stages = 0;
long   bench_n_B_iter   = 0;
double bench_t_C        = 0.0;
long   bench_n_C_iter   = 0;
double bench_t_guess    = 0.0;
double bench_t_vcycle   = 0.0;
long   bench_n_vcycle   = 0;

double ts_diff(struct timespec a, struct timespec b)
{
    return (b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) * 1e-9;
}

void bench_reset(void)
{
    bench_t_pre = 0; bench_n_pre = 0;
    bench_t_pos = 0; bench_n_pos = 0;
    bench_t_B = 0; bench_n_B_stages = 0; bench_n_B_iter = 0;
    bench_t_C = 0; bench_n_C_iter = 0;
    bench_t_guess = 0;
    bench_t_vcycle = 0; bench_n_vcycle = 0;
}

void bench_accumulate_into(ResultEntry *r)
{
    if (!r) return;
    r->n_vcycle   += bench_n_vcycle;
    r->t_vcycle   += bench_t_vcycle;
    r->t_pre      += bench_t_pre;
    r->n_pre      += bench_n_pre;
    r->t_pos      += bench_t_pos;
    r->n_pos      += bench_n_pos;
    r->t_B        += bench_t_B;
    r->n_B_stages += bench_n_B_stages;
    r->n_B_iter   += bench_n_B_iter;
    r->t_C        += bench_t_C;
    r->n_C_iter   += bench_n_C_iter;
    r->t_guess    += bench_t_guess;
}


void bench_print(void)
{
    printf("\n--- BENCHMARK ---\n");
    printf("Tempo pra gerar o chute: %.6f s", bench_t_guess);

    if (cfg_C == 1 || cfg_C == 2) 
    {
        double med = bench_n_C_iter > 0 ? bench_t_C / bench_n_C_iter : 0.0;
        printf("  (media por iteracao de C: %.9f s, %ld iteracoes)", med, bench_n_C_iter);
    }

    printf("\n");
    printf("Tempo gasto com v-cycles (loop principal): %.6f s\n", bench_t_vcycle);
    printf("Pre-smoothing : total %.6f s  (%ld chamadas)  media por chamada %.9f s\n",
           bench_t_pre, bench_n_pre, bench_n_pre > 0 ? bench_t_pre / bench_n_pre : 0.0);
    printf("Pos-smoothing : total %.6f s  (%ld chamadas)  media por chamada %.9f s\n",
           bench_t_pos, bench_n_pos, bench_n_pos > 0 ? bench_t_pos / bench_n_pos : 0.0);
    printf("B (coarse): total %.6f s  (%ld etapas)  media por etapa %.9f s\n",
           bench_t_B, bench_n_B_stages,
           bench_n_B_stages > 0 ? bench_t_B / bench_n_B_stages : 0.0);

    if (cfg_B >= 4) 
    {
        printf("B iterativo: %ld iteracoes totais  media por iteracao %.9f s\n",
               bench_n_B_iter, bench_n_B_iter > 0 ? bench_t_B / bench_n_B_iter : 0.0);
    }

    printf("V-cycle: %ld ciclos  media por ciclo %.9f s\n",
           bench_n_vcycle, bench_n_vcycle > 0 ? bench_t_vcycle / bench_n_vcycle : 0.0);
}

