#ifndef LAB_BENCHMARK_H
#define LAB_BENCHMARK_H
#include <time.h>
#include "resultados.h"
extern double bench_t_pre, bench_t_pos, bench_t_B, bench_t_C, bench_t_guess, bench_t_vcycle;
extern long bench_n_pre, bench_n_pos, bench_n_B_stages, bench_n_B_iter, bench_n_C_iter, bench_n_vcycle;
double ts_diff(struct timespec a, struct timespec b);
void bench_reset(void);
void bench_accumulate_into(ResultEntry *r);
void bench_print(void);
#endif
