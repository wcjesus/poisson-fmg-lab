#ifndef LAB_RESULTADOS_H
#define LAB_RESULTADOS_H
typedef struct {
    int a, b, c;
    double e1, e2, p;
    double t_total;
    long n_vcycle;
    double t_vcycle;
    double t_pre;   long n_pre;
    double t_pos;   long n_pos;
    double t_B;     long n_B_stages; long n_B_iter;
    double t_C;     long n_C_iter;
    double t_guess;
} ResultEntry;
#define MAX_RESULTS 100
#endif
