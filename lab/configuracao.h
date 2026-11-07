#ifndef LAB_CONFIGURACAO_H
#define LAB_CONFIGURACAO_H
// Macros originais 
#define K 1
#define EXP 1 
#define OMEGA (4.0/5.0)   // peso padrao do Jacobi 
#define OMEGAGS 1.0        // peso padrao do Gauss-Seidel 
#define NU1 1 //pre suavizacoes
#define NU2 1 //pos suavizacoes
#define MAX_ITER_GROSSO 100000   // max iteracoes no nivel grosso (B ou C iterativos) 
#define MAX_CICLOS 100000        // max ciclos de V-cycle no loop global 
#define TOLERANCIA_C 1e-2       // tol default do bloco C 
#define TOLERANCIA_B 1e-10       // tol default do bloco B iterativo 
#define T 1e-7                   // tol global default 
#define NUM 32

extern int cfg_A, cfg_B, cfg_C;
extern double cfg_K, cfg_omega_A, cfg_omega_B, cfg_omega_C;
extern int cfg_NU1, cfg_NU2;
extern double cfg_tol_global, cfg_tol_B, cfg_tol_C;
extern int cfg_N;
const char *nome_A(int a);
const char *nome_B(int b);
const char *nome_C(int c);
#endif
