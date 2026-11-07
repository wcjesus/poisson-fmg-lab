#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "configuracao.h"
#include "resultados.h"
#include "benchmark.h"
#include "multigrid.h"
#include "direto.h"

#ifndef VIDEO
static ResultEntry g_results[MAX_RESULTS];
static int g_nresults = 0;

static void imprime_overview(void)
{
    printf("=================================================================================================================\n"); 
    printf("|                                               OVERVIEW FINAL                                                  |\n");
    printf("=================================================================================================================\n");
    printf("|      A      |      B      |      C      |    p   |  Vcyc  | t/cic(ms)| tPre(ms) | tPos(ms) | Ncoar  |   t(s)  |\n");
    printf("=================================================================================================================\n");

    for (int k = 0; k < g_nresults; k++) {
        ResultEntry *r = &g_results[k];
        double t_cic_ms = r->n_vcycle > 0 ? 1000.0 * r->t_vcycle / r->n_vcycle : 0.0;
        double t_pre_ms = r->n_pre    > 0 ? 1000.0 * r->t_pre    / r->n_pre    : 0.0;
        double t_pos_ms = r->n_pos    > 0 ? 1000.0 * r->t_pos    / r->n_pos    : 0.0;
        long   n_coarse = (r->b >= 4) ? r->n_B_iter : r->n_B_stages;

        printf("|%-13s|%-13s|%-13s|%-8.4f|%-8ld|%-10.4f|%-10.4f|%-10.4f|%-8ld|%-9.3f|\n",
               nome_A(r->a), nome_B(r->b), nome_C(r->c),
               r->p,
               r->n_vcycle, t_cic_ms, t_pre_ms, t_pos_ms,
               n_coarse, r->t_total);
    printf("-----------------------------------------------------------------------------------------------------------------\n");
    }
    
}

#endif

static void executa_e_calcula_p(ResultEntry *out,int numero)
{
    double h1, h2, e1, e2, p;
    struct timespec t0, t1;

#ifdef VIDEO
    // O video acompanha uma malha, sem repetir o FMG para estimar p. 
    Multigrid(numero, cfg_tol_global, NULL, NULL, out);//funcao em multigrid.c 
    return;
#endif

    clock_gettime(CLOCK_MONOTONIC, &t0);
    if (Multigrid(numero,     cfg_tol_global, &h1, &e1, out) != 0) return;//malha, lida da linha de comando; tolerancia globa; passo(será modificado); erro1(vai ser calculado); resultentry 
    if (numero > INT_MAX/2 || Multigrid(2 * numero, cfg_tol_global, &h2, &e2, out) != 0) return;//execucao com malha com dobro de pontos pra calcular p
    clock_gettime(CLOCK_MONOTONIC, &t1);

    p = log(e1 / e2) / log(h1 / h2);//ordem de convergência

    printf("\n ORDEM DE CONVERGENCIA (norma do maximo) \n");
    printf("A = %s\n", nome_A(cfg_A));
    printf("B = %s\n", nome_B(cfg_B));
    printf("C = %s\n", nome_C(cfg_C));
    printf("N1=%d  h1=%.6e  e1=%.6e\n", cfg_N,   h1, e1);
    printf("N2=%d  h2=%.6e  e2=%.6e\n", 2*cfg_N, h2, e2);
    printf("p  = %.6f\n\n", p);

    if (out) {
        out->a = cfg_A;
        out->b = cfg_B;
        out->c = cfg_C;
        out->e1 = e1;
        out->e2 = e2;
        out->p  = p;
        out->t_total = ts_diff(t0, t1);
    }
}



int main(int argc, char** argv)
{
#ifndef VIDEO
    int escolher=0;
#endif
    //teste para validar se os argumentos foram corretamente passados
    if (argc < 2) {
        fprintf(stderr, "Uso: %s N [modo]\nN = 16, 32, 64, 128, ...\n", argv[0]);
        return 1;
    }

    printf("\nLABORATORIO MULTIGRID\n");
#ifdef VIDEO
    printf("Modo video: escolha A, B e C para acompanhar uma unica malha.\n");
#else
    if( (argc>2 && atoi(argv[2]) == 1) || argc == 2)
    {
    	printf("Escolher configuracao? (1=sim / 0=nao): \n");
    	if (scanf("%d", &escolher) != 1) return 1;
    }
#endif
    // MODO ESGOTAR COMBINACOES
#ifndef VIDEO
    if (!escolher || (argc >2 && atoi(argv[2]) == 0))
    {
        printf("\nModo esgotar: todas as combinacoes de A, B, C com macros padrao.\n");
        g_nresults = 0;
        
	//loop que alterna entre todas as combinações
        for (int a = 1; a <= 2; a++)//métodos A(suavizadores)
        	for (int b = 1; b <= 5; b++)//métodos B(solucionadores)
       			 for (int c = 1; c <= 5; c++)//métodos C(inicializadores)
        			{
					
                            printf("%d/50\n", (a - 1) * 25 + (b - 1) * 5 + c);//imprime quantas configurações ja fez
            				cfg_A = a; cfg_B = b; cfg_C = c;
            				cfg_K       = K;//macro, definido em configuracao.h. Representa a frequência da solução
            				cfg_omega_A = (a == 2) ? OMEGAGS : OMEGA; //A sendo GS, usa o peso do Gauss Seidel em vez do de Jacobi
            				cfg_omega_B = (b == 5) ? OMEGAGS : OMEGA; //se o método for Jacobi, ai sim usa o peso dele. Macros em configuracao.h
            				cfg_omega_C = (c == 2) ? OMEGAGS : OMEGA;
            				cfg_NU1     = NU1;//macro, número de pre-suavizações. definido em configuraco.h
            				cfg_NU2     = NU2;// pos suavizacoes
            				cfg_tol_global = T;
            				cfg_tol_B   = TOLERANCIA_B;
            				cfg_tol_C   = TOLERANCIA_C;
            				cfg_N       = NUM;
            				printf("\n\n##################################################################\n");
            				printf("### COMBINACAO: A=%s  B=%s  C=%s%*s###\n",
       					nome_A(a),nome_B(b),nome_C(c),
       					63 - (int)strlen("### COMBINACAO: A=")
           				   - (int)strlen(nome_A(a))
           				   - (int)strlen("  B=")
         				   - (int)strlen(nome_B(b))
        				   - (int)strlen("  C=")
       					   - (int)strlen(nome_C(c)),
       					"");
            				printf("##################################################################\n");

            				ResultEntry *slot = NULL;
            				if (g_nresults < MAX_RESULTS) 
					{
                				slot = &g_results[g_nresults++];
                				memset(slot, 0, sizeof(*slot));
            				}
            				executa_e_calcula_p(slot,atoi(argv[1])); //funcao de convergencia. em 
					
        			}

        imprime_overview();
direto_libera_fatores();
return 0;

    }

#endif
    // MODO INTERATIVO 

    int usar_padrao;//variavel logica, usar padrao ou nao.
    printf("\nDeseja usar padrao? (1=sim / 0=nao): ");
    if (scanf("%d", &usar_padrao) != 1) return 1;

    // ---- A ---- 
    printf("\n A: suavizador (iterativo) \n");
    printf("  1 - Jacobi\n");
    printf("  2 - Gauss-Seidel\n");
    printf("Escolha A: ");

    //se A nao for Jacobi ou Gauss Seidel...
    if (scanf("%d", &cfg_A) != 1 || cfg_A < 1 || cfg_A > 2) {
        printf("Opcao invalida.\n"); return 1;
    }
    //se escolheu usar padrão. Escolhe o peso de A entre Jacobi e GS
    if (usar_padrao)
        cfg_omega_A = (cfg_A == 2) ? OMEGAGS : OMEGA;

    //se nao vai usar padrao, pega da entrada
    if (!usar_padrao) {
        printf("Peso A (omega_A): ");
        scanf("%lf", &cfg_omega_A);
        printf("numero de pre suavizacoes (NU1): ");
        scanf("%d", &cfg_NU1);
        printf("numero de pos suavizacoes (NU2): ");
        scanf("%d", &cfg_NU2);
    }

    // ---- B ---- 
    printf("\n B: solver da malha mais grossa \n");
    printf("  1 - LU            (direto)\n");
    printf("  2 - LDL           (direto)\n");
    printf("  3 - Gauss-Jordan  (direto)\n");
    printf("  4 - Jacobi        (iterativo)\n");
    printf("  5 - Gauss-Seidel  (iterativo)\n");
    printf("Escolha B: ");
    if (scanf("%d", &cfg_B) != 1 || cfg_B < 1 || cfg_B > 5) {
        printf("Opcao invalida.\n"); return 1;
    }

    int B_iter = (cfg_B == 4 || cfg_B == 5); //B_iter começa como variavel logica. se B for iterativo (opcoes 4 e 5) ele comeca como 1.
    if (usar_padrao && B_iter) //B é iterativo e vai usar padrao? escolhe peso entre jacobi e seidel
        cfg_omega_B = (cfg_B == 5) ? OMEGAGS : OMEGA;

    if (!usar_padrao && B_iter) {// B iterativo mas nao vai usar padrao? le do usuario
        printf("Peso B (omega_B): ");
        scanf("%lf", &cfg_omega_B);
        printf("tolerancia de B: ");
        scanf("%lf", &cfg_tol_B);
    }

    // ---- C ---- 
    printf("\n C: aproximacao inicial do FMG \n");
    printf("  1 - Jacobi        (iterativo)\n");
    printf("  2 - Gauss-Seidel  (iterativo)\n");
    printf("  3 - LU            (direto)\n");
    printf("  4 - LDL           (direto)\n");
    printf("  5 - Gauss-Jordan  (direto)\n");
    printf("Escolha C: ");
    if (scanf("%d", &cfg_C) != 1 || cfg_C < 1 || cfg_C > 5) {
        printf("Opcao invalida.\n"); return 1;
    }

    int C_iter = (cfg_C == 1 || cfg_C == 2);//var. logica
    if (usar_padrao && C_iter) 
        cfg_omega_C = (cfg_C == 2) ? OMEGAGS : OMEGA;

    if (!usar_padrao && C_iter) {
        printf("Peso C (omega_C): ");
        scanf("%lf", &cfg_omega_C);
        printf("tolerancia de C: ");
        scanf("%lf", &cfg_tol_C);
    }

    //  Parametros globais 
    if (!usar_padrao) {
        printf("\nfrequencia K: ");//macros em configuracao.h
        scanf("%lf", &cfg_K);
        printf("tolerancia global: ");
        scanf("%lf", &cfg_tol_global);
    }

    printf("\nCONFIGURACAO EM USO\n");//mostra os parâmetros que está usando
    printf("A = %s  (omega_A = %g, NU1=%d, NU2=%d)\n",
           nome_A(cfg_A), cfg_omega_A, cfg_NU1, cfg_NU2);
    printf("B = %s", nome_B(cfg_B));
    if (B_iter) printf("  (omega_B = %g, tol_B = %.3e)", cfg_omega_B, cfg_tol_B);
    else        printf("  (direto)");
    printf("\n");
    printf("C = %s", nome_C(cfg_C));
    if (C_iter) printf("  (omega_C = %g, tol_C = %.3e)", cfg_omega_C, cfg_tol_C);
    else        printf("  (direto)");
    printf("\n");
    printf("K = %g  tolerância global = %.3e  N = %d\n", cfg_K, cfg_tol_global, cfg_N);
    printf("MÁXIMO DE ITERAÇÃO GROSSO = %d  MÁXIMO DE CICLOS = %d\n", MAX_ITER_GROSSO, MAX_CICLOS);

    executa_e_calcula_p(NULL, atoi(argv[1]));
    return 0;
}

