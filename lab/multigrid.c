#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"
#include "resultados.h"
#include "benchmark.h"
#include "matriz.h"
#include "erro.h"
#include "problema.h"
#include "operadores.h"
#include "residuo.h"
#include "suavizacao.h"
#include "vcycle.h"
#include "direto.h"
#include "multigrid.h"

#ifdef VIDEO
// Captura somente amostras de v; o renderizador 3D fica no video.sh. 
static int video_numero_quadro = 0;

static void video_salva_quadro(double **v, int n, const char *etapa,
                              int ciclo, double residuo)
{
    char caminho[128];
    video_numero_quadro++;
    snprintf(caminho, sizeof(caminho), "video_dados/quadro_%06d.dat", video_numero_quadro);
    FILE *dados = fopen(caminho, "w");
    if (!dados) { perror(caminho); exit(EXIT_FAILURE); }

    // Amostragem evita arquivos enormes e nao modifica a solucao calculada. 
    int amostras = n < 96 ? n : 96;
    for (int a = 0; a <= amostras; a++) {
        int i = a * n / amostras;
        double x = (double)i / n;
        for (int b = 0; b <= amostras; b++) {
            int j = b * n / amostras;
            double y = (double)j / n;
            double analitica = uex(x, y);
            fprintf(dados, "%.17g %.17g %.17g %.17g %.17g\n",
                    x, y, v[i][j], analitica, fabs(v[i][j] - analitica));
        }
        fputc('\n', dados);
    }
    if (fclose(dados) != 0) { perror(caminho); exit(EXIT_FAILURE); }

    FILE *indice = fopen("video_dados/quadros.tsv",
                         video_numero_quadro == 1 ? "w" : "a");
    if (!indice) { perror("video_dados/quadros.tsv"); exit(EXIT_FAILURE); }
    fprintf(indice, "%06d\t%s\t%d\t%d\t%.12e\n",
            video_numero_quadro, etapa, n, ciclo, residuo);
    if (fclose(indice) != 0) { perror("video_dados/quadros.tsv"); exit(EXIT_FAILURE); }

    if (video_numero_quadro == 1) {
        FILE *escala = fopen("video_dados/escala.txt", "w");
        if (!escala) { perror("video_dados/escala.txt"); exit(EXIT_FAILURE); }
        fprintf(escala, "%.17g\n", 1.2 * (double)EXP * EXP);
        if (fclose(escala) != 0) { perror("video_dados/escala.txt"); exit(EXIT_FAILURE); }
    }
}
#endif

int Multigrid(int n, double custom_tol, double *out_h, double *out_erro_inf,
              ResultEntry *accum)
{
    int nx, ny, i, j, ciclo;
    double hx, hy, tol;
    double **u, **f, **r;
    double residuo, erro, erro_inf, tempo;
    struct timespec inicio, fim, tc0, tc1, tg0, tg1;

    bench_reset();

    nx = n; ny = n;
    if (nx <= 0 || ny <= 0) { 
	    printf("Erro: nx,ny>0.\n"); 
	    return 1; }
    if (nx != ny)           { 
	    printf("Erro: nx!=ny.\n"); 
	    return 1; }
    if (nx < 16)            { 
	    printf("Erro: nx>=16.\n"); 
	    return 1; }

    int temp_n = nx;
    while (temp_n > 16) {
        if (temp_n % 2 != 0) { 
		printf("Erro: nx multiplo de 2^n ate 16.\n"); 
		return 1; }
        temp_n /= 2;
    }
    if (temp_n != 16) { 
	    printf("Erro: nx deve ser pot. 2 vezes 16.\n"); 
	    return 1; }

    hx = 1.0 / nx;
    hy = 1.0 / ny;
    if (custom_tol > 0.0) tol = custom_tol; 
    else tol = cfg_tol_global;

    u = aloca_matriz(nx, ny);
    f = aloca_matriz(nx, ny);
    r = aloca_matriz(nx, ny);
    if (!u || !f || !r) { 
	    libera_matriz(u,nx); 
	    libera_matriz(f,nx); 
	    libera_matriz(r,nx); 
	    return 1; }

//inicializacao:

    for (i = 0; i <= nx; i++) 
        for (j = 0; j <= ny; j++) 
	{
            u[i][j] = 0.0; //chute inicial de zeros
            f[i][j] = func(i * hx, j * hy); //funcao em problema.c
        }

    //dirichlet:

    //paredes de baixo e de cima
    for (i = 0; i <= nx; i++) 
    { 
	    u[i][0] = uex(i*hx,0.0); 
	    u[i][ny] = uex(i*hx,1.0);
    }

    //esquerda e direita
    for (j = 0; j <= ny; j++) 
    { 
	    u[0][j] = uex(0.0,j*hy); 
	    u[nx][j] = uex(1.0,j*hy); 
    }

//residuo inicial. funcao em residuo.c
    calcula_residuo(u, f, r, nx, ny, hx);
    residuo = norma_residuo(r, nx, ny, hx, hy);

    clock_gettime(CLOCK_MONOTONIC, &inicio);

    printf("\nMULTIGRID - FMG\n");
    printf("nx = %d  ny = %d\n", nx, ny);
    printf("A = %s  (omega_A = %g)\n", nome_A(cfg_A), cfg_omega_A);

    printf("B = %s", nome_B(cfg_B));
    if (cfg_B >= 4) printf("  (omega_B = %g, tol_B = %.3e)", cfg_omega_B, cfg_tol_B);
    else            printf("  (direto)");
    printf("\n");

    printf("C = %s", nome_C(cfg_C));
    if (cfg_C <= 2) printf("  (omega_C = %g, tol_C = %.3e)", cfg_omega_C, cfg_tol_C);
    else            printf("  (direto)");
    printf("\n");

    printf("K = %g  NU1 = %d  NU2 = %d\n", cfg_K, cfg_NU1, cfg_NU2);
    printf("tolerância global = %.3e\n", tol);
    if (cfg_B >= 4 || cfg_C <= 2)
        printf("máximo de iterações no nível grosso = %d\n", MAX_ITER_GROSSO);
    printf("máximo de ciclos = %d\n", MAX_CICLOS);
    printf("\n");

//niveis será o numero de vezes que precisará engrossar a malha

    int niveis = 1;//comeca em 1
    int n_nivel = nx;
    while (n_nivel > 16) 
    { 
	    n_nivel /= 2; //vai cortando o numero de pontos ate chegar em 16(nivel coarse)
	    niveis++; 
    }

    int *N = malloc(niveis * sizeof(int));//numero de pontos em certo nivel
    if (!N) 
    { 
	    libera_matriz(u,nx); 
	    libera_matriz(f,nx); 
	    libera_matriz(r,nx); 
	    return 1; 
    }

    N[0] = nx;//no nivel 0, nx pontos

    for (i = 1; i < niveis; i++) N[i] = N[i - 1] / 2;//no nivel seguinte, metade do nivel de cima

    double ***f_niveis = malloc(niveis * sizeof(double **));

    if (!f_niveis) 
    { 
	    free(N); 
	    libera_matriz(u,nx); 
	    libera_matriz(f,nx); 
	    libera_matriz(r,nx); 
	    return 1; 
    }

    f_niveis[0] = f;//funcao nos niveis. no zero foi preenchida

    for (int nivel = 1; nivel < niveis; nivel++) 
    {
        int n_fina = N[nivel - 1], n_grossa = N[nivel];
        f_niveis[nivel] = aloca_matriz(n_grossa, n_grossa);
        if (!f_niveis[nivel]) 
	{
            for (int k = 1; k < nivel; k++) libera_matriz(f_niveis[k], N[k]);
	    free(f_niveis); 
	    free(N);
            libera_matriz(u,nx);
	    libera_matriz(f,nx); 
	    libera_matriz(r,nx);
            return 1;
        }
        restringe(f_niveis[nivel - 1], f_niveis[nivel], n_fina, n_fina);//restringe com FW a funcao ao nivel seguinte. funcao em operadores.c
    }

    double ***r_niveis      = malloc(niveis * sizeof(double **));
    double ***e_niveis      = malloc(niveis * sizeof(double **));
    double ***u_novo_niveis = malloc(niveis * sizeof(double **));

    for (int k = 0; k < niveis; k++) 
    {
        r_niveis[k]      = aloca_matriz(N[k], N[k]);
        e_niveis[k]      = aloca_matriz(N[k], N[k]);
        u_novo_niveis[k] = aloca_matriz(N[k], N[k]);

        if (!r_niveis[k] || !e_niveis[k] || !u_novo_niveis[k]) {
            for (int l = 0; l <= k; l++) {
                libera_matriz(r_niveis[l], N[l]);
                libera_matriz(e_niveis[l], N[l]);
                libera_matriz(u_novo_niveis[l], N[l]);
            }

            free(r_niveis); 
	    free(e_niveis); 
	    free(u_novo_niveis);
            
	    for (int l = 1; l < niveis; l++) libera_matriz(f_niveis[l], N[l]);
            free(f_niveis); 
	    free(N);
            libera_matriz(u,nx); 
	    libera_matriz(f,nx); 
	    libera_matriz(r,nx);
            return 1;
        }
    }

    

    int nivel = niveis - 1; 
    int n_atual = N[nivel];
    double h_atual = 1.0 / n_atual;

    double **u_atual = aloca_matriz(n_atual, n_atual);
    double **u_novo  = aloca_matriz(n_atual, n_atual);
    double **r_atual = aloca_matriz(n_atual, n_atual);
    if (!u_atual || !u_novo || !r_atual) {
        libera_matriz(u_atual,n_atual); 
	libera_matriz(u_novo,n_atual); 
	libera_matriz(r_atual,n_atual);
        for (int k = 0; k < niveis; k++) {
            libera_matriz(r_niveis[k], N[k]);
            libera_matriz(e_niveis[k], N[k]);
            libera_matriz(u_novo_niveis[k], N[k]);
            if (k > 0) libera_matriz(f_niveis[k], N[k]);
        }
        free(r_niveis); free(e_niveis); free(u_novo_niveis); free(f_niveis); free(N);
        libera_matriz(u,nx); libera_matriz(f,nx); libera_matriz(r,nx);
        
        return 1;
    }

    zera_matriz(u_atual, n_atual, n_atual);

    for (i = 0; i <= n_atual; i++) 
    {
        u_atual[i][0]       = uex(i*h_atual, 0.0);
        u_atual[i][n_atual] = uex(i*h_atual, 1.0);
    }

    for (j = 0; j <= n_atual; j++) 
    {
        u_atual[0][j]       = uex(0.0, j*h_atual);
        u_atual[n_atual][j] = uex(1.0, j*h_atual);
    }

    // (C) APROX. INICIAL DO FMG 
    clock_gettime(CLOCK_MONOTONIC, &tg0);
    int iter_grosso = 0;
    double residuo_grosso = 0.0;
    int C_iter = (cfg_C == 1 || cfg_C == 2);//var logica. C eh iterativo?

    clock_gettime(CLOCK_MONOTONIC, &tc0);

    if (C_iter)
    {
        calcula_residuo(u_atual, f_niveis[nivel], r_atual, n_atual, n_atual, h_atual);
        residuo_grosso = norma_residuo(r_atual, n_atual, n_atual, h_atual, h_atual);

        
       while (residuo_grosso > cfg_tol_C && iter_grosso < MAX_ITER_GROSSO)
        {
    //C 
            double **res = smooth_step(u_atual, u_novo, f_niveis[nivel],
                               n_atual, n_atual, h_atual, cfg_C, cfg_omega_C);
            if (res == u_novo) 
	    {
		    double **t = u_atual; 
		    u_atual = u_novo; 
		    u_novo = t; 
	    }

            calcula_residuo(u_atual, f_niveis[nivel], r_atual, n_atual, n_atual, h_atual);
            residuo_grosso = norma_residuo(r_atual, n_atual, n_atual, h_atual, h_atual);
#ifdef PRINT
            printf("Resolvendo em %dx%d (C=%s)...resíduo %d: %.12e\n",
                   n_atual, n_atual, nome_C(cfg_C), iter_grosso + 1, residuo_grosso);
#endif
            iter_grosso++;
        }
    }
    else
    {
        //C 
        if      (cfg_C == 3) LU         (u_atual, f_niveis[nivel], n_atual, n_atual, h_atual);
        else if (cfg_C == 4) LDL        (u_atual, f_niveis[nivel], n_atual, n_atual, h_atual);
        else if (cfg_C == 5) gaussjordan(u_atual, f_niveis[nivel], n_atual, n_atual, h_atual);
        printf("(C=%s) resolvido diretamente em %dx%d.\n", nome_C(cfg_C), n_atual, n_atual);
    }

    clock_gettime(CLOCK_MONOTONIC, &tc1);
    bench_t_C += ts_diff(tc0, tc1);
    bench_n_C_iter += iter_grosso;
#ifdef VIDEO
    video_salva_quadro(u_atual, n_atual, "C", 0, residuo_grosso);
#endif


    //agora que resolveu na malha grossa, vai voltar à malha original fazendo uma interpolacao e um vciclo até chegar
    while (nivel > 0)//enquanto eu nao voltei a malha original...
    {
        int n_grosso = N[nivel];//pontos nesse nivel
        int n_fino   = N[nivel - 1];//póntos no de cima
        double h_fino = 1.0 / n_fino;//passo de cima
        double **u_fino = aloca_matriz(n_fino, n_fino);
        if (!u_fino) {
            libera_matriz(u_atual, n_grosso); libera_matriz(u_novo, n_grosso); libera_matriz(r_atual, n_grosso);
            for (int k = 0; k < niveis; k++) {
                libera_matriz(r_niveis[k], N[k]); libera_matriz(e_niveis[k], N[k]);
                libera_matriz(u_novo_niveis[k], N[k]);
                if (k > 0) libera_matriz(f_niveis[k], N[k]);
            }
            free(r_niveis); free(e_niveis); free(u_novo_niveis); free(f_niveis); free(N);
            libera_matriz(u,nx); libera_matriz(f,nx); libera_matriz(r,nx);
            
            return 1;
        }
        prolonga(u_atual, u_fino, n_fino, n_fino);//interpolaçao bilinear para passar a malha de cima (codigo em operadores)
        for (int ib = 0; ib <= n_fino; ib++) {
            u_fino[ib][0] = uex(ib*h_fino, 0.0);
            u_fino[ib][n_fino] = uex(ib*h_fino, 1.0);
            u_fino[0][ib] = uex(0.0, ib*h_fino);
            u_fino[n_fino][ib] = uex(1.0, ib*h_fino);
        }
#ifdef VIDEO
        video_salva_quadro(u_fino, n_fino, "interpolacao", 0, 0.0);
#endif
        v_cycle(u_fino, f_niveis[nivel - 1], n_fino, n_fino, h_fino,
                f_niveis, r_niveis, e_niveis, u_novo_niveis, nivel - 1);//um vcicle (codigo em vcycle)
#ifdef VIDEO
        video_salva_quadro(u_fino, n_fino, "FMG", 0, 0.0);
#endif
        libera_matriz(u_atual, n_grosso);
        libera_matriz(u_novo,  n_grosso);
        u_atual = u_fino;
        u_novo  = NULL;
        n_atual = n_fino;
        nivel--;//(ai subiu um nivel)
    }

    clock_gettime(CLOCK_MONOTONIC, &tg1);
    bench_t_guess += ts_diff(tg0, tg1);

    for (i = 0; i <= nx; i++)
        for (j = 0; j <= ny; j++)
            u[i][j] = u_atual[i][j];

    calcula_residuo(u, f, r, nx, ny, hx);
    residuo = norma_residuo(r, nx, ny, hx, hy);

#ifdef VIDEO
    video_salva_quadro(u, nx, "chute_FMG", 0, residuo);
#endif
    ciclo = 0;

    //ja inicializou. agora é vciclo até convergir
    while (residuo > tol && ciclo < MAX_CICLOS)
    {
        struct timespec tv0, tv1;
        clock_gettime(CLOCK_MONOTONIC, &tv0);
        v_cycle(u, f, nx, ny, hx, f_niveis, r_niveis, e_niveis, u_novo_niveis, 0);
        clock_gettime(CLOCK_MONOTONIC, &tv1);
        bench_t_vcycle += ts_diff(tv0, tv1);
        bench_n_vcycle++;

        calcula_residuo(u, f, r, nx, ny, hx);
        residuo = norma_residuo(r, nx, ny, hx, hy);//o vciclo é recursivo. entao aqui é onde decide se sai do loop ou nao
#ifdef PRINT
        printf("Ciclo %d resíduo: %.12e\n", ciclo, residuo);
#endif
        ciclo++;
#ifdef VIDEO
        // Primeiros ciclos, depois potencias de dois: custo de visualizacao limitado. 
        if (ciclo <= 4 || (ciclo & (ciclo - 1)) == 0)
            video_salva_quadro(u, nx, "ciclo", ciclo, residuo);
#endif
    }
#ifdef VIDEO
    //sempre registrar o estado final, mesmo fora da amostragem de ciclos. 
    video_salva_quadro(u, nx, "final", ciclo, residuo);
#endif

    clock_gettime(CLOCK_MONOTONIC, &fim);
    tempo = ts_diff(inicio, fim);
#ifdef PLOT

    /*
      Salva os parametros efetivamente utilizados pelo programa.
     
      cfg_K pode ter sido alterado pelo usuario durante a execucao,
      portanto nao devemos recuperar K novamente do #define no script.
     */
    FILE *cfgplot = fopen("plot_config.dat", "w");

    if (cfgplot != NULL)
    {
        fprintf(cfgplot, "%.17g %d\n", cfg_K, EXP);
        fclose(cfgplot);
    }

    FILE *gnu = fopen("grafico.dat", "w");

    if (gnu != NULL)
    {
        for (int inteiro = 0; inteiro < nx; inteiro++)
        {
            for (int jota = 0; jota < ny; jota++)
            {
                double x = (double)inteiro / nx;
                double y = (double)jota / ny;

                fprintf(gnu,
                        "%.17g %.17g %.17g %.17g\n",
                        x,
                        y,
                        u[inteiro][jota],
                        fabs(u[inteiro][jota] - uex(x, y)));
            }

            fprintf(gnu, "\n");
        }

        fclose(gnu);
    }

#endif

    erro     = norma_erro(u, nx, ny, hx, hy);
    erro_inf = norma_infinito_erro(u, nx, ny, hx, hy);

    printf("\nRESULTADOS\n");
    printf("A = %s  B = %s  C = %s\n", nome_A(cfg_A), nome_B(cfg_B), nome_C(cfg_C));
    printf("Ciclos V-cycle apos FMG: %d\n", ciclo);
    printf("Iteracoes em C: %d\n", iter_grosso);
    printf("Residuo final: %.12e\n", residuo);
    printf("Norma 2 do erro: %.12e\n", erro);
    printf("Norma infinito do erro: %.12e\n", erro_inf);
    printf("Tempo de execucao: %.6f s\n", tempo);
    if (residuo < tol) printf("Convergencia atingida.\n");
    else               printf("Maximo de ciclos atingido.\n");

    bench_print();
    bench_accumulate_into(accum);

    if (out_h)        *out_h = hx;
    if (out_erro_inf) *out_erro_inf = erro_inf;

    libera_matriz(u_atual, nx);
    libera_matriz(r_atual, 16);
    for (int k = 0; k < niveis; k++) {
        libera_matriz(r_niveis[k], N[k]);
        libera_matriz(e_niveis[k], N[k]);
        libera_matriz(u_novo_niveis[k], N[k]);
        if (k > 0) libera_matriz(f_niveis[k], N[k]);
    }
    free(r_niveis); free(e_niveis); free(u_novo_niveis); free(f_niveis); free(N);
    libera_matriz(u,nx); libera_matriz(f,nx); libera_matriz(r,nx);
    
    return 0;
}

