#define _USE_MATH_DEFINES

#define K 2
#define OMEGA (4.0/5.0)
#define OMEGAGS 1.0
#define EXP 8
#define TOLERANCIA 1e-8
#define MAX_ITER 1000000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/*
    Problema:

        -Delta u = f,     em Omega = (0,1)x(0,1)
               u = 0,     em dOmega

    Solucao exata:

        u(x,y) = sin(K*pi*x) sin(K*pi*y)

    Logo:

        f(x,y) = 2 K^2 pi^2 sin(K*pi*x) sin(K*pi*y)

    Discretizacao:

        (4u_ij - u_(i-1,j) - u_(i+1,j)
               - u_(i,j-1) - u_(i,j+1))/h^2 = f_ij
*/




double **aloca_matriz(int n)
{
    double **A = malloc((n + 1) * sizeof(double *));

    if (A == NULL)
        return NULL;

    for (int i = 0; i <= n; i++)
    {
        A[i] = malloc((n + 1) * sizeof(double));

        if (A[i] == NULL)
        {
            for (int j = 0; j < i; j++)
                free(A[j]);

            free(A);
            return NULL;
        }
    }

    return A;
}


void libera_matriz(double **A, int n)
{
    if (A == NULL)
        return;

    for (int i = 0; i <= n; i++)
        free(A[i]);

    free(A);
}


double soma (double x)
{
	double s=0.0;
	int freq = K;

	for (int j=1; j<EXP; j++)
	{
		s += sin(freq*M_PI*x);
		freq*=K;
	}
	return s;
}

double soma2 (double x)
{
	double s=0.0;
	int freq =K;

	for (int i=1; i<EXP; i++)
	{
		s+= (double)freq*freq*sin(freq*M_PI*x);
		freq *=K;
	}
	return s;
}

double uex(double x, double y)
{
    return (soma(x) * soma(y));
}

double func(double x, double y)
{
    return M_PI*M_PI*(soma2(x)*soma(y)+soma(x)*soma2(y));
}




void inicializa(double **u, double **f, int n, double h)
{
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            double x = i * h;
            double y = j * h;

            u[i][j] = 0.0;
            f[i][j] = func(x, y);
        }
    }

    /*
        Condicoes de Dirichlet.

        Neste problema sao zero, mas escrevemos usando uex()
        para deixar explicito de onde vem o contorno.
    */

    for (int i = 0; i <= n; i++)
    {
        u[i][0] = uex(i*h, 0.0);
        u[i][n] = uex(i*h, 1.0);
    }

    for (int j = 0; j <= n; j++)
    {
        u[0][j] = uex(0.0, j*h);
        u[n][j] = uex(1.0, j*h);
    }
}




double norma_residuo(double **u, double **f, int n, double h)
{
    double soma = 0.0;

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            double Au =
                (4.0 * u[i][j]
                - u[i-1][j]
                - u[i+1][j]
                - u[i][j-1]
                - u[i][j+1]) / (h*h);

            double r = f[i][j] - Au;

            soma += r * r;
        }
    }

    return sqrt(h * h * soma);
}




double norma_erro(double **u, int n, double h)
{
    double soma = 0.0;

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            double erro =
                u[i][j] - uex(i*h, j*h);

            soma += erro * erro;
        }
    }

    return sqrt(h * h * soma);
}


double norma_infinito_erro(double **u, int n, double h)
{
    double maximo = 0.0;

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            double erro =
                fabs(u[i][j] - uex(i*h, j*h));

            if (erro > maximo)
                maximo = erro;
        }
    }

    return maximo;
}




int jacobi_ponderado(double **u, double **f, int n, double h)
{
    double **u_novo = aloca_matriz(n);

    if (u_novo == NULL)
    {
        printf("Erro de alocacao.\n");
        exit(1);
    }

    /*
        Copia inicialmente a matriz inteira para que
        o contorno fique definido em u_novo.
    */

    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            u_novo[i][j] = u[i][j];


    int iter = 0;
    double residuo = norma_residuo(u, f, n, h);


    while (residuo > TOLERANCIA && iter < MAX_ITER)
    {
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < n; j++)
            {
                double u_jacobi =
                    0.25 *
                    (h*h*f[i][j]
                    + u[i-1][j]
                    + u[i+1][j]
                    + u[i][j-1]
                    + u[i][j+1]);

                u_novo[i][j] =
                    (1.0 - OMEGA) * u[i][j]
                    + OMEGA * u_jacobi;
            }
        }


        /*
            O estado novo passa a ser o estado corrente.
        */

        for (int i = 1; i < n; i++)
            for (int j = 1; j < n; j++)
                u[i][j] = u_novo[i][j];


        iter++;

        residuo = norma_residuo(u, f, n, h);
	printf("Iteração %d resíduo: %.12e\n", iter, residuo);
    }


    libera_matriz(u_novo, n);

    return iter;
}




int gauss_seidel_ponderado(double **u, double **f, int n, double h)
{
    int iter = 0;

    double residuo =
        norma_residuo(u, f, n, h);


    while (residuo > TOLERANCIA && iter < MAX_ITER)
    {
        /*
            Gauss-Seidel escreve diretamente em u.

            Durante a varredura, os valores u[i-1][j]
            e u[i][j-1] ja pertencem a iteracao nova.
        */

        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < n; j++)
            {
                double u_gs =
                    0.25 *
                    (h*h*f[i][j]
                    + u[i-1][j]
                    + u[i+1][j]
                    + u[i][j-1]
                    + u[i][j+1]);


                /*
                    Gauss-Seidel ponderado:

                    u_novo =
                        (1-omega) u_velho
                        + omega u_GS

                    omega = 1:
                        Gauss-Seidel convencional.

                    0 < omega < 1:
                        sub-relaxacao.

                    1 < omega < 2:
                        sobre-relaxacao (SOR).
                */

                u[i][j] =
                    (1.0 - OMEGAGS) * u[i][j]
                    + OMEGAGS * u_gs;
            }
        }


        iter++;

        residuo =
            norma_residuo(u, f, n, h);
	printf("Iteracao: %d, resíduo: %.12e\n", iter, residuo);
    }


    return iter;
}




int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("Uso: %s N\n", argv[0]);
        printf("Exemplo: %s 64\n", argv[0]);
        return 1;
    }


    int n = atoi(argv[1]);

    if (n < 2)
    {
        printf("Erro: N deve ser maior ou igual a 2.\n");
        return 1;
    }


    double h = 1.0 / n;


    double **u = aloca_matriz(n);
    double **f = aloca_matriz(n);


    if (u == NULL || f == NULL)
    {
        printf("Erro de alocacao.\n");

        libera_matriz(u, n);
        libera_matriz(f, n);

        return 1;
    }


    inicializa(u, f, n, h);


    struct timespec inicio, fim;

    clock_gettime(CLOCK_MONOTONIC, &inicio);


    


    int iteracoes = jacobi_ponderado(u, f, n, h);

    // int iteracoes = gauss_seidel_ponderado(u, f, n, h); 


    clock_gettime(CLOCK_MONOTONIC, &fim);


    double tempo =
        (fim.tv_sec - inicio.tv_sec)
        + (fim.tv_nsec - inicio.tv_nsec) * 1e-9;


    double residuo =
        norma_residuo(u, f, n, h);

    double erro2 =
        norma_erro(u, n, h);

    double erro_inf =
        norma_infinito_erro(u, n, h);


    
    printf("\nRESULTADOS\n");
    

    printf("N                = %d\n", n);
    printf("h                = %.12e\n", h);
    printf("K                = %.6f\n", K);
    printf("omega            = %.6f\n", OMEGA);

    printf("iteracoes        = %d\n", iteracoes);

    printf("residuo final    = %.12e\n", residuo);
    printf("norma 2 do erro  = %.12e\n", erro2);
    printf("norma inf. erro  = %.12e\n", erro_inf);

    printf("tempo            = %.6f s\n", tempo);

    if (residuo <= TOLERANCIA)
        printf("convergencia     = atingida\n");
    else
        printf("convergencia     = NAO atingida\n");

    


    libera_matriz(u, n);
    libera_matriz(f, n);


    return 0;
}
